// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"path/filepath"
	"time"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/msolid"

	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
)

// Global holds global data
var Global struct {

	// constants
	LogPrefix string // extra string to prefix log file

	// multiprocessing data
	Rank     int   // my rank in distributed cluster
	Nproc    int   // number of processors
	Root     bool  // am I root? (i.e. myrank == 0)
	Distr    bool  // distributed simulation with more than one mpi processor
	Verbose  bool  // verbose == root
	WspcStop []int // stop flags [nprocs]
	WspcInum []int // workspace of integer numbers [nprocs]

	// simulation, materials, meshes and convenience variables
	Sim    *inp.Simulation // simulation data
	Ndim   int             // space dimension
	Dirout string          // directory for output of results
	Fnkey  string          // filename key; e.g. mysim.sim => mysim
	Enc    string          // encoder; e.g. "gob" or "json"
	Stat   bool            // save residuals in summary
	LogBcs bool            // log essential and ptnatural boundary conditions
	Debug  bool            // debug flag

	// auxiliar structures
	DynCoefs *DynCoefs    // dynamic coefficients
	HydroSt  *HydroStatic // computes hydrostatic states

	// time control
	Time float64 // curent simulation time

	// domains and summary
	Domains []*Domain // all domains
	Summary *Summary  // summary structure
	Solver  FEsolver  // finite element method solver; e.g. implicit, Richardson extrapolation, etc.

	// for debugging
	DebugKb func(d *Domain, it int) // debug Kb callback function
}

// Start initialises Global data and starts logging
//  Input:
//    simfilepath   -- .sim filename including full path
//    erasefiles    -- do erase previous results files
//    allowParallel -- allow parallel execution; otherwise,
//                     run in serial mode regardless whether MPI is on or not
//  Output:
//    startisok -- success
func Start(simfilepath string, erasefiles, verbose, allowParallel bool) (startisok bool) {

	// multiprocessing data
	Global.Rank = 0
	Global.Nproc = 1
	Global.Root = true
	Global.Distr = false
	if mpi.IsOn() && allowParallel {
		Global.Rank = mpi.Rank()
		Global.Nproc = mpi.Size()
		Global.Root = Global.Rank == 0
		Global.Distr = Global.Nproc > 1
	}
	Global.Verbose = verbose
	if !Global.Root {
		Global.Verbose = false
	}
	Global.WspcStop = make([]int, Global.Nproc)
	Global.WspcInum = make([]int, Global.Nproc)

	// simulation and convenience variables
	dir := filepath.Dir(simfilepath)
	fn := filepath.Base(simfilepath)
	Global.Sim = inp.ReadSim(dir, fn, Global.LogPrefix, erasefiles)
	LogErrCond(Global.Sim == nil, "ReadSim failed\n")
	if Stop() {
		return
	}
	Global.Ndim = Global.Sim.Ndim
	Global.Dirout = Global.Sim.Data.DirOut
	Global.Fnkey = Global.Sim.Data.FnameKey
	Global.Enc = Global.Sim.Data.Encoder
	Global.Stat = Global.Sim.Data.Stat
	Global.LogBcs = Global.Sim.Data.LogBcs
	Global.Debug = Global.Sim.Data.Debug

	// linear solver name
	if mpi.IsOn() && allowParallel {
		if mpi.Size() > 1 {
			Global.Sim.LinSol.Name = "mumps"
		}
	} else {
		Global.Sim.LinSol.Name = "umfpack"
	}

	// fix show residual flag
	if !Global.Root {
		Global.Sim.Data.ShowR = false
	}

	// auxiliar structures
	Global.DynCoefs = new(DynCoefs)
	if !Global.DynCoefs.Init(&Global.Sim.Solver) {
		return
	}
	Global.HydroSt = new(HydroStatic)
	Global.HydroSt.Init()

	// plot functions
	if Global.Sim.PlotF != nil {
		if Global.Root {
			Global.Sim.Functions.PlotAll(Global.Sim.PlotF, Global.Dirout, Global.Fnkey)
		}
		if Global.Verbose {
			io.Pfyel("\nfunctions plotted\n")
		}
		return
	}

	// domains and summary
	Global.Domains = make([]*Domain, 0)
	Global.Summary = nil
	Global.Solver = nil

	// success
	return true
}

// RunAll runs FE simulation, after allocating domains and summary. All active stages are run.
func RunAll() (runisok bool) {

	// clean up
	defer CleanUp()

	// benchmarking
	cputime := time.Now()
	defer func() {
		if Global.Verbose && !Global.Debug {
			io.Pf("\nfinal time = %v\n", Global.Time)
			io.Pfblue2("cpu time   = %v\n", time.Now().Sub(cputime))
		}
	}()

	// allocate domains, summary, solver
	if !Alloc(true) {
		return
	}

	// solve for all stages
	if !SolveAllStages(true) {
		return
	}

	// success
	return true
}

// Alloc allocates domains, summary and solver structures
func Alloc(withSummary bool) (allocisok bool) {

	// alloc domains
	for _, reg := range Global.Sim.Regions {
		dom := NewDomain(reg, Global.Distr)
		if dom == nil {
			break
		}
		Global.Domains = append(Global.Domains, dom)
	}
	if Stop() {
		return
	}

	// alloc summary
	if withSummary {
		Global.Summary = new(Summary)
	}

	// alloc solver
	if alloc, ok := solverallocators[Global.Sim.Solver.Type]; ok {
		Global.Solver = alloc()
	} else {
		LogErrCond(true, "cannot find solver type=%q. e.g. {imp, exp, rex} => implicit, explicit, Richardson extrapolation", Global.Sim.Solver.Type)
		return
	}

	// success
	return true
}

// SolveAllStages solve problem for all stages
//  Input:
//   output -- perform output of results and log material models
func SolveAllStages(output bool) (ok bool) {

	// current time
	Global.Time = 0.0

	// loop over stages
	for stgidx, stg := range Global.Sim.Stages {

		// set stage
		if !SetStage(stgidx) {
			return
		}

		// initialise solution vectors
		if !InitSolution(stgidx, false) {
			return
		}

		// log models
		if output {
			mconduct.LogModels()
			mreten.LogModels()
			mporous.LogModels()
			msolid.LogModels()
		}

		// skip stage?
		if stg.Skip {
			continue
		}

		// time loop
		if !Global.Solver.Run(stg) {
			return
		}
	}
	return true
}

// SolveOneStage solves one stage that was already set
//  Input:
//   stgidx  -- stage index (in Global.Sim.Stages)
//   zeroSol -- zero vectors in domains.Sol
func SolveOneStage(stgidx int, zeroSol bool) (ok bool) {

	// current time and stage
	Global.Time = 0.0
	stg := Global.Sim.Stages[stgidx]

	// initialise solution vectors
	if !InitSolution(stgidx, zeroSol) {
		return
	}

	// time loop
	if !Global.Solver.Run(stg) {
		return
	}

	// success
	return true
}

// SetStage sets stage for all domains
//  Input:
//   stgidx -- stage index (in Global.Sim.Stages)
func SetStage(stgidx int) (ok bool) {
	for _, d := range Global.Domains {
		if LogErrCond(!d.SetStage(stgidx, Global.Sim.Stages[stgidx], Global.Distr), "SetStage failed") {
			break
		}
		d.Sol.T = Global.Time
	}
	if Stop() {
		return
	}
	return true
}

// InitSolution initialises solution vectors (Y, dYdt, internal values such as States.Sig, etc.)
// in all domains for all nodes and all elements
//  Input:
//   stgidx  -- stage index (in Global.Sim.Stages)
//   zeroSol -- zero vectors in domains.Sol
func InitSolution(stgidx int, zeroSol bool) (ok bool) {
	for _, d := range Global.Domains {
		if LogErrCond(!d.SetIniVals(Global.Sim.Stages[stgidx], zeroSol), "InitSolution failed") {
			break
		}
	}
	if Stop() {
		return
	}
	return true
}

// CleanUp cleans memory and flush log
func CleanUp() {

	// flush log
	inp.FlushLog()

	// domains: clean memory
	for _, d := range Global.Domains {
		d.End()
	}

	// summary: save file
	if Global.Summary != nil {
		Global.Summary.Save()
	}
}

// AllocSetAndInit allocates domains, summary, solver and sets stage # stgidx
// and initial values in all domains. It also returns the first domain.
//  Input:
//   stgidx      -- stage index
//   withSummary -- also allocate summary
//   readSummary -- reads summary back from previous calculation
//  Output:
//   dom -- first domain; i.e. dom := Global.Domains[0];
//   sum -- summary, if withSummary is true
func AllocSetAndInit(stgidx int, withSummary, readSummary bool) (dom *Domain, sum *Summary, ok bool) {

	// allocate domain and others
	if !Alloc(withSummary) {
		return
	}

	// set stage
	if !SetStage(stgidx) {
		return
	}

	// set initial solution vectors
	if !InitSolution(stgidx, false) {
		return
	}

	// read summary
	if readSummary {
		if !Global.Summary.Read(Global.Dirout, Global.Fnkey) {
			return
		}
	}

	// success
	return Global.Domains[0], Global.Summary, true
}

// factory ////////////////////////////////////////////////////////////////////////////////////////

// FEsolver implements the actual solver (time loop)
type FEsolver interface {
	Run(stg *inp.Stage) (ok bool)
}

// solverallocators holds all available solvers
var solverallocators = make(map[string]func() FEsolver)

// auxiliary //////////////////////////////////////////////////////////////////////////////////////

// debug_print_p_results print results
func debug_print_p_results(d *Domain) {
	io.Pf("\ntime = %23.10f\n", d.Sol.T)
	for _, v := range d.Msh.Verts {
		n := d.Vid2node[v.Id]
		eqpl := n.GetEq("pl")
		var pl float64
		if eqpl >= 0 {
			pl = d.Sol.Y[eqpl]
		}
		if math.Abs(pl) < 1e-13 {
			pl = 0
		}
		io.Pf("%3d : pl=%23.10v\n", v.Id, pl)
	}
}

// debug_print_up_results print results
func debug_print_up_results(d *Domain) {
	io.Pf("\ntime = %23.10f\n", d.Sol.T)
	for _, v := range d.Msh.Verts {
		n := d.Vid2node[v.Id]
		eqpl := n.GetEq("pl")
		equx := n.GetEq("ux")
		equy := n.GetEq("uy")
		var pl, ux, uy float64
		if eqpl >= 0 {
			pl = d.Sol.Y[eqpl]
		}
		if equx >= 0 {
			ux = d.Sol.Y[equx]
		}
		if equy >= 0 {
			uy = d.Sol.Y[equy]
		}
		if math.Abs(pl) < 1e-13 {
			pl = 0
		}
		if math.Abs(ux) < 1e-13 {
			ux = 0
		}
		if math.Abs(uy) < 1e-13 {
			uy = 0
		}
		io.Pf("%3d : pl=%23.10v ux=%23.10f uy=%23.10f\n", v.Id, pl, ux, uy)
	}
}

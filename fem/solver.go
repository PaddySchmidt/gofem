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

	// time control
	Time    float64 // curent simulation time
	TimeOut float64 // time for output
	TimeIdx int     // time output index

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

	// for debugging
	DebugKb func(d *Domain, it int) // debug Kb callback function
}

// End must be called and the end to flush log file
func End() {
	inp.FlushLog()
}

// Start initialises 'global' and starts logging
func Start(simfilepath string, erasefiles, verbose bool) (startisok bool) {

	// multiprocessing data
	Global.Rank = 0
	Global.Nproc = 1
	Global.Root = true
	Global.Distr = false
	if mpi.IsOn() {
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

	// success
	return true
}

// Run runs FE simulation
func Run() (runisok bool) {

	// current time and output time
	Global.Time = 0.0
	Global.TimeOut = 0.0
	Global.TimeIdx = 0

	// plot functions
	if Global.Sim.PlotF != nil {
		if Global.Root {
			Global.Sim.Functions.PlotAll(Global.Sim.PlotF, Global.Dirout, Global.Fnkey)
		}
		io.PfRed("\nfunctions plotted => simulation not ran\n")
		return
	}

	// alloc domains
	var domains []*Domain
	for _, reg := range Global.Sim.Regions {
		dom := NewDomain(reg, Global.Distr)
		if dom == nil {
			break
		}
		domains = append(domains, dom)
	}
	if Stop() {
		return
	}

	// make sure to call linear solver clean up routines upon exit
	defer func() {
		for _, d := range domains {
			if !d.InitLSol {
				d.LinSol.Clean()
			}
		}
	}()

	// summary of outputs; e.g. with output times
	cputime := time.Now()
	var sum Summary
	sum.OutTimes = []float64{Global.Time}
	defer func() {
		sum.Save()
		if Global.Verbose && !Global.Debug {
			io.Pf("\nfinal time = %v\n", Global.Time)
			io.Pfblue2("cpu time   = %v\n", time.Now().Sub(cputime))
		}
	}()

	// alloc solver
	var solver FEsolver
	if alloc, ok := solverallocators[Global.Sim.Solver.Type]; ok {
		solver = alloc(domains, &sum)
	} else {
		LogErrCond(true, "cannot find solver type=%q. e.g. {imp, exp, rex} => implicit, explicit, Richardson extrapolation", Global.Sim.Solver.Type)
		return
	}

	// loop over stages
	for stgidx, stg := range Global.Sim.Stages {

		// time for output
		Global.TimeOut = Global.Time + stg.Control.DtoFunc.F(Global.Time, nil)

		// set stage
		for _, d := range domains {
			if LogErrCond(!d.SetStage(stgidx, Global.Sim.Stages[stgidx], Global.Distr), "SetStage failed") {
				break
			}
			d.Sol.T = Global.Time
			if !d.Out(Global.TimeIdx) {
				break
			}
		}
		if Stop() {
			return
		}
		Global.TimeIdx += 1

		// log models
		mconduct.LogModels()
		mreten.LogModels()
		mporous.LogModels()
		msolid.LogModels()

		// skip stage?
		if stg.Skip {
			continue
		}

		// time loop
		if !solver.Run(stg) {
			return
		}
	}
	return true
}

// factory ////////////////////////////////////////////////////////////////////////////////////////

// FEsolver implements the actual solver (time loop)
type FEsolver interface {
	Run(stg *inp.Stage) (ok bool)
}

// solverallocators holds all available solvers
var solverallocators = make(map[string]func([]*Domain, *Summary) FEsolver)

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

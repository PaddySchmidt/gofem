// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package fem contains elements and solvers for running simulations using the finite element method
package fem

import (
	"math"
	"time"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
)

// function to debug global Jacobian matrix
type DebugKb_t func(d *Domain, it int)

// FEsolver implements the actual solver (time loop)
type FEsolver interface {
	Run(tf float64, dtFunc, dtoFunc fun.Func, verbose bool, dbgKb DebugKb_t) (err error)
}

// solverallocators holds all available solvers
var solverallocators = make(map[string]func(doms []*Domain, sum *Summary, dc *DynCoefs) FEsolver)

// FEM holds all data for a simulation using the finite element method
type FEM struct {
	Sim     *inp.Simulation // simulation data
	Summary *Summary        // summary structure
	DynCfs  *DynCoefs       // coefficients for dynamics/transient simulations
	HydSta  *HydroStatic    // function to compute hydrostatic state
	Domains []*Domain       // all domains
	Solver  FEsolver        // finite element method solver; e.g. implicit, Richardson extrapolation, etc.
	DebugKb DebugKb_t       // debug Kb callback function
	Proc    int             // processor id
	Verbose bool            // show messages
}

// Clean cleans memory allocated by FEM
func (o *FEM) Clean() {
	for _, d := range o.Domains {
		d.Clean()
	}
}

// NewFEM returns a new FEM structure
//  Input:
//   simfilepath   -- simulation (.sim) filename including full path
//   alias         -- word to be appended to simulation key; e.g. when running multiple FE solutions
//   erasePrev     -- erase previous results files
//   readSummary   -- ready summary of previous simulation
//   allowParallel -- allow parallel execution; otherwise, run in serial mode regardless whether MPI is on or not
//   verbose       -- show messages
func NewFEM(simfilepath, alias string, erasePrev, readSummary, allowParallel, verbose bool) (o *FEM) {

	// new FEM object
	o = new(FEM)

	// read input data
	o.Sim = inp.ReadSim(simfilepath, alias, erasePrev)
	if o.Sim == nil {
		chk.Panic("cannot ready simulation input data")
	}

	// read summary of previous simulation
	if readSummary {
		o.Summary = new(Summary)
		err := o.Summary.Read(o.Sim.DirOut, o.Sim.Key, o.Sim.EncType)
		if err != nil {
			chk.Panic("cannot ready summary:\n%v", err)
		}
	}

	// multiprocessing data
	nproc, distr := 1, false
	if mpi.IsOn() {
		if allowParallel {
			o.Proc = mpi.Rank()
			nproc = mpi.Size()
			distr = nproc > 1
			if nproc > 1 {
				o.Sim.LinSol.Name = "mumps"
			}
		}
	} else {
		o.Sim.LinSol.Name = "umfpack"
	}
	o.Verbose = verbose && (o.Proc == 0)

	// auxiliary structures
	o.DynCfs = new(DynCoefs)
	o.DynCfs.Init(&o.Sim.Solver)
	o.HydSta = new(HydroStatic)
	o.HydSta.Init(o.Sim.WaterLevel, o.Sim.WaterRho0, o.Sim.WaterBulk, o.Sim.Gravity.F(0, nil))

	// allocate domains
	o.Domains = NewDomains(o.Sim, o.DynCfs, o.HydSta, o.Proc, nproc, distr)

	// allocate solver
	if alloc, ok := solverallocators[o.Sim.Solver.Type]; ok {
		o.Solver = alloc(o.Domains, o.Summary, o.DynCfs)
	} else {
		chk.Panic("cannot find solver type named %q", o.Sim.Solver.Type)
	}
	return
}

// Run runs FE simulation
func (o *FEM) Run() (err error) {

	// plot functions
	if o.Sim.PlotF != nil {
		if o.Proc == 0 {
			o.Sim.Functions.PlotAll(o.Sim.PlotF, o.Sim.DirOut, o.Sim.Key)
		}
		if o.Verbose {
			io.Pfyel("\nfunctions plotted\n")
		}
		return
	}

	// benchmarking
	cputime := time.Now()
	defer func() {
		if o.Verbose {
			io.Pf("\nfinal time = %v\n", o.Domains[0].Sol.T)
			io.Pfblue2("cpu time   = %v\n", time.Now().Sub(cputime))
		}
	}()

	// loop over stages
	for stgidx, stg := range o.Sim.Stages {

		// skip stage?
		if stg.Skip {
			continue
		}

		// set stage
		err = o.SetStage(stgidx)
		if err != nil {
			return
		}

		// initialise solution vectors
		err = o.ZeroStage(stgidx, true)
		if err != nil {
			return
		}

		// time loop
		err = o.Solver.Run(stg.Control.Tf, stg.Control.DtFunc, stg.Control.DtoFunc, o.Verbose, o.DebugKb)
		if err != nil {
			return
		}
	}
	return
}

// SetStage sets stage for all domains
//  Input:
//   stgidx -- stage index (in o.Sim.Stages)
func (o *FEM) SetStage(stgidx int) (err error) {
	for _, d := range o.Domains {
		err = d.SetStage(stgidx)
		if err != nil {
			return
		}
	}
	return
}

// ZeroStage zeroes solution varaibles; i.e. it initialises solution vectors (Y, dYdt, internal
// values such as States.Sig, etc.) in all domains for all nodes and all elements
//  Input:
//   stgidx  -- stage index (in o.Sim.Stages)
//   zeroSol -- zero vectors in domains.Sol
func (o *FEM) ZeroStage(stgidx int, zeroSol bool) (err error) {
	for _, d := range o.Domains {
		err = d.SetIniVals(stgidx, zeroSol)
		if err != nil {
			return
		}
	}
	return
}

/*
// SolveOneStage solves one stage that was already set
//  Input:
//   stgidx  -- stage index (in o.Sim.Stages)
//   zeroSol -- zero vectors in domains.Sol
func (o *FEM) SolveOneStage(stgidx int, zeroSol bool) (ok bool) {

	// current time and stage
	o.Time = 0.0
	stg := o.Sim.Stages[stgidx]

	// initialise solution vectors
	if !o.InitSolution(stgidx, zeroSol) {
		return
	}

	// time loop
	if !o.Solver.Run(stg) {
		return
	}

	// success
	return true
}


// CleanUp cleans memory and flush log
func (o *FEM) CleanUp() {

	// flush log
	io.WriteFile(io.Sf("%s/%s_p%d.log", o.Data.DirOut, o.Data.LogPrefix+o.Data.FnameKey, o.Rank))

	// domains: clear memory
	for _, d := range o.Domains {
		d.End()
	}

	// summary: save file
	if o.Summary != nil {
		o.Summary.Save()
	}
}

// AllocSetAndInit allocates domains, summary, solver and sets stage # stgidx
// and initial values in all domains. It also returns the first domain.
//  Input:
//   stgidx      -- stage index
//   withSummary -- also allocate summary
//   readSummary -- reads summary back from previous calculation
//  Output:
//   dom -- first domain; i.e. dom := o.Domains[0];
//   sum -- summary, if withSummary is true
func (o *FEM) AllocSetAndInit(stgidx int, withSummary, readSummary bool) (dom *Domain, sum *Summary, ok bool) {

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
	}

	// success
	return o.Domains[0], o.Summary, true
}

*/

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

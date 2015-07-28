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
	Nproc   int             // number of processors
	Proc    int             // processor id
	Verbose bool            // show messages
}

// NewFEM returns a new FEM structure
//  Input:
//   simfilepath   -- simulation (.sim) filename including full path
//   alias         -- word to be appended to simulation key; e.g. when running multiple FE solutions
//   erasePrev     -- erase previous results files
//   saveSummary   -- save summary
//   readSummary   -- ready summary of previous simulation
//   allowParallel -- allow parallel execution; otherwise, run in serial mode regardless whether MPI is on or not
//   verbose       -- show messages
func NewFEM(simfilepath, alias string, erasePrev, saveSummary, readSummary, allowParallel, verbose bool) (o *FEM) {

	// new FEM object
	o = new(FEM)

	// read input data
	o.Sim = inp.ReadSim(simfilepath, alias, erasePrev)
	if o.Sim == nil {
		chk.Panic("cannot ready simulation input data")
	}

	// read summary of previous simulation
	if saveSummary || readSummary {
		o.Summary = new(Summary)
	}
	if readSummary {
		err := o.Summary.Read(o.Sim.DirOut, o.Sim.Key, o.Sim.EncType)
		if err != nil {
			chk.Panic("cannot ready summary:\n%v", err)
		}
	}

	// multiprocessing data
	o.Nproc = 1
	distr := false
	if mpi.IsOn() {
		if allowParallel {
			o.Proc = mpi.Rank()
			o.Nproc = mpi.Size()
			distr = o.Nproc > 1
			if distr {
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
	o.Domains = NewDomains(o.Sim, o.DynCfs, o.HydSta, o.Proc, o.Nproc, distr)

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

	// loop over stages
	cputime := time.Now()
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

	// message
	if o.Verbose {
		io.Pf("\n\n")
		if len(o.Domains) > 0 {
			if o.Domains[0].Sol != nil {
				io.Pf("\nfinal time = %v\n", o.Domains[0].Sol.T)
			}
		}
		io.Pflmag("cpu time   = %v\n", time.Now().Sub(cputime))
	}

	// save summary
	if o.Summary != nil {
		err = o.Summary.Save(o.Sim.DirOut, o.Sim.Key, o.Sim.EncType, o.Nproc, o.Proc, o.Verbose)
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

// SolveOneStage solves one stage that was already set
//  Input:
//   stgidx    -- stage index (in o.Sim.Stages)
//   zerostage -- zero vectors in domains.Sol => call ZeroStage
func (o *FEM) SolveOneStage(stgidx int, zerostage bool) (err error) {

	// clean memory allocated by domains
	defer func() {
		for _, d := range o.Domains {
			d.Clean()
		}
	}()

	// zero stage
	if zerostage {
		err = o.ZeroStage(stgidx, true)
		if err != nil {
			return
		}
	}

	// run
	stg := o.Sim.Stages[stgidx]
	err = o.Solver.Run(stg.Control.Tf, stg.Control.DtFunc, stg.Control.DtoFunc, o.Verbose, o.DebugKb)
	return
}

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

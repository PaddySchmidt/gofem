// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/mpi"
)

// SolverImplicit solves FEM problem using an implicit procedure (with Newthon-Raphson method)
type SolverImplicit struct {
}

// set factory
func init() {
	solverallocators["imp"] = func() FEsolver {
		solver := new(SolverImplicit)
		return solver
	}
}

func (o *SolverImplicit) Run(stg *inp.Stage) (runisok bool) {

	// auxiliary
	md := 1.0    // time step multiplier if divergence control is on
	ndiverg := 0 // number of steps diverging

	// time control
	t := Global.Time
	tf := stg.Control.Tf
	tout := Global.TimeOut
	tidx := Global.TimeIdx
	defer func() {
		Global.Time = t
		Global.TimeOut = tout
		Global.TimeIdx = tidx
	}()

	// time loop
	var Δt, Δtout float64
	var lasttimestep bool
	for t < tf {

		// check for continued divergence
		if LogErrCond(ndiverg >= Global.Sim.Solver.NdvgMax, "continuous divergence after %d steps reached", ndiverg) {
			return
		}

		// time increment
		Δt = stg.Control.DtFunc.F(t, nil) * md
		if t+Δt >= tf {
			Δt = tf - t
			lasttimestep = true
		}
		if Δt < Global.Sim.Solver.DtMin {
			if md < 1 {
				LogErrCond(true, "Δt increment is too small: %g < %g", Δt, Global.Sim.Solver.DtMin)
				return false
			}
			return true
		}

		// dynamic coefficients
		if LogErr(Global.DynCoefs.CalcBoth(Δt), "cannot compute dynamic coefficients") {
			return
		}

		// time update
		t += Δt
		for _, d := range Global.Domains {
			d.Sol.T = t
		}
		Δtout = stg.Control.DtoFunc.F(t, nil)

		// message
		if Global.Verbose {
			if !Global.Sim.Data.ShowR && !Global.Debug {
				io.PfWhite("%30.15f\r", t)
			}
		}

		// for all domains
		docontinue := false
		for _, d := range Global.Domains {

			// backup solution if divergence control is on
			if Global.Sim.Solver.DvgCtrl {
				d.backup()
			}

			// run iterations
			diverging, ok := run_iterations(t, Δt, d)
			if !ok {
				return
			}

			// restore solution and reduce time step if divergence control is on
			if Global.Sim.Solver.DvgCtrl {
				if diverging {
					if Global.Verbose {
						io.Pfred(". . . iterations diverging (%2d) . . .\n", ndiverg+1)
					}
					d.restore()
					t -= Δt
					d.Sol.T = t
					md *= 0.5
					ndiverg += 1
					docontinue = true
					break
				}
				ndiverg = 0
				md = 1.0
			}
		}
		if docontinue {
			continue
		}

		// perform output
		if t >= tout || lasttimestep {
			if Global.Summary != nil {
				Global.Summary.OutTimes = append(Global.Summary.OutTimes, t)
			}
			for _, d := range Global.Domains {
				//if true {
				if false {
					debug_print_p_results(d)
				}
				if false {
					debug_print_up_results(d)
				}
				if !d.Out(tidx) {
					break
				}
			}
			if Stop() {
				return
			}
			tout += Δtout
			tidx += 1
		}
	}

	// success
	return true
}

// run_iterations solves the nonlinear problem
func run_iterations(t, Δt float64, d *Domain) (diverging, ok bool) {

	// zero accumulated increments
	la.VecFill(d.Sol.ΔY, 0)

	// calculate global starred vectors and interpolate starred variables from nodes to integration points
	if LogErr(d.star_vars(Δt), "cannot compute starred variables") {
		return
	}

	// auxiliary variables
	var it int
	var largFb, largFb0, Lδu float64
	var prevFb, prevLδu float64

	// message
	if Global.Sim.Data.ShowR {
		io.Pf("\n%13s%4s%23s%23s\n", "t", "it", "largFb", "Lδu")
		defer func() {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}()
	}

	// iterations
	for it = 0; it < Global.Sim.Solver.NmaxIt; it++ {

		// assemble right-hand side vector (fb) with negative of residuals
		la.VecFill(d.Fb, 0)
		for _, e := range d.Elems {
			if !e.AddToRhs(d.Fb, d.Sol) {
				break
			}
		}
		if Stop() {
			return
		}

		// join all fb
		if Global.Distr {
			mpi.AllReduceSum(d.Fb, d.Wb) // this must be done here because there might be nodes sharing boundary conditions
		}

		// point natural boundary conditions; e.g. concentrated loads
		d.PtNatBcs.AddToRhs(d.Fb, t)

		// essential boundary conditioins; e.g. constraints
		d.EssenBcs.AddToRhs(d.Fb, d.Sol)

		// debug
		if Global.Debug {
			//la.PrintVec("fb", d.Fb[:d.Ny], "%13.10f ", false)
			//panic("stop")
		}

		// find largest absolute component of fb
		largFb = la.VecLargest(d.Fb, 1)

		// save residual
		if Global.Stat {
			if Global.Summary != nil {
				Global.Summary.Resids.Append(it == 0, largFb)
			}
		}

		// check largFb value
		if it == 0 {
			// store largest absolute component of fb
			largFb0 = largFb
		} else {
			// check convergence on Lf0
			if largFb < Global.Sim.Solver.FbTol*largFb0 { // converged on fb
				break
			}
			// check convergence on fb_min
			if largFb < Global.Sim.Solver.FbMin { // converged with smallest value of fb
				break
			}
		}

		// check divergence on fb
		if it > 1 && Global.Sim.Solver.DvgCtrl {
			if largFb > prevFb {
				diverging = true
				break
			}
		}
		prevFb = largFb

		// assemble Jacobian matrix
		do_asm_fact := (it == 0 || !Global.Sim.Data.CteTg)
		if do_asm_fact {

			// assemble element matrices
			d.Kb.Start()
			for _, e := range d.Elems {
				if !e.AddToKb(d.Kb, d.Sol, it == 0) {
					break
				}
			}
			if Stop() {
				return
			}

			// debug
			if Global.DebugKb != nil {
				Global.DebugKb(d, it)
			}

			// join A and tr(A) matrices into Kb
			if Global.Root {
				d.Kb.PutMatAndMatT(&d.EssenBcs.A)
			}

			// initialise linear solver
			if d.InitLSol {
				if LogErr(d.LinSol.InitR(d.Kb, Global.Sim.LinSol.Symmetric, Global.Sim.LinSol.Verbose, Global.Sim.LinSol.Timing), "cannot initialise linear solver") {
					return
				}
				d.InitLSol = false
			}

			// perform factorisation
			LogErr(d.LinSol.Fact(), "factorisation")
			if Stop() {
				return
			}
		}

		// debug
		//KK := d.Kb.ToMatrix(nil).ToDense()
		//la.PrintMat("KK", KK, "%20.10f", false)
		//panic("stop")

		// solve for wb := δyb
		LogErr(d.LinSol.SolveR(d.Wb, d.Fb, false), "solve")
		if Stop() {
			return
		}

		// debug
		if Global.Debug {
			//la.PrintVec("wb", d.Wb[:d.Ny], "%13.10f ", false)
		}

		// update primary variables (y)
		for i := 0; i < d.Ny; i++ {
			d.Sol.Y[i] += d.Wb[i]  // y += δy
			d.Sol.ΔY[i] += d.Wb[i] // ΔY += δy
		}
		if !Global.Sim.Data.Steady {
			for _, I := range d.T1eqs {
				d.Sol.Dydt[I] = Global.DynCoefs.β1*d.Sol.Y[I] - d.Sol.Psi[I]
			}
			for _, I := range d.T2eqs {
				d.Sol.Dydt[I] = Global.DynCoefs.α4*d.Sol.Y[I] - d.Sol.Chi[I]
				d.Sol.D2ydt2[I] = Global.DynCoefs.α1*d.Sol.Y[I] - d.Sol.Zet[I]
			}
		}

		// update Lagrange multipliers (λ)
		for i := 0; i < d.Nlam; i++ {
			d.Sol.L[i] += d.Wb[d.Ny+i] // λ += δλ
		}

		// backup / restore
		if it == 0 {
			// create backup copy of all secondary variables
			for _, e := range d.ElemIntvars {
				e.BackupIvs(false)
			}
		} else {
			// recover last converged state from backup copy
			for _, e := range d.ElemIntvars {
				e.RestoreIvs(false)
			}
		}

		// update secondary variables
		for _, e := range d.Elems {
			if !e.Update(d.Sol) {
				break
			}
		}
		if Stop() {
			return
		}

		// compute RMS norm of δu and check convegence on δu
		Lδu = la.VecRmsErr(d.Wb[:d.Ny], Global.Sim.Solver.Atol, Global.Sim.Solver.Rtol, d.Sol.Y[:d.Ny])

		// message
		if Global.Sim.Data.ShowR {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}

		// stop if converged on δu
		if Lδu < Global.Sim.Solver.Itol {
			break
		}

		// check divergence on Lδu
		if it > 1 && Global.Sim.Solver.DvgCtrl {
			if Lδu > prevLδu {
				diverging = true
				break
			}
		}
		prevLδu = Lδu
	}

	// check if iterations diverged
	if it == Global.Sim.Solver.NmaxIt {
		io.PfMag("max number of iterations reached: it = %d\n", it)
		return
	}

	// success
	ok = true
	return
}

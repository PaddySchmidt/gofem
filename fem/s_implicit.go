// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/mpi"
)

// SolverImplicit solves FEM problem using an implicit procedure (with Newthon-Raphson method)
type SolverImplicit struct {
	doms []*Domain
	sum  *Summary
	dc   *DynCoefs
}

// set factory
func init() {
	solverallocators["imp"] = func(doms []*Domain, sum *Summary, dc *DynCoefs) FEsolver {
		solver := new(SolverImplicit)
		solver.doms = doms
		solver.sum = sum
		solver.dc = dc
		return solver
	}
}

func (o *SolverImplicit) Run(tf float64, dtFunc, dtoFunc fun.Func, verbose bool, dbgKb DebugKb_t) (err error) {

	// auxiliary
	md := 1.0    // time step multiplier if divergence control is on
	ndiverg := 0 // number of steps diverging

	// control
	t := o.doms[0].Sol.T
	dat := o.doms[0].Sim.Solver
	tout := t + dtoFunc.F(t, nil)

	// first output
	if o.sum != nil {
		err = o.sum.SaveDomains(t, o.doms, false)
		if err != nil {
			return chk.Err("cannot save results:\n%v", err)
		}
	}

	// time loop
	var Δt float64
	var lasttimestep bool
	for t < tf {

		// check for continued divergence
		if ndiverg >= dat.NdvgMax {
			return chk.Err("continuous divergence after %d steps reached", ndiverg)
		}

		// time increment
		Δt = dtFunc.F(t, nil) * md
		if t+Δt >= tf {
			lasttimestep = true
		}
		if Δt < dat.DtMin {
			if md < 1 {
				return chk.Err("Δt increment is too small: %g < %g", Δt, dat.DtMin)
			}
		}
		t += Δt

		// dynamic coefficients
		err = o.dc.CalcBoth(Δt)
		if err != nil {
			return chk.Err("cannot compute dynamic coefficients")
		}

		// message
		if verbose {
			if !dat.ShowR {
				io.PfWhite("%30.15f\r", t)
			}
		}

		// for all domains
		docontinue := false
		for _, d := range o.doms {

			// backup solution if divergence control is on
			if dat.DvgCtrl {
				d.backup()
			}

			// run iterations
			d.Sol.T = t
			diverging, err := run_iterations(t, Δt, d, o.dc, o.sum, dbgKb)
			if err != nil {
				return chk.Err("run_iterations failed:\n%v", err)
			}

			// restore solution and reduce time step if divergence control is on
			if dat.DvgCtrl {
				if diverging {
					if verbose {
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
			if o.sum != nil {
				err = o.sum.SaveDomains(t, o.doms, false)
				if err != nil {
					return chk.Err("cannot save results:\n%v", err)
				}
			}
			tout += dtoFunc.F(t, nil)
		}
	}
	return
}

// run_iterations solves the nonlinear problem
func run_iterations(t, Δt float64, d *Domain, dc *DynCoefs, sum *Summary, dbgKb DebugKb_t) (diverging bool, err error) {

	// zero accumulated increments
	la.VecFill(d.Sol.ΔY, 0)

	// calculate global starred vectors and interpolate starred variables from nodes to integration points
	if !d.Sim.Data.Steady {

		// compute starred vectors
		for _, I := range d.T1eqs {
			d.Sol.Psi[I] = dc.β1*d.Sol.Y[I] + dc.β2*d.Sol.Dydt[I]
		}
		for _, I := range d.T2eqs {
			d.Sol.Zet[I] = dc.α1*d.Sol.Y[I] + dc.α2*d.Sol.Dydt[I] + dc.α3*d.Sol.D2ydt2[I]
			d.Sol.Chi[I] = dc.α4*d.Sol.Y[I] + dc.α5*d.Sol.Dydt[I] + dc.α6*d.Sol.D2ydt2[I]
		}

		// set internal starred variables
		for _, e := range d.Elems {
			err = e.InterpStarVars(d.Sol)
			if err != nil {
				err = chk.Err("cannot compute starred variables:\n%v", err)
				return
			}
		}
	}

	// auxiliary variables
	var it int
	var largFb, largFb0, Lδu float64
	var prevFb, prevLδu float64
	dat := d.Sim.Solver

	// message
	if dat.ShowR {
		io.Pf("\n%13s%4s%23s%23s\n", "t", "it", "largFb", "Lδu")
		defer func() {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}()
	}

	// iterations
	for it = 0; it < dat.NmaxIt; it++ {

		// assemble right-hand side vector (fb) with negative of residuals
		la.VecFill(d.Fb, 0)
		for _, e := range d.Elems {
			err = e.AddToRhs(d.Fb, d.Sol)
			if err != nil {
				return
			}
		}

		// join all fb
		if d.Distr {
			mpi.AllReduceSum(d.Fb, d.Wb) // this must be done here because there might be nodes sharing boundary conditions
		}

		// point natural boundary conditions; e.g. concentrated loads
		d.PtNatBcs.AddToRhs(d.Fb, t)

		// essential boundary conditioins; e.g. constraints
		d.EssenBcs.AddToRhs(d.Fb, d.Sol)

		// find largest absolute component of fb
		largFb = la.VecLargest(d.Fb, 1)

		// save residual
		if d.Sim.Data.Stat {
			if sum != nil {
				sum.Resids.Append(it == 0, largFb)
			}
		}

		// check largFb value
		if it == 0 {
			// store largest absolute component of fb
			largFb0 = largFb
		} else {
			// check convergence on Lf0
			if largFb < dat.FbTol*largFb0 { // converged on fb
				break
			}
			// check convergence on fb_min
			if largFb < dat.FbMin { // converged with smallest value of fb
				break
			}
		}

		// check divergence on fb
		if it > 1 && dat.DvgCtrl {
			if largFb > prevFb {
				diverging = true
				break
			}
		}
		prevFb = largFb

		// assemble Jacobian matrix
		do_asm_fact := (it == 0 || !dat.CteTg)
		if do_asm_fact {

			// assemble element matrices
			d.Kb.Start()
			for _, e := range d.Elems {
				err = e.AddToKb(d.Kb, d.Sol, it == 0)
				if err != nil {
					return
				}
			}

			// debug
			if dbgKb != nil {
				dbgKb(d, it)
			}

			// join A and tr(A) matrices into Kb
			if d.Proc == 0 {
				d.Kb.PutMatAndMatT(&d.EssenBcs.A)
			}

			// initialise linear solver
			if d.InitLSol {
				err = d.LinSol.InitR(d.Kb, d.Sim.LinSol.Symmetric, d.Sim.LinSol.Verbose, d.Sim.LinSol.Timing)
				if err != nil {
					err = chk.Err("cannot initialise linear solver:\n%v", err)
					return
				}
				d.InitLSol = false
			}

			// perform factorisation
			err = d.LinSol.Fact()
			if err != nil {
				err = chk.Err("factorisation failed:\n%v", err)
				return
			}
		}

		// solve for wb := δyb
		err = d.LinSol.SolveR(d.Wb, d.Fb, false)
		if err != nil {
			err = chk.Err("solve failed:%v\n", err)
			return
		}

		// update primary variables (y)
		for i := 0; i < d.Ny; i++ {
			d.Sol.Y[i] += d.Wb[i]  // y += δy
			d.Sol.ΔY[i] += d.Wb[i] // ΔY += δy
		}
		if !d.Sim.Data.Steady {
			for _, I := range d.T1eqs {
				d.Sol.Dydt[I] = dc.β1*d.Sol.Y[I] - d.Sol.Psi[I]
			}
			for _, I := range d.T2eqs {
				d.Sol.Dydt[I] = dc.α4*d.Sol.Y[I] - d.Sol.Chi[I]
				d.Sol.D2ydt2[I] = dc.α1*d.Sol.Y[I] - d.Sol.Zet[I]
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
			err = e.Update(d.Sol)
			if err != nil {
				break
			}
		}

		// compute RMS norm of δu and check convegence on δu
		Lδu = la.VecRmsErr(d.Wb[:d.Ny], dat.Atol, dat.Rtol, d.Sol.Y[:d.Ny])

		// message
		if dat.ShowR {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}

		// stop if converged on δu
		if Lδu < dat.Itol {
			break
		}

		// check divergence on Lδu
		if it > 1 && dat.DvgCtrl {
			if Lδu > prevLδu {
				diverging = true
				break
			}
		}
		prevLδu = Lδu
	}

	// check if iterations diverged
	if it == dat.NmaxIt {
		io.PfMag("max number of iterations reached: it = %d\n", it)
		return
	}
	return
}

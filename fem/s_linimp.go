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

// SolverLinearImplicit solves **linear** FEM problem using an implicit procedure
// (with Newthon-Raphson method)
type SolverLinearImplicit struct {
	dom *Domain
	sum *Summary
	dc  *DynCoefs
}

// set factory of solvers
func init() {
	solverallocators["lin-imp"] = func(doms []*Domain, sum *Summary, dc *DynCoefs) FEsolver {
		if len(doms) != 1 {
			chk.Panic("SolverLinearImplicit works with one domain only")
		}
		solver := new(SolverLinearImplicit)
		solver.dom = doms[0]
		solver.sum = sum
		solver.dc = dc
		return solver
	}
}

func (o *SolverLinearImplicit) Run(tf float64, dtFunc, dtoFunc fun.Func, verbose bool, notused DebugKb_t) (err error) {

	// control
	t := o.dom.Sol.T
	tout := t + dtoFunc.F(t, nil)
	steady := o.dom.Sim.Data.Steady

	// first output
	if o.sum != nil {
		err = o.sum.SaveDomains(t, []*Domain{o.dom}, false)
		if err != nil {
			return chk.Err("cannot save results:\n%v", err)
		}
	}

	// auxiliary variables
	Y := o.dom.Sol.Y
	ψ := o.dom.Sol.Psi
	ζ := o.dom.Sol.Zet
	χ := o.dom.Sol.Chi
	dydt := o.dom.Sol.Dydt
	d2ydt2 := o.dom.Sol.D2ydt2

	// time loop
	first := true
	var Δt float64
	var lasttimestep bool
	for t < tf {

		// time increment
		Δt = dtFunc.F(t, nil)
		if t+Δt >= tf {
			lasttimestep = true
		}
		t += Δt

		// update time variable in solution array
		o.dom.Sol.T = t
		o.dom.Sol.Dt = Δt

		// dynamic coefficients
		if !steady {
			err = o.dc.CalcBoth(Δt)
			if err != nil {
				return chk.Err("cannot compute dynamic coefficients")
			}
		}

		// message
		if verbose {
			io.PfWhite("%30.15f\r", t)
		}

		// calculate global starred vectors and interpolate starred variables from nodes to integration points
		if !steady {

			// compute starred vectors
			for _, I := range o.dom.T1eqs {
				ψ[I] = o.dc.β1*Y[I] + o.dc.β2*dydt[I]
			}
			for _, I := range o.dom.T2eqs {
				ζ[I] = o.dc.α1*Y[I] + o.dc.α2*dydt[I] + o.dc.α3*d2ydt2[I]
				χ[I] = o.dc.α4*Y[I] + o.dc.α5*dydt[I] + o.dc.α6*d2ydt2[I]
			}

			// set internal starred variables
			for _, e := range o.dom.Elems {
				err = e.InterpStarVars(o.dom.Sol)
				if err != nil {
					err = chk.Err("cannot compute starred variables:\n%v", err)
					return
				}
			}
		}

		// solve linear problem
		err := solve_linear_problem(t, o.dom, o.dc, o.sum, first)
		if err != nil {
			return chk.Err("solve_linear_problem failed:\n%v", err)
		}
		first = false

		// update velocity and acceleration
		if !steady {
			for _, I := range o.dom.T1eqs {
				dydt[I] = o.dc.β1*Y[I] - ψ[I]
			}
			for _, I := range o.dom.T2eqs {
				dydt[I] = o.dc.α4*Y[I] - χ[I]
				d2ydt2[I] = o.dc.α1*Y[I] - ζ[I]
			}
		}

		// perform output
		if t >= tout || lasttimestep {
			if o.sum != nil {
				err = o.sum.SaveDomains(t, []*Domain{o.dom}, false)
				if err != nil {
					return chk.Err("cannot save results:\n%v", err)
				}
			}
			tout += dtoFunc.F(t, nil)
		}
	}
	return
}

// solve_linear_problem solves the linear problem
func solve_linear_problem(t float64, d *Domain, dc *DynCoefs, sum *Summary, first bool) (err error) {

	// assemble right-hand side vector (fb) with **negative** of residuals
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

	// assemble and factorise Jacobian matrix just once
	if first {

		// assemble element matrices
		d.Kb.Start()
		for _, e := range d.Elems {
			err = e.AddToKb(d.Kb, d.Sol, true)
			if err != nil {
				return
			}
		}

		// join A and tr(A) matrices into Kb
		if d.Proc == 0 {
			d.Kb.PutMatAndMatT(&d.EssenBcs.A)
		}

		// initialise linear solver (just once)
		if d.InitLSol {
			err = d.LinSol.InitR(d.Kb, d.Sim.LinSol.Symmetric, d.Sim.LinSol.Verbose, d.Sim.LinSol.Timing)
			if err != nil {
				err = chk.Err("cannot initialise linear solver:\n%v", err)
				return
			}
			d.InitLSol = false
		}

		// perform factorisation (always if not CteTg)
		err = d.LinSol.Fact()
		if err != nil {
			err = chk.Err("factorisation failed:\n%v", err)
			return
		}
	}

	// solve for wb
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

	// update Lagrange multipliers (λ)
	for i := 0; i < d.Nlam; i++ {
		d.Sol.L[i] += d.Wb[d.Ny+i] // λ += δλ
	}

	// update secondary variables
	for _, e := range d.Elems {
		err = e.Update(d.Sol)
		if err != nil {
			break
		}
	}
	return
}

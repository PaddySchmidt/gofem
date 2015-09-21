// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// RichardsonExtrap solves FEM problem implicitely and with Richardson's extrapolation
type RichardsonExtrap struct {

	// input
	doms []*Domain
	sum  *Summary
	dc   *DynCoefs

	// variables after big step
	Y_big []float64 // primary variables

	// auxiliary variables
	nsteps   int  // total number of steps
	naccept  int  // number of accepted steps
	nreject  int  // number of rejected steps
	ngustaf  int  // number of Gustaffson's corrections
	reject   bool // do reject this step?
	laststep bool // is last step?

	// divergence control
	ndiverg   int  // number of diverging steps (in a row)
	prevdiv   bool // previous step was diverging
	diverging bool // current step is diverging

	// time loop
	Δt    float64 // time step
	Δtcpy float64 // copy of Δt for divergence control
}

// set factory
func init() {
	solverallocators["rex"] = func(doms []*Domain, sum *Summary, dc *DynCoefs) FEsolver {
		solver := new(RichardsonExtrap)
		solver.doms = doms
		solver.sum = sum
		solver.dc = dc
		solver.Init()
		return solver
	}
}

func (o *RichardsonExtrap) Init() {

	// check
	if len(o.doms) != 1 {
		chk.Panic("RichardsonExtrap works with one domain only for now")
	}

	// auxiliary variables
	o.nsteps = 0
	o.naccept = 0
	o.nreject = 0
	o.ngustaf = 0
	o.reject = false
	o.laststep = false

	// divergence control
	o.ndiverg = 0
	o.prevdiv = false
	o.diverging = false
}

func (o *RichardsonExtrap) Run(tf float64, dtFunc, dtoFunc fun.Func, verbose bool, dbgKb DebugKb_t) (err error) {

	// constants
	dat := o.doms[0].Sim.Solver
	atol := dat.REatol
	rtol := dat.RErtol
	mmin := dat.REmmin
	mmax := dat.REmmax
	mfac := dat.REmfac

	// control
	t := o.doms[0].Sol.T
	tout := t + dtoFunc.F(t, nil)
	steady := o.doms[0].Sim.Data.Steady

	// first output
	if o.sum != nil {
		err = o.sum.SaveDomains(t, o.doms, false)
		if err != nil {
			return chk.Err("cannot save results:\n%v", err)
		}
	}

	// domain and variables
	d := o.doms[0]
	o.Y_big = make([]float64, d.Ny)

	// time loop
	o.Δt = dtFunc.F(t, nil)
	o.Δtcpy = o.Δt
	var ΔtOld, rerrOld float64
	for t < tf {

		// check for continued divergence
		if o.ndiverg >= dat.NdvgMax {
			return chk.Err("continuous divergence after %d steps reached", o.ndiverg)
		}

		// check time increment
		if o.Δt < dat.DtMin {
			return chk.Err("Δt increment is too small: %g < %g", o.Δt, dat.DtMin)
		}

		// dynamic coefficients
		if !steady {
			err = o.dc.CalcBoth(o.Δt)
			if err != nil {
				return chk.Err("cannot compute dynamic coefficients:\n%v", err)
			}
		}

		// check for maximum number of substeps
		o.nsteps += 1
		if o.nsteps >= dat.REnssmax {
			return chk.Err("RE: max number of steps reached: %d", o.nsteps)
		}

		// backup domain
		d.backup()

		// single step with Δt
		d.Sol.T = t + o.Δt
		d.Sol.Dt = o.Δt
		o.diverging, err = run_iterations(t+o.Δt, o.Δt, d, o.dc, o.sum, dbgKb)
		if err != nil {
			return chk.Err("single step with Δt: run_iterations failed:\n%v", err)
		}
		if dat.DvgCtrl {
			if o.divergence_control(d, "big step", verbose) {
				continue
			}
		}

		// save intermediate state
		for i := 0; i < d.Ny; i++ {
			o.Y_big[i] = d.Sol.Y[i]
		}

		// restore initial state
		d.restore()

		// 1st halved step
		d.Sol.T = t + o.Δt/2.0
		d.Sol.Dt = o.Δt / 2.0
		o.diverging, err = run_iterations(t+o.Δt/2.0, o.Δt/2.0, d, o.dc, o.sum, dbgKb)
		if err != nil {
			return chk.Err("1st halved step: run_iterations failed:\n%v", err)
		}
		if dat.DvgCtrl {
			if o.divergence_control(d, "1st half step", verbose) {
				continue
			}
		}

		// 2nd halved step
		d.Sol.T = t + o.Δt
		d.Sol.Dt = o.Δt
		o.diverging, err = run_iterations(t+o.Δt, o.Δt/2.0, d, o.dc, o.sum, dbgKb)
		if err != nil {
			return chk.Err("2nd halved step: run_iterations failed:\n%v", err)
		}
		if dat.DvgCtrl {
			if o.divergence_control(d, "2nd half step", verbose) {
				continue
			}
		}

		// Richardson's extrapolation error
		rerr := la.VecRmsError(d.Sol.Y, o.Y_big, atol, rtol, d.Sol.Y) / 3.0

		// step size change
		m := utl.Min(mmax, utl.Max(mmin, mfac*math.Pow(1.0/rerr, 1.0/2.0)))
		ΔtNew := m * o.Δt

		// accepted
		if rerr < 1.0 {

			// update variables
			o.naccept += 1
			t += o.Δt
			d.Sol.T = t

			// output
			if verbose {
				if !dat.ShowR {
					io.PfWhite("%30.15f\r", t)
				}
			}
			if t >= tout || o.laststep {
				if o.sum != nil {
					err = o.sum.SaveDomains(t, o.doms, false)
					if err != nil {
						return chk.Err("cannot save results:\n%v", err)
					}
				}
				tout += dtoFunc.F(t, nil)
			}

			// reached final time
			if o.laststep {
				if verbose {
					io.Pfgreen("\n\nRichardson extrapolation succeeded\n")
				}
				return
			}

			// predictive controller of Gustafsson
			if !dat.REnogus {
				if o.naccept > 1 {
					m = mfac * (o.Δt / ΔtOld) * math.Sqrt(1.0/rerr) * math.Sqrt(rerrOld/rerr)
					if m*o.Δt < ΔtNew {
						o.ngustaf += 1
					}
					ΔtNew = utl.Min(ΔtNew, m*o.Δt)
				}
				ΔtOld = o.Δt
				rerrOld = utl.Max(0.9, rerr) // 1e-2
			}

			// next step size
			if o.reject { // do not alow Δt to grow if previous was a reject
				ΔtNew = utl.Min(o.Δt, ΔtNew)
			}
			o.reject = false
			o.Δt = ΔtNew
			if t+ΔtNew-tf >= 0.0 {
				o.laststep = true
				o.Δt = tf - t
			}

			// rejected
		} else {

			// restore state
			d.restore()

			// set flags
			o.nreject += 1
			o.reject = true
			o.laststep = false

			// next step size
			o.Δt = ΔtNew
			if t+o.Δt > tf {
				o.Δt = tf - t
			}
		}
	}
	return
}

func (o *RichardsonExtrap) divergence_control(d *Domain, name string, verbose bool) (docontinue bool) {
	if o.diverging {
		if verbose {
			io.Pfred(". . . %s: iterations diverging (%2d) . . .\n", name, o.ndiverg+1)
		}
		d.restore()
		o.Δtcpy = o.Δt
		o.Δt *= 0.5
		o.ndiverg += 1
		o.nreject += 1
		o.laststep = false
		o.prevdiv = true
		return true
	}
	if o.prevdiv && false {
		o.Δt = o.Δtcpy
		o.ndiverg = 0
		o.prevdiv = false
	}
	return false
}

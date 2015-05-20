// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

type RichardsonExtrap struct {

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

func (o *RichardsonExtrap) Init(d *Domain, Dt fun.Func) {

	// backup variables
	o.Y_big = make([]float64, d.Ny)

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

	// time loop
	o.Δt = Dt.F(0, nil)
	o.Δtcpy = o.Δt
}

func (o *RichardsonExtrap) Run(d *Domain, s *Summary, DtOut fun.Func, time *float64, tf, tout float64, tidx *int) (ok bool) {

	// constants
	atol := Global.Sim.Solver.REatol
	rtol := Global.Sim.Solver.RErtol
	mmin := Global.Sim.Solver.REmmin
	mmax := Global.Sim.Solver.REmmax
	mfac := Global.Sim.Solver.REmfac

	// time loop
	t := *time
	defer func() { *time = t }()
	var ΔtOld, rerrOld float64
	for t < tf {

		// check for continued divergence
		if LogErrCond(o.ndiverg >= Global.Sim.Solver.NdvgMax, "continuous divergence after %d steps reached", o.ndiverg) {
			return
		}

		// check time increment
		if LogErrCond(o.Δt < Global.Sim.Solver.DtMin, "Δt increment is too small: %g < %g", o.Δt, Global.Sim.Solver.DtMin) {
			return
		}

		// compute dynamic coefficients
		if LogErr(Global.DynCoefs.CalcBoth(o.Δt), "cannot compute dynamic coefficients") {
			return
		}

		// check for maximum number of substeps
		o.nsteps += 1
		if LogErrCond(o.nsteps >= Global.Sim.Solver.REnssmax, "RE: max number of steps reached: %d", o.nsteps) {
			return
		}

		// backup domain
		d.backup()

		// single step with Δt
		d.Sol.T = t + o.Δt
		o.diverging, ok = run_iterations(t+o.Δt, o.Δt, d, s)
		if !ok {
			return
		}
		if Global.Sim.Solver.DvgCtrl {
			if o.divergence_control(d, "big step") {
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
		o.diverging, ok = run_iterations(t+o.Δt/2.0, o.Δt/2.0, d, s)
		if !ok {
			break
		}
		if Global.Sim.Solver.DvgCtrl {
			if o.divergence_control(d, "1st half step") {
				continue
			}
		}

		// 2nd halved step
		d.Sol.T = t + o.Δt
		o.diverging, ok = run_iterations(t+o.Δt, o.Δt/2.0, d, s)
		if !ok {
			break
		}
		if Global.Sim.Solver.DvgCtrl {
			if o.divergence_control(d, "2nd half step") {
				continue
			}
		}

		// Richardson's extrapolation error
		rerr := la.VecRmsError(d.Sol.Y, o.Y_big, atol, rtol, d.Sol.Y) / 3.0

		// step size change
		m := min(mmax, max(mmin, mfac*math.Pow(1.0/rerr, 1.0/2.0)))
		ΔtNew := m * o.Δt

		// accepted
		if rerr < 1.0 {

			// update variables
			o.naccept += 1
			t += o.Δt
			d.Sol.T = t

			io.Pfgreen("accepted t = %v\n", t)

			// output
			if Global.Verbose {
				if !Global.Sim.Data.ShowR && !Global.Debug {
					io.PfWhite("time     = %g\r", t)
				}
			}
			if true {
				//if t >= tout || o.laststep {
				s.OutTimes = append(s.OutTimes, t)
				if !d.Out(*tidx) {
					return
				}
				io.Pforan("doing output with t = %v\n", t)
				tout += DtOut.F(t, nil)
				*tidx += 1
			}

			// reached final time
			if o.laststep {
				return true
			}

			// predictive controller of Gustafsson
			if !Global.Sim.Solver.REnogus {
				if o.naccept > 1 {
					m = mfac * (o.Δt / ΔtOld) * math.Sqrt(1.0/rerr) * math.Sqrt(rerrOld/rerr)
					if m*o.Δt < ΔtNew {
						o.ngustaf += 1
					}
					ΔtNew = min(ΔtNew, m*o.Δt)
				}
				ΔtOld = o.Δt
				rerrOld = max(0.9, rerr) // 1e-2
			}

			// next step size
			if o.reject { // do not alow Δt to grow if previous was a reject
				ΔtNew = min(o.Δt, ΔtNew)
			}
			o.reject = false
			o.Δt = ΔtNew
			if t+ΔtNew-tf >= 0.0 {
				o.laststep = true
				o.Δt = tf - t
			}

			// rejected
		} else {

			io.Pfyel("rejected t = %v\n", t)

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
	return true
}

func (o *RichardsonExtrap) divergence_control(d *Domain, name string) (docontinue bool) {
	if o.diverging {
		if Global.Verbose {
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

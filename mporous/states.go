// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

// State holds state variables for porous media with liquid and gas
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type State struct {
	A_ns0 float64 // 1 initial partial fraction of solids
	A_sl  float64 // 2 liquid saturation
	A_ρL  float64 // 3 real (intrinsic) density of liquid
	A_ρG  float64 // 4 real (intrinsic) density of gas
	A_Δpc float64 // 5 step increment of capillary pressure
	A_wet bool    // 6 wetting flag
}

// GetCopy returns a copy of State
func (o State) GetCopy() *State {
	return &State{
		o.A_ns0, // 1
		o.A_sl,  // 2
		o.A_ρL,  // 3
		o.A_ρG,  // 4
		o.A_Δpc, // 5
		o.A_wet, // 6
	}
}

// Set sets this State with another State
func (o *State) Set(s *State) {
	o.A_ns0 = s.A_ns0 // 1
	o.A_sl = s.A_sl   // 2
	o.A_ρL = s.A_ρL   // 3
	o.A_ρG = s.A_ρG   // 4
	o.A_Δpc = s.A_Δpc // 5
	o.A_wet = s.A_wet // 6
}

// LsVars hold data for liquid-solid computations
type LsVars struct {
	ρl, ρ, p, Cpl, Cvs                      float64
	dρdpl, dpdpl, dCpldpl, dCvsdpl, dklrdpl float64
	dρldusM, dρdusM, dCpldusM               float64
}

// CalcLs calculates variables for liquid-solid simulations
func (o Model) CalcLs(res *LsVars, sta *State, pl, divus float64, derivs bool) (err error) {

	// auxiliary
	ns0 := sta.A_ns0
	sl := sta.A_sl
	ρL := sta.A_ρL
	Cl := o.Cl
	ρS := o.RhoS0

	// n variables; Eqs (13) and (28) of [1]
	ns := (1.0 - divus) * ns0
	nf := 1.0 - ns
	nl := nf * sl

	// ρ variables; Eq (13) of [1]
	ρs := ns * ρS
	res.ρl = nl * ρL
	res.ρ = res.ρl + ρs

	// capillary pressure and pore-fluid pressure
	pc := -pl
	res.p = pl * sl // Eq. (16) of [1]

	// moduli
	Ccb, e := o.Ccb(sta, pc)
	if e != nil {
		return e
	}
	res.Cpl = nf * (sl*Cl - ρL*Ccb) // Eq (32a) of [1]
	res.Cvs = sl * ρL               // Eq (32b) of [1]

	// derivatives
	if derivs {

		// Ccd
		Ccd, e := o.Ccd(sta, pc)
		if e != nil {
			return e
		}

		// derivatives w.r.t pl
		res.dρdpl = nf * (sl*Cl - ρL*Ccb)        // Eq (A.9) of [1]
		res.dpdpl = sl + pc*Ccb                  // Eq (A.11) of [1]
		res.dCpldpl = nf * (ρL*Ccd - 2.0*Ccb*Cl) // Eq (A.2) of[1]
		res.dCvsdpl = sl*Cl - Ccb*ρL             // Eq (A.4) of [1]
		res.dklrdpl = -o.Cnd.DklrDsl(sl) * Ccb   // Eq (A.7) of [1]

		// derivatives w.r.t us (multipliers only)
		res.dρldusM = sl * ρL * ns0
		res.dρdusM = (sl*ρL - ρS) * ns0       // Eq (A.10) of [1]
		res.dCpldusM = (sl*Cl - ρL*Ccb) * ns0 // Eq (A.3) of [1]
	}
	return
}

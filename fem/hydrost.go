// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
)

// HydroStatic computes water pressure (pl) and intrinsic liquid density (ρL)
// based on the following model
//
//    ρL = ρL0 + Cl・pl   thus   dρL/dpl = Cl
//
//    Z(z) = zmax + T・(z - zmax)   with 0 ≤ T ≤ 1
//    dZ   = (z - zmax)・dT
//    dpl  = ρL(pl)・g・(-dZ)
//    dpl  = ρL(pl)・g・(zmax - z)・dT
//    Δz   = zmax - z
//
//            / dpl/dT \   / ρL(pl)・g・Δz \
//    dY/dT = |        | = |               |
//            \ dρL/dT /   \ Cl・dpl/dT    /
//
type HydroStatic struct {
	zwater float64
	ρL0    float64
	Cl     float64
	g      float64
	fcn    ode.Cb_fcn
	Jac    ode.Cb_jac
	sol    ode.ODE
}

// Init initialises this structure
func (o *HydroStatic) Init(waterLevel, waterRho0, waterBulk, g float64) {

	// basic data
	o.zwater = waterLevel
	o.ρL0 = waterRho0
	o.Cl = o.ρL0 / waterBulk
	o.g = g

	// x := {pl, ρL}
	o.fcn = func(f []float64, x float64, y []float64, args ...interface{}) error {
		Δz := args[0].(float64)
		//ρL := o.ρL0
		ρL := y[1]
		f[0] = ρL * o.g * Δz // dpl/dT
		f[1] = o.Cl * f[0]   // dρL/dT
		return nil
	}

	o.Jac = func(dfdy *la.Triplet, x float64, y []float64, args ...interface{}) error {
		if dfdy.Max() == 0 {
			dfdy.Init(2, 2, 4)
		}
		Δz := args[0].(float64)
		dfdy.Start()
		dfdy.Put(0, 0, 0)
		dfdy.Put(0, 1, o.g*Δz)
		dfdy.Put(1, 0, 0)
		dfdy.Put(1, 1, o.Cl*o.g*Δz)
		return nil
	}

	silent := true
	o.sol.Init("Radau5", 2, o.fcn, o.Jac, nil, nil, silent)
	o.sol.Distr = false // must be sure to disable this; otherwise it causes problems in parallel runs
}

// Calc computes pressure and density
func (o HydroStatic) Calc(z float64) (pl, ρL float64, err error) {
	Δz := o.zwater - z
	y := []float64{0, o.ρL0} // pl0, ρL0
	err = o.sol.Solve(y, 0, 1, 1, false, Δz)
	if err != nil {
		err = chk.Err("HydroStatic failed when calculating pressure using ODE solver: %v", err)
		return
	}
	return y[0], y[1], nil
}

// SetHydroSt sets the initial state to a hydrostatic condition
func (o *Domain) SetHydroSt(stg *inp.Stage) (err error) {

	// set Sol
	ndim := o.Sim.Ndim
	for _, n := range o.Nodes {
		z := n.Vert.C[ndim-1]
		dof := n.GetDof("pl")
		if dof != nil {
			pl, _, err := o.HydSta.Calc(z)
			if err != nil {
				return chk.Err("hydrost: cannot compute pl:\n%v", err)
			}
			o.Sol.Y[dof.Eq] = pl
		}
	}

	// set elements
	for _, e := range o.ElemIntvars {

		// build map with pressures @ ips
		coords := e.Ipoints()
		nip := len(coords)
		pl := make([]float64, nip)
		ρL := make([]float64, nip)
		for i := 0; i < nip; i++ {
			z := coords[i][ndim-1]
			pl[i], ρL[i], err = o.HydSta.Calc(z)
			if err != nil {
				return chk.Err("hydrost: cannot compute pl and ρL:\n%v", err)
			}
		}
		ivs := map[string][]float64{"pl": pl, "ρL": ρL}

		// set element's states
		err = e.SetIniIvs(o.Sol, ivs)
		if err != nil {
			return chk.Err("hydrost: element's internal values setting failed:\n%v", err)
		}
	}
	return
}

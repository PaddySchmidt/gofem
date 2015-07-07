// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_rjoint01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("rjoint01. curved line in 3D")

	// initialisation
	defer End()
	if !Start("data/rjoint01.sim", true, chk.Verbose) {
		tst.Errorf("Start failed\n")
		return
	}

	// callback to check consistent tangent operators
	eid := 2 // rjoint element
	if true {
		defer rjoint_DebugKb(&testKb{
			tst: tst, eid: eid, tol: 1e-8, verb: chk.Verbose,
			ni: -1, nj: -1, itmin: 1, itmax: -1, tmin: -1, tmax: -1,
		})()
	}

	// run simulation
	if !Run() {
		tst.Errorf("Run failed\n")
		return
	}

	// plot
	//if true {
	if false {

		// allocate domain
		sum := ReadSum(Global.Dirout, Global.Fnkey)
		dom := NewDomain(Global.Sim.Regions[0], false)
		if !dom.SetStage(0, Global.Sim.Stages[0], false) {
			tst.Errorf("SetStage failed\n")
			return
		}
		ele := dom.Elems[eid].(*Rjoint)
		ipd := ele.OutIpsData()

		// load results from file
		n := len(sum.OutTimes)
		mτ := make([]float64, n)
		ωpb := make([]float64, n)
		for i, _ := range sum.OutTimes {
			if !dom.In(sum, i, true) {
				tst.Errorf("cannot read solution\n")
				return
			}
			for _, dat := range ipd {
				res := dat.Calc(dom.Sol)
				mτ[i] = -res["tau"]
				ωpb[i] = res["ompb"]
			}
		}

		// plot
		plt.Plot(ωpb, mτ, "'b-', marker='o', clip_on=0")
		plt.Gll("$\\bar{\\omega}_p$", "$-\\tau$", "")
		plt.Show()
	}
}

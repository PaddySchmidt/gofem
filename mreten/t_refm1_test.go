// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_refm1a(tst *testing.T) {

	doplot := false
	//doplot := true
	chk.PrintTitle("refm1a")

	mdl := GetModel("testsim", "mat1", "ref-m1", false)
	mdl.Init(mdl.GetPrms(true))

	pc0 := -5.0
	sl0 := 1.0
	pcf := 20.0
	nptsA := 41
	nptsB := 11

	if doplot {
		plt.Reset()
		Plot(mdl, pc0, sl0, pcf, nptsA, "'b.-'", "'r+-'", "ref-m1_drying")
	}

	tolCc := 1e-17
	tolD1a, tolD1b := 1e-11, 1e-11
	tolD2a, tolD2b := 1e-12, 1e-10
	Check(tst, mdl, pc0, sl0, pcf, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0}, 1e-7, doplot)

	slf, err := Update(mdl, pc0, sl0, pcf-pc0)
	if err != nil {
		tst.Errorf("update failed: %v\n", err)
		return
	}

	if doplot {
		Plot(mdl, pcf, slf, pc0, nptsA, "'b*-'", "'r+:'", "ref-m1_wetting")
	}

	tolD1b = 1e-4
	tolD2a, tolD2b = 1e-11, 1e-10
	Check(tst, mdl, pcf, slf, pc0, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0}, 1e-7, doplot)

	if doplot {
		PlotEnd(true)
	}
}

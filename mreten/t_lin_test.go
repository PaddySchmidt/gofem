// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_lin01(tst *testing.T) {

	doplot := false
	//doplot := true
	chk.PrintTitle("lin01")

	mdl := GetModel("testsim", "mat1", "lin", false)
	mdl.Init(mdl.GetPrms(true))

	pc0 := -1.0
	sl0 := 1.0
	pcf := 3.0
	nptsA := 11
	nptsB := 11

	if doplot {
		plt.Reset()
		Plot(mdl, pc0, sl0, pcf, nptsA, "'b.-'", "'r+-'", "lin")
	}

	tolCc := 1e-13
	tolD1a, tolD1b := 1e-13, 1e-17
	tolD2a, tolD2b := 1e-13, 1e-17
	Check(tst, mdl, pc0, sl0, pcf, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0.2}, 1e-7, doplot)

	if doplot {
		PlotEnd(true)
	}
}

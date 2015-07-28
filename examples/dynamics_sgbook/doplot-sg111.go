// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"math"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func solution_uy(t, ta float64) float64 {
	if t < ta {
		return 0.441*math.Sin(math.Pi*t/ta) - 0.216*math.Sin(2.0*math.Pi*t/ta)
	}
	return -0.432 * math.Sin(2*math.Pi*(t-ta))
}

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "sg111", ".sim", true)

	// results
	out.Start(filename, 0, 0)
	out.Define("tip", out.N{1})
	out.LoadResults(nil)

	// plot FEM results
	out.Plot("t", "uy", "tip", plt.Fmt{C: "r", Ls: "None", M: ".", L: "gofem"}, -1)

	// analytical solution
	tAna := utl.LinSpace(0, 5, 101)
	uyAna := make([]float64, len(tAna))
	for i, t := range tAna {
		uyAna[i] = solution_uy(t, 1.0)
	}

	// save
	plt.SetForPng(0.8, 400, 200)
	out.Draw("/tmp", fnkey+".png", false, func(i, j, n int) {
		plt.Plot(tAna, uyAna, "'g-', clip_on=0, label='analytical'")
	})
}

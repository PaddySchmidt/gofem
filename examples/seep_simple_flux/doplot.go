// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.Args0toFilename("d2-simple-flux", ".sim", true)

	// start analysis process
	out.Extrap = []string{"nwlx", "nwly"}
	out.Start(filename, 0, 0)

	// define entities
	out.Define("top-middle", out.At{5, 3})
	out.Define("section-A", out.N{-1})
	out.Define("section-B", out.Along{{0, 0}, {10, 0}})

	// load results
	out.LoadResults(nil)

	// compute water discharge along section-A
	nwlx_TM := out.GetRes("ex_nwlx", "top-middle", -1)
	Q := out.Integrate("ex_nwlx", "section-A", "y", -1)
	io.PfYel("Q = %g m³/s [answer: 0.0003]\n", Q)

	// plot
	kt := len(out.Times) - 1
	out.Splot("")
	out.Plot("pl", "y", "section-A", plt.Fmt{L: "t=0"}, 0)
	out.Plot("pl", "y", "section-A", plt.Fmt{L: io.Sf("t=%g", out.Times[kt])}, kt)
	out.Splot("")
	out.Plot("x", "pl", "section-B", plt.Fmt{L: "t=0"}, 0)
	out.Plot("x", "pl", "section-B", plt.Fmt{L: io.Sf("t=%g", out.Times[kt])}, kt)
	out.Splot("")
	out.Plot("t", nwlx_TM, "top-middle", plt.Fmt{}, -1)
	out.Csplot.Ylbl = "$n_{\\ell}\\cdot w_{\\ell x}$"

	// save
	plt.SetForPng(1.5, 400, 200)
	out.Draw("/tmp", "seep_simple_flux_"+fnkey+".png", false, func(i, j, nplots int) {
		if i == 2 && j == 1 {
			plt.Plot([]float64{0, 10}, []float64{10, 9}, "'k--'")
		}
	})
}

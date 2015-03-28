// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfn := "d2-simple-flux"
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// start analysis process
	out.Extrap = []string{"nwlx", "nwly"}
	out.Start(simfn, 0, 0)

	// define entities
	out.Define("top-middle", out.At{5, 3})
	out.Define("section-A", out.N{-1})
	out.Define("section-B", out.Along{{0, 0}, {10, 0}})

	// load results
	out.LoadResults(nil)

	// compute water discharge along section-A
	nwlx_TM := out.GetRes("ex_nwlx", "top-middle", -1)
	nwlx := out.GetRes("ex_nwlx", "section-A", -1)
	_, y, _ := out.GetXYZ("ex_nwlx", "section-A")
	io.Pforan("y = %v\n", y)
	io.Pforan("nwlx = %v\n", nwlx)
	Q := num.Trapz(y, nwlx)
	io.PfYel("Q = %g m³/s [answer: 0.0003]\n", Q)
	/*
		ids, _ := out.GetIds("section-A")
		nv := len(ids)
		s, f := make([]float64, nv), make([]float64, nv)
		for i, id := range ids {
			s[i] = out.Dom.Msh.Verts[id].C[1]
			f[i] = out.ExVals[id]["nwlx"]
			io.Pforan("%2d : y=%.2f nwlx=%v\n", id, s[i], f[i])
		}
		Q := num.Trapz(s, f)
		io.PfYel("Q = %g m³/s [answer: 0.0003]\n", Q)
	*/

	// plot
	kt := len(out.T) - 1
	out.Splot("")
	out.Plot("pl", "y", "section-A", plt.Fmt{L: "t=0"}, 0)
	out.Plot("pl", "y", "section-A", plt.Fmt{L: io.Sf("t=%g", out.T[kt])}, kt)
	out.Splot("")
	out.Plot("x", "pl", "section-B", plt.Fmt{L: "t=0"}, 0)
	out.Plot("x", "pl", "section-B", plt.Fmt{L: io.Sf("t=%g", out.T[kt])}, kt)
	out.Splot("")
	out.Plot("t", nwlx_TM, "top-middle", plt.Fmt{}, -1)
	out.Csplot.Ylbl = "$n_{\\ell}\\cdot w_{\\ell x}$"

	// show
	if true {
		out.Draw("", "", true, func(i, j, nplots int) {
			if i == 2 && j == 1 {
				plt.Plot([]float64{0, 10}, []float64{10, 9}, "'k--'")
			}
		})
	}
}

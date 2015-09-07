// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "onepulse-qua9co.sim", ".sim", true)

	// start analysis process
	out.Start(filename, 0, 0)

	// define entities
	out.Define("A B C D E", out.N{-5, -4, -3, -2, -1})
	out.Define("a b c d e", out.P{{15, 8}, {13, 8}, {8, 8}, {4, 8}, {0, 0}})
	out.Define("left-side", out.Along{{0, 0}, {0, 3}})

	// load results
	out.LoadResults(nil)

	// styles
	me := 10
	S := []plt.Fmt{
		plt.Fmt{C: "b", M: "*", Me: me},
		plt.Fmt{C: "g", M: "o", Me: me},
		plt.Fmt{C: "m", M: "x", Me: me},
		plt.Fmt{C: "orange", M: "+", Me: me},
		plt.Fmt{C: "r", M: "^", Me: me},
	}

	// pl
	out.Splot("liquid pressure")
	for i, l := range []string{"A", "B", "C", "D", "E"} {
		out.Plot("t", "pl", l, S[i], -1)
	}

	// uy
	out.Splot("displacements")
	for i, l := range []string{"A", "B", "C", "D", "E"} {
		out.Plot("t", "uy", l, S[i], -1)
	}

	// sl
	out.Splot("liquid saturation")
	for i, l := range []string{"a", "b", "c", "d", "e"} {
		out.Plot("t", "sl", l, S[i], -1)
	}

	// sy
	out.Splot("stresses")
	for i, l := range []string{"a", "b", "c", "d", "e"} {
		out.Plot("t", "sy", l, S[i], -1)
	}

	// pl
	out.Splot("pressure along column")
	I, _ := utl.GetITout(out.Times, []float64{0, 1000, 2000, 3000, 4000}, 0.1)
	for _, i := range I {
		t := out.Times[i]
		out.Plot("pl", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
	}

	// empty
	out.Splot("")

	// save
	plt.SetForPng(1.2, 800, 200)
	out.Draw("/tmp", "up_3mcolumn_dessication_"+fnkey+".png", false, nil)
}

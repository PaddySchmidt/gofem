// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.Args0toFilename("rjoint01", ".sim", true)

	// results
	out.Start(filename, 0, 0)
	out.Define("p0", out.P{{-3, 0}}) // -3=tag, first ip
	out.Define("p1", out.P{{-3, 1}}) // -3=tag, second ip
	out.Define("p2", out.P{{-3, 2}}) // -3=tag, third ip
	out.LoadResults(nil)
	mtau0 := out.GetRes("tau", "p0", -1)
	mtau1 := out.GetRes("tau", "p1", -1)
	mtau2 := out.GetRes("tau", "p2", -1)
	for i := 0; i < len(out.Times); i++ {
		mtau0[i] *= -1
		mtau1[i] *= -1
		mtau2[i] *= -1
	}

	// plot FEM results
	out.Plot("ompb", mtau0, "p0", plt.Fmt{C: "r", Ls: "-", M: "."}, -1)
	out.Plot("ompb", mtau1, "p1", plt.Fmt{C: "g", Ls: "-", M: "."}, -1)
	out.Plot("ompb", mtau2, "p2", plt.Fmt{C: "b", Ls: "-", M: "."}, -1)
	out.Csplot.Ylbl = "$-\\tau$"

	// save
	plt.SetForPng(0.8, 400, 200)
	out.Draw("/tmp", fnkey+".png", false, nil)
}

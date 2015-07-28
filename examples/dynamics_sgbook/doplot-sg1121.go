// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"encoding/json"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "sg1121", ".sim", true)

	// results
	out.Start(filename, 0, 0)
	out.Define("A", out.N{30})
	out.LoadResults(nil)

	// plot FEM results
	out.Plot("t", "uy", "A", plt.Fmt{C: "k", Ls: "-", L: "gofem"}, -1)

	// old results
	b, err := io.ReadFile("cmp/sg1121gofemold.json")
	if err != nil {
		io.PfRed("cannot read comparison file\n")
		return
	}
	var gofemold struct {
		Time, Uy30 []float64
	}
	err = json.Unmarshal(b, &gofemold)
	if err != nil {
		io.PfRed("cannot unmarshal comparison file\n")
		return
	}

	// mechsys results
	_, res, err := io.ReadTable("cmp/sg1121mechsysN30.cmp")
	if err != nil {
		io.PfRed("cannot read mechsys comparison file\n")
		return
	}

	// save
	plt.SetForPng(0.8, 400, 200)
	out.Draw("/tmp", fnkey+".png", false, func(i, j, n int) {
		plt.Plot(gofemold.Time, gofemold.Uy30, "'r-', lw=2, label='gofemOld'")
		plt.Plot(res["Time"], res["uy"], "'b-', label='mechsys'")
	})
}

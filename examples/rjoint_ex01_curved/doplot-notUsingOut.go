// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "rjoint01", ".sim", true)

	// fem
	if !fem.Start(filename, false, false, false) {
		io.PfRed("Start failed\n")
		return
	}
	dom, sum, ok := fem.AllocSetAndInit(0, true, true)
	if !ok {
		io.PfRed("AllocSetAndInit failed\n")
		return
	}

	// rjoint element
	eid := 2
	ele := dom.Elems[eid].(*fem.Rjoint)
	ipd := ele.OutIpsData()

	// load results from file
	n := len(sum.OutTimes)
	mtau0 := make([]float64, n)
	mtau1 := make([]float64, n)
	mtau2 := make([]float64, n)
	ompb0 := make([]float64, n)
	ompb1 := make([]float64, n)
	ompb2 := make([]float64, n)
	for i, _ := range sum.OutTimes {
		if !dom.In(sum, i, true) {
			io.PfRed("cannot read solution\n")
			return
		}
		res0 := ipd[0].Calc(dom.Sol)
		res1 := ipd[1].Calc(dom.Sol)
		res2 := ipd[2].Calc(dom.Sol)
		mtau0[i] = -res0["tau"]
		mtau1[i] = -res1["tau"]
		mtau2[i] = -res2["tau"]
		ompb0[i] = res0["ompb"]
		ompb1[i] = res1["ompb"]
		ompb2[i] = res2["ompb"]
	}

	// plot
	plt.SetForPng(0.8, 400, 200)
	plt.Plot(ompb0, mtau0, "'r-', marker='.', label='p0', clip_on=0")
	plt.Plot(ompb1, mtau1, "'g-', marker='.', label='p1', clip_on=0")
	plt.Plot(ompb2, mtau2, "'b-', marker='.', label='p2', clip_on=0")
	plt.Gll("$\\bar{\\omega}_p$", "$-\\tau$", "")
	plt.SaveD("/tmp", fnkey+".png")
}

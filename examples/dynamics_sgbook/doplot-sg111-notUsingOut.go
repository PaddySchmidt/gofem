// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"math"

	"github.com/cpmech/gofem/fem"
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
	filename, fnkey := io.Args0toFilename("sg111", ".sim", true)

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

	// selected node and dof index
	nidx := 1
	didx := 1

	// gofem
	ntout := len(sum.OutTimes)
	uy := make([]float64, ntout)
	for tidx, _ := range sum.OutTimes {

		// read results from file
		if !dom.In(sum, tidx, true) {
			io.PfRed("plot_spo751: cannot read solution\n")
			return
		}

		// collect results for load versus time plot
		nod := dom.Nodes[nidx]
		eq := nod.Dofs[didx].Eq
		uy[tidx] = dom.Sol.Y[eq]

		// check
		if math.Abs(dom.Sol.T-sum.OutTimes[tidx]) > 1e-14 {
			io.PfRed("output times do not match time in solution array\n")
			return
		}
	}

	// plot fem results
	plt.SetForPng(0.8, 400, 200)
	plt.Plot(sum.OutTimes, uy, "'ro-', clip_on=0, label='gofem'")

	// analytical solution
	tAna := utl.LinSpace(0, 5, 101)
	uyAna := make([]float64, len(tAna))
	for i, t := range tAna {
		uyAna[i] = solution_uy(t, 1.0)
	}
	plt.Plot(tAna, uyAna, "'g-', clip_on=0, label='analytical'")

	// save
	plt.Gll("$t$", "$u_y$", "")
	plt.SaveD("/tmp", fnkey+".png")
}

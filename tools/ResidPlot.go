// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
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

func read_summary(simfn string) (*utl.DblSlist, string) {
	if simfn == "" {
		return nil, ""
	}
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}
	out.Start(simfn, 0, 0)
	return &out.Sum.Resids, io.FnKey(simfn)
}

func count_iters(resid *utl.DblSlist) (N []float64) {
	P := resid.Ptrs
	for i := 0; i < len(P)-1; i++ {
		n := P[i+1] - P[i]
		N = append(N, float64(n-1))
	}
	return
}

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data
	simfnA, fnkA := io.ArgToFilename(0, "o2elastCO", ".sim", true)
	skip := io.ArgToInt(1, 0)
	simfnB, fnkB := io.ArgToFilename(2, "", ".sim", false)
	labelA := io.ArgToString(3, "")
	labelB := io.ArgToString(4, "")

	// print input data
	io.Pf("\n%s\n", io.ArgsTable(
		"simulation filename", "simfnA", simfnA,
		"number of initial increments to skip", "skip", skip,
		"simulation filename for comparison", "simfnB", simfnB,
		"label for histogram", "labelA", labelA,
		"label for histogram", "labelB", labelB,
	))

	// read residuals
	residA, fnkA := read_summary(simfnA)
	residB, fnkB := read_summary(simfnB)

	// residuals: it => residuals
	io.Pf("\nResiduals A\n")
	io.Pf("============\n")
	residA.Print("%10.2e")
	if simfnB != "" {
		io.Pf("\nResiduals B\n")
		io.Pf("============\n")
		residB.Print("%10.2e")
	}
	io.Pf("\n")

	// plot convergence curves
	plot_conv_curve(fnkA, skip, residA)
	if simfnB != "" {
		plot_conv_curve(fnkB, skip, residB)
	}

	// plot histogram
	io.Pf("\n")
	X := [][]float64{count_iters(residA)}
	labels := []string{fnkA}
	if labelA != "" {
		labels[0] = labelA
	}
	if simfnB != "" {
		X = append(X, count_iters(residB))
		labels = append(labels, fnkB)
		if labelB != "" {
			labels[1] = labelB
		}
	}
	plt.Reset()
	plt.SetForEps(0.75, 300)
	plt.Hist(X, labels, "")
	plt.Gll("number of iterations", "counts", "")
	plt.SaveD("/tmp", "gofem_residplot_"+fnkA+"_"+fnkB+"_hist.eps")
}

func plot_conv_curve(fnk string, skip int, resid *utl.DblSlist) {
	R := resid.Vals
	P := resid.Ptrs
	plt.Reset()
	plt.SetForEps(0.75, 250)
	for i := 0; i < len(P)-1; i++ {
		if i >= skip {
			n := P[i+1] - P[i]
			x := make([]float64, n)
			y := make([]float64, n)
			k := 0
			for j := P[i]; j < P[i+1]; j++ {
				x[k] = float64(k)
				y[k] = math.Log10(R[j])
				k += 1
			}
			plt.Plot(x, y, "")
		}
	}
	plt.Gll("iteration index", "$\\mathrm{log_{10}}(R)$", "")
	plt.SaveD("/tmp", "gofem_residplot_"+fnk+"_curves.eps")
}

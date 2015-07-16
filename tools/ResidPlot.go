// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"flag"
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

	// input data
	simfnA := "o2elastCO"
	skip := 0
	simfnB := ""
	labelA := ""
	labelB := ""
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfnA = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		skip = io.Atoi(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		simfnB = flag.Arg(2)
	}
	if len(flag.Args()) > 3 {
		labelA = flag.Arg(3)
	}
	if len(flag.Args()) > 4 {
		labelB = flag.Arg(4)
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfnA = %30s // simulation filename\n", simfnA)
	io.Pf("  skip   = %30d // number of initial increments to skip\n", skip)
	io.Pf("  simfnB = %30s // simulation filename for comparison\n", simfnB)
	io.Pf("  labelA = %30s // label for histogram\n", labelA)
	io.Pf("  labelB = %30s // label for histogram\n", labelB)
	io.Pf("\n")

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

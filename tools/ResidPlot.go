// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"math"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfnA := "o2elastCO"
	skip := 0
	simfnB := ""
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

	// check extension
	if io.FnExt(simfnA) == "" {
		simfnA += ".sim"
	}
	fnkA := io.FnKey(simfnA)

	// B file
	var fnkB string
	if simfnB != "" {
		if io.FnExt(simfnB) == "" {
			simfnB += ".sim"
		}
		fnkB = io.FnKey(simfnB)
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfnA = %30s // simulation filename\n", simfnA)
	io.Pf("  skip   = %30d // number of initial increments to skip\n", skip)
	io.Pf("  simfnB = %30s // simulation filename for comparison\n", simfnB)
	io.Pf("\n")

	// start analysis process
	out.Start(simfnA, 0, 0)

	// residuals: it => residuals
	io.Pf("\nResiduals\n")
	io.Pf("==========\n")
	R := out.Sum.Resids.Vals
	P := out.Sum.Resids.Ptrs
	out.Sum.Resids.Print("%10.2e")

	// plot convergence curves
	plot_conv_curve(fnkA, skip, R, P)
	if simfnB != "" {
		plot_conv_curve(fnkB, skip, R, P)
	}

	// plot histrogram
	plt.Reset()

}

func plot_conv_curve(fnk string, skip int, R []float64, P []int) {
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

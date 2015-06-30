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
	simfn := "o2Elast"
	skip := 0
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		skip = io.Atoi(flag.Arg(1))
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}
	fnk := io.FnKey(simfn)

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfn = %30s // simulation filename\n", simfn)
	io.Pf("  skip  = %30s // number of initial increments to skip\n", skip)
	io.Pf("\n")

	// start analysis process
	out.Start(simfn, 0, 0)

	// residuals: it => residuals
	io.Pf("\nResiduals\n")
	io.Pf("==========\n")
	R := out.Sum.Resids.Vals
	P := out.Sum.Resids.Ptrs
	out.Sum.Resids.Print("%10.2e")

	// plot convergence curves
	plt.SetForEps(0.75, 300)
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

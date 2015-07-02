// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	// test title
	//verbose()
	chk.PrintTitle("plot01")

	// constants
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	simfn := "p02.sim"

	// run FE simulation
	defer fem.End()
	if !fem.Start(datadir+simfn, true, chk.Verbose) {
		chk.Panic("cannot start FE simulation")
	}
	if !fem.Run() {
		chk.Panic("cannot run FE simulation")
	}

	// start post-processing
	Start(datadir+simfn, 0, 0)

	// define entities
	Define("A", N{1})
	Define("B", At{0, 1})
	//Define("left", Along{{0, 0}, {0, 10}})
	Define("left", AlongY{0}) // 0 => x_cte

	// load results
	LoadResults(nil)

	plA := GetRes("pl", "A", 0)

	Splot("liquid pressure")
	Plot("t", "pl", "B", plt.Fmt{C: "b", M: "."}, -1)
	Plot("t", plA, "A", plt.Fmt{C: "r", M: "."}, -1)

	Splot("")
	Plot("pl", "pl", "A", plt.Fmt{C: "k", M: "o"}, -1)

	Splot("")
	Plot("pl", "y", "left", plt.Fmt{C: "b", M: "o"}, 0)
	Plot("pl", "y", "left", plt.Fmt{C: "g", M: "o"}, -1)

	Splot("")
	io.Pforan("T = %v\n", Times)
	last := len(Times) - 1
	Plot("y", "pl", "left", plt.Fmt{C: "b", M: "o", L: io.Sf("t=%g", Times[0])}, 0)
	Plot("y", "pl", "left", plt.Fmt{C: "m", M: "*", Lw: 2, L: io.Sf("t=%g", Times[last])}, -1)

	//Draw("", "", true, nil)
	Draw("", "", false, nil)
}

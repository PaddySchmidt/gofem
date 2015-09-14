// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
)

func Test_platehole01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("platehole01")

	if chk.Verbose {

		var sol PlateHole
		sol.Init([]*fun.Prm{
			&fun.Prm{N: "r", V: 1.0},
			&fun.Prm{N: "E", V: 1e3},
			&fun.Prm{N: "nu", V: 0.3},
			&fun.Prm{N: "qnV", V: 0.0},
			&fun.Prm{N: "qnH", V: 10.0},
		})

		L := 4.0
		npts := 101

		plt.SetForEps(1.2, 455)
		sol.PlotStress(1, L, npts)
		plt.SaveD("/tmp/gofem", "ana_platehole01.eps")
	}
}

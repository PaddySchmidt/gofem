// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gemlab"
	"github.com/cpmech/gosl/io"
)

func main() {

	// dimensions
	lx, ly := 10.0, 3.0

	// define structured mesh data
	var gd gemlab.InData
	gd.Nparts = 4
	gd.Sregs = &gemlab.Sregs{
		Tags: []int{-1},
		Nxs:  []int{10},
		Nys:  []int{3},
		Nzs:  []int{3},
		Points: [][]float64{
			{0, 0}, {lx, 0}, {lx, ly}, {0, ly},
		},
		Conn:  [][]int{{0, 1, 2, 3}},
		Btags: [][]int{{-20, -11, -21, -10}},
	}

	// tag vertices along line (middle vertical line)
	gd.VtagsL = &gemlab.VtagsL{
		Tags: []int{-1},
		Xxa:  [][]float64{{lx / 2.0, 0}},
		Xxb:  [][]float64{{lx / 2.0, ly}},
	}

	fnk := "d2-coarse"
	if err := gemlab.Generate(fnk, &gd); err != nil {
		io.PfRed("%v\n", err.Error())
	}
}

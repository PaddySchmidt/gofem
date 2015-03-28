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
	lx, ly, lz := 10.0, 3.0, 3.0

	// define structured mesh data
	var gd gemlab.InData
	gd.Nparts = 4
	gd.Sregs = &gemlab.Sregs{
		Tags: []int{-1},
		Nxs:  []int{10},
		Nys:  []int{3},
		Nzs:  []int{3},
		Points: [][]float64{
			{0, 0, 0}, {lx, 0, 0}, {lx, ly, 0}, {0, ly, 0},
			{0, 0, lz}, {lx, 0, lz}, {lx, ly, lz}, {0, ly, lz},
		},
		Conn:  [][]int{{0, 1, 2, 3, 4, 5, 6, 7}},
		Btags: [][]int{{-10, -11, -20, -21, -30, -31}},
	}

	fnk := "d3-coarse"
	if err := gemlab.Generate(fnk, &gd); err != nil {
		io.PfRed("%v\n", err.Error())
	}
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_shape01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("shape01")

	r := []float64{0, 0, 0}

	verb := true
	for name, shape := range factory {

		io.Pfyel("--------------------------------- %-6s---------------------------------\n", name)

		// check S
		tol := 1e-17
		if name == "tri10" {
			tol = 1e-14
		}
		CheckShape(tst, shape, tol, verb)

		// check Sf
		tol = 1e-18
		CheckShapeFace(tst, shape, tol, verb)

		// check dSdR
		tol = 1e-14
		if name == "lin5" || name == "lin4" || name == "tri10" || name == "qua12" || name == "qua16" {
			tol = 1e-10
		}
		if name == "tri15" {
			tol = 1e-9
		}
		CheckDSdR(tst, shape, r, tol, verb)

		io.PfGreen("OK\n")
	}
}

func Test_shape02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("shape02")

	xmat := [][]float64{
		{10, 13, 13, 10},
		{8, 8, 9, 9},
	}
	dx, dy := 3.0, 1.0
	dr, ds := 2.0, 2.0
	r := []float64{0, 0, 0}
	shape := factory["qua4"]
	shape.CalcAtIp(xmat, r, true)
	io.Pforan("J = %v\n", shape.J)
	chk.Scalar(tst, "J", 1e-17, shape.J, (dx/dr)*(dy/ds))

	tol := 1e-14
	verb := true
	x := []float64{12.0, 8.5}
	CheckDSdx(tst, shape, xmat, x, tol, verb)
}

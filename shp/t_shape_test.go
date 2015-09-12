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
	chk.PrintTitle("Test shape01")

	verb := true
	for name, shape := range factory {

		io.Pfyel("--------------------------------- %-6s---------------------------------\n", name)

		// check S
		tol := 1e-17
		if name == "tri10" {
			tol = 1e-14
		}
		shape.Check_S(tst, tol, verb)

		// check dSdR
		tol = 1e-14
		if name == "lin5" || name == "lin4" || name == "tri10" || name == "qua12" || name == "qua16" {
			tol = 1e-10
		}
		if name == "tri15" {
			tol = 1e-9
		}
		shape.Check_dSdR(tst, tol, verb)

		// check face vertices
		tol = 1e-18
		shape.Check_Sf(tst, tol, verb)

		io.PfGreen("OK\n")
	}
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/io"
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data
	matOld := io.ArgToString(0, "matOld.mat")
	matNew := io.ArgToString(1, "matNew.mat")
	convSymb := io.ArgToBool(2, true)
	io.Pf("\n%s\n", io.ArgsTable(
		"old material filename", "matOld", matOld,
		"new material filenamen", "matNew", matNew,
		"do convert symbols", "convSymb", convSymb,
	))

	// convert old => new
	inp.MatfileOld2New("", matNew, matOld, convSymb)
	io.Pf("conversion successful\n")
	io.Pfblue2("file <matNew.mat> created\n")
}

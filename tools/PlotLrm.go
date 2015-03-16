// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"flag"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/io"
)

func main() {

	// input data
	simfn := "elast.sim"
	matname := "lrm1"
	pcmax := 30.0
	npts := 101

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		matname = flag.Arg(1)
	}
	if len(flag.Args()) > 2 {
		pcmax = io.Atof(flag.Arg(2))
	}
	if len(flag.Args()) > 3 {
		npts = io.Atoi(flag.Arg(3))
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfn   = %30s // simulation filename\n", simfn)
	io.Pf("  matname = %30s // material name\n", matname)
	io.Pf("  pcmax   = %30v // max pc\n", pcmax)
	io.Pf("  npts    = %30v // number of points\n", npts)
	io.Pf("\n")

	// load simulation
	sim := inp.ReadSim("", simfn, false)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.Mdb.Get(matname)
	if mat == nil {
		io.PfRed("cannot get material\n")
		return
	}
	io.Pforan("mat = %v\n", mat)

	// get and initialise model
	mdl := mreten.GetModel(simfn, matname, mat.Model, false)
	if mdl == nil {
		io.PfRed("cannot allocate model\n")
		return
	}
	mdl.Init(mat.Prms)

	// plot model
	err := mreten.Plot(mdl, 0, 1, pcmax, npts, "'b-'", "", matname)
	if err != nil {
		io.PfRed("plot failed:\n%v\n", err)
		return
	}
	mreten.PlotEnd(true)
}

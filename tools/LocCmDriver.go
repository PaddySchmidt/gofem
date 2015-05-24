// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"flag"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/io"
)

func main() {

	// input data
	dir := "data"
	simfn := "smp-coarse.sim"
	matname := "M.8.4.3-smp"
	pathfn := "path1.pat"

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		dir = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		simfn = flag.Arg(1)
	}
	if len(flag.Args()) > 2 {
		matname = flag.Arg(2)
	}
	if len(flag.Args()) > 3 {
		pathfn = flag.Arg(3)
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}
	if io.FnExt(pathfn) == "" {
		pathfn += ".pat"
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  dir     = %30s // directory with .sim and .pat files\n", dir)
	io.Pf("  simfn   = %30s // simulation filename\n", simfn)
	io.Pf("  matname = %30s // material name\n", matname)
	io.Pf("  pathfn  = %30v // path filename\n", pathfn)
	io.Pf("\n")

	// load simulation
	sim := inp.ReadSim(dir, simfn, "cmd_", false)
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
	mdl, _ := msolid.GetModel(simfn, matname, mat.Model, false)
	if mdl == nil {
		io.PfRed("cannot allocate model\n")
		return
	}
	ndim := 3
	pstress := false
	mdl.Init(ndim, pstress, mat.Prms)
}

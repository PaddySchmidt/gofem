// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/fem"
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
	simfile, _ := io.ArgToFilename(0, "simfile.sim", true)
	zmin := io.ArgToFloat(1, 0.0)
	zmax := io.ArgToFloat(2, 3.0)
	npts := io.ArgToInt(3, 11)
	io.Pf("\n%s\n", io.ArgsTable(
		"simulation filename", "simfile", simfile,
		"min elevation", "zmin", zmin,
		"max elevation", "zmax", zmax,
		"number of points", "npts", npts,
	))

	// sim file
	sim := inp.ReadSim("", simfile, false)
	if sim == nil {
		io.PfRed("cannot read sim file\n")
		return
	}

	// layer
	var lay fem.GeoLayer
	lay.Zmin = zmin
	lay.Zmax = zmax
	lay.Cl = sim.WaterRho0 / sim.WaterBulk
	//if !lay.ReadPorousParameters(sim.Regions[0],
	// TODO

}

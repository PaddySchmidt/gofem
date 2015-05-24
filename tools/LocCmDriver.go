// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"encoding/json"
	"flag"
	"path"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/io"
)

type Input struct {
	Dir     string
	SimFn   string
	MatName string
	PathFn  string
	PlotSet []string
	Eps     bool
}

func (o *Input) PostProcess() {
	if len(o.PlotSet) == 0 {
		o.PlotSet = msolid.PlotSet6
	}
}

func (o Input) String() (l string) {
	l += "\nInput data\n"
	l += "==========\n"
	l += io.Sf("directory with .sim and .pat files : Dir     = %v\n", o.Dir)
	l += io.Sf("simulation filename                : SimFn   = %v\n", o.SimFn)
	l += io.Sf("material name                      : MatName = %v\n", o.MatName)
	l += io.Sf("path filename                      : PathFn  = %v\n", o.PathFn)
	l += io.Sf("plot set                           : PlotSet = %q\n", o.PlotSet)
	l += io.Sf("generate .eps instead of .png      : Eps     = %v\n", o.Eps)
	l += "\n"
	return
}

func main() {

	// input data file
	inpfn := "data/loccmdrv1.inp"
	flag.Parse()
	if len(flag.Args()) > 0 {
		inpfn = flag.Arg(0)
	}
	if io.FnExt(inpfn) == "" {
		inpfn += ".inp"
	}

	// read and parse input data
	var in Input
	b, err := io.ReadFile(inpfn)
	if err != nil {
		io.PfRed("cannot read %s\n", inpfn)
		return
	}
	err = json.Unmarshal(b, &in)
	if err != nil {
		io.PfRed("cannot parse %s\n", inpfn)
		return
	}
	in.PostProcess()

	// print input data
	io.Pf("%v\n", in)

	// load simulation
	sim := inp.ReadSim(in.Dir, in.SimFn, "cmd_", false)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.Mdb.Get(in.MatName)
	if mat == nil {
		io.PfRed("cannot get material\n")
		return
	}
	//io.Pfcyan("mat = %v\n", mat)

	// get and initialise model
	mdl, _ := msolid.GetModel(in.SimFn, in.MatName, mat.Model, false)
	if mdl == nil {
		io.PfRed("cannot allocate model\n")
		return
	}
	ndim := 3
	pstress := false
	mdl.Init(ndim, pstress, mat.Prms)
	//io.Pforan("mdl = %v\n", mdl)

	// load path
	var pth msolid.Path
	err = pth.ReadJson(ndim, path.Join(in.Dir, in.PathFn))
	if err != nil {
		io.PfRed("cannot read path file %v\n", err)
		return
	}
	//io.PfYel("pth = %v\n", pth)

	// driver
	var drv msolid.Driver
	drv.InitWithModel(ndim, mdl)

	// run
	err = drv.Run(&pth)
	if err != nil {
		io.PfRed("driver: Run failed: %v\n", err)
		return
	}

	for _, sta := range drv.Res {
		io.Pforan("sta = %v\n", sta.Sig)
	}
	io.Pfcyan("eps = %v\n", drv.Eps)

	// plot
	var plr msolid.Plotter
	plr.SetFig(false, in.Eps, 1.0, 400, "/tmp", "cmd_"+in.SimFn)
	//plr.SetModel(mdl)
	plr.Plot(in.PlotSet, drv.Res, drv.Eps, true, true)
}

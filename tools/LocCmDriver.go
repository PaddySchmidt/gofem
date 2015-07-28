// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"encoding/json"
	"path"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/tsr"
)

type Input struct {
	Dir     string
	SimFn   string
	MatName string
	PathFn  string
	PlotSet []string
	FigEps  bool
	FigProp float64
	FigWid  float64

	// derived
	inpfn string
}

func (o *Input) PostProcess() {
	if len(o.PlotSet) == 0 {
		o.PlotSet = msolid.PlotSet6
	}
	if o.FigProp < 0.1 {
		o.FigProp = 1.0
	}
	if o.FigWid < 10 {
		o.FigWid = 400
	}
}

func (o Input) String() (l string) {
	l = io.ArgsTable(
		"input filename", "inpfn", o.inpfn,
		"directory with .sim and .pat files", "Dir", o.Dir,
		"simulation filename", "SimFn", o.SimFn,
		"material name", "MatName", o.MatName,
		"path filename", "PathFn", o.PathFn,
		"plot set", "PlotSet", io.Sf("%v", o.PlotSet),
		"fig: generate .eps instead of .png", "FigEps", o.FigEps,
		"fig: proportion of figure", "FigProp", o.FigProp,
		"fig: width  of figure", "FigWid", o.FigWid,
	)
	return
}

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data file
	var in Input
	in.inpfn, _ = io.ArgToFilename(0, "data/loccmdrv1", ".inp", true)

	// read and parse input data
	b, err := io.ReadFile(in.inpfn)
	if err != nil {
		io.PfRed("cannot read %s\n", in.inpfn)
		return
	}
	err = json.Unmarshal(b, &in)
	if err != nil {
		io.PfRed("cannot parse %s\n", in.inpfn)
		return
	}
	in.PostProcess()

	// print input table
	io.Pf("%v\n", in)

	// load simulation
	sim := inp.ReadSim(in.Dir+"/"+in.SimFn, "", false, 0)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.MatParams.Get(in.MatName)
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
		io.Pfred("driver: Run failed: %v\n", err)
	}

	// plot
	//if false {
	if true {
		var plr msolid.Plotter
		plr.SetFig(false, in.FigEps, in.FigProp, in.FigWid, "/tmp", "cmd_"+in.SimFn)
		var epm msolid.EPmodel
		if m, ok := mdl.(msolid.EPmodel); ok {
			plr.SetModel(m)
			epm = m
		}
		if epm != nil {
			//plr.Phi = epm.Get_phi()
			b := epm.Get_bsmp()
			epm.Set_bsmp(0)
			plr.YsClr0 = "magenta"
			plr.Plot(in.PlotSet, drv.Res, nil, true, false)
			epm.Set_bsmp(b)
		}
		plr.YsClr0 = "green"
		plr.Plot(in.PlotSet, drv.Res, drv.Eps, false, true)
	}

	// plot ys
	if false {
		//if true {
		plt.Reset()
		m := mdl.(*msolid.SmpInvs)
		φ := m.Get_phi()
		σcCte := 10.0
		M := tsr.Phi2M(φ, "oct")
		rmin, rmax := 0.0, 1.28*M*σcCte
		nr, nα := 31, 81
		//nr,   nα   := 31, 1001
		npolarc := true
		simplec := false
		only0 := false
		grads := false
		showpts := false
		ferr := 10.0
		tsr.PlotOct("fig_isofun02.png", σcCte, rmin, rmax, nr, nα, φ, m.Isof.Fa, m.Isof.Ga,
			npolarc, simplec, only0, grads, showpts, true, true, ferr)
	}
}

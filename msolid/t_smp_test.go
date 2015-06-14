// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

func Test_smp01(tst *testing.T) {

	defer func() {
		if err := recover(); err != nil {
			io.Pfred("error = %v\n", err)
		}
	}()

	//verbose()
	chk.PrintTitle("smp01")

	E, ν := 1500.0, 0.25
	K := Calc_K_from_Enu(E, ν)
	G := Calc_G_from_Enu(E, ν)
	io.Pforan("K = %v\n", K)
	io.Pforan("G = %v\n", G)

	// allocate driver
	ndim, pstress := 2, false
	simfnk, modelname := "test", "smp"
	var drv Driver
	err := drv.Init(simfnk, modelname, ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "c", V: 1},
		&fun.Prm{N: "phi", V: 20},
		&fun.Prm{N: "a", V: -1},
		&fun.Prm{N: "b", V: 0.0},
		&fun.Prm{N: "bet", V: 1},
		&fun.Prm{N: "eps", V: 1e-3},
		&fun.Prm{N: "le", V: 1},
		&fun.Prm{N: "pr", V: 1.0},
		&fun.Prm{N: "G0", V: G},
		&fun.Prm{N: "K0", V: K},
		&fun.Prm{N: "p0", V: 0.0},
		&fun.Prm{N: "ev0", V: 0.0},
		&fun.Prm{N: "rtyp", V: 1.0},
		&fun.Prm{N: "r", V: 1.0},
		&fun.Prm{N: "pe", V: 10.0},
	})

	// set flags
	drv.CheckD = true
	//drv.CheckD = false
	drv.TolD = 1e-3
	drv.VerD = io.Verbose // verbose
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// model
	smp := drv.model.(*SmpInvs)

	// path
	p0 := 0.0
	DP := []float64{2, -1, -5}
	DQ := []float64{6, 0, 0}
	//DP := []float64{-4}
	//DQ := []float64{8}
	nincs := 1
	niout := 1
	noise := 0.0
	var pth Path
	//if true {
	if false {
		err = pth.SetPQstrain(ndim, nincs, niout, K, G, p0, DP, DQ, noise)
		if err != nil {
			tst.Errorf("test failed: %v\n", err)
			return
		}
	} else {
		pth.Sx = []float64{-1}
		pth.Sy = []float64{-2}
		pth.Sz = []float64{-1}
		pth.Ex = []float64{0, 0, 0.001, -0.004}
		pth.Ey = []float64{0, -0.006, 0, -0.002}
		pth.Ez = []float64{0, 0, 0.001, -0.004}
		//pth.Ex = []float64{0, -0.0033333333333333335, -0.0028333333333333335}
		//pth.Ey = []float64{0, -0.0033333333333333335, -0.0028333333333333335}
		//pth.Ez = []float64{0, -0.0033333333333333335, -0.005333333333333334}
		pth.UseS = []int{0, 0}
		pth.UseE = []int{0, 1}
		pth.Init(ndim)
	}

	// run
	err = drv.Run(&pth)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// plot
	//if true {
	if false {
		var plr Plotter
		prop := 2.0
		plr.SetFig(false, false, prop, 400, "/tmp", "test_smp01")
		plr.SetModel(smp)
		//n := len(drv.Res)
		plr.Pmin = -smp.HE.pt
		plr.Pmax = 5
		plr.PqLims = []float64{plr.Pmin, plr.Pmax, -6, 6}
		plr.UsePmin = true
		plr.UsePmax = true
		plr.PreCor = drv.PreCor
		//plr.Plot(PlotSet8, drv.Res, drv.Eps, true, true)
		plr.Plot([]string{"p,q,ys", "p,ev", "s3,s1,ys", "oct,ys"}, drv.Res, drv.Eps, true, true)
		//plr.Plot([]string{"log(p),ev", "p,q,ys"}, drv.Res, drv.Eps, true, true)
	}
}

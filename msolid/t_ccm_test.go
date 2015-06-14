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

func Test_ccm01(tst *testing.T) {

	defer func() {
		if err := recover(); err != nil {
			io.Pfred("error = %v\n", err)
		}
	}()

	verbose()
	chk.PrintTitle("ccm01")

	E, ν := 1500.0, 0.25
	K := Calc_K_from_Enu(E, ν)
	G := Calc_G_from_Enu(E, ν)
	io.Pforan("K = %v\n", K)
	io.Pforan("G = %v\n", G)

	// allocate driver
	pr := 1.0
	ndim, pstress := 2, false
	simfnk, modelname := "test", "ccm"
	var drv Driver
	err := drv.Init(simfnk, modelname, ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "phi", V: 25},
		&fun.Prm{N: "Mfix", V: 1},
		&fun.Prm{N: "c", V: 1},
		&fun.Prm{N: "lam", V: 0.1},
		&fun.Prm{N: "ocr", V: 1},
		&fun.Prm{N: "kap", V: 0.05},
		&fun.Prm{N: "kapb", V: 0.01},
		&fun.Prm{N: "G0", V: G},
		&fun.Prm{N: "pr", V: pr},
		&fun.Prm{N: "le", V: 0},
		&fun.Prm{N: "K0", V: K},
	})
	drv.CheckD = true
	//drv.CheckD = false
	drv.TolD = 1e-4
	drv.VerD = io.Verbose // verbose
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// model
	ccm := drv.model.(*CamClayMod)

	// path
	p0 := 0.0
	DP := []float64{10, 1}
	DQ := []float64{0, 3}
	nincs := 1
	niout := 1
	noise := 0.0
	var pth Path
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
		pth.Ex = []float64{0, 0, 0.001, -0.004}  //, -0.005}
		pth.Ey = []float64{0, -0.005, 0, -0.002} //, -0.005}
		pth.Ez = []float64{0, 0, 0.001, -0.004}  //, -0.005}
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
		plr.Pr = pr
		prop := 2.0
		plr.SetFig(false, false, prop, 400, "/tmp", "test_ccm01")
		plr.SetModel(ccm)
		n := len(drv.Res)
		plr.Pmin = -ccm.HE.pt
		plr.Pmax = drv.Res[n-1].Alp[0]
		qmax := ccm.CS.Mcs * (plr.Pmax - plr.Pmin)
		plr.PqLims = []float64{plr.Pmin, plr.Pmax, -qmax, qmax}
		plr.UsePmin = true
		plr.UsePmax = true
		plr.PreCor = drv.PreCor
		//plr.Plot(PlotSet8, drv.Res, drv.Eps, true, true)
		plr.Plot([]string{"log(p),ev", "p,q,ys"}, drv.Res, drv.Eps, true, true)
	}
}

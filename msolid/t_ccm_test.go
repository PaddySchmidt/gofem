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
		&fun.Prm{N: "kapb", V: 0},
		&fun.Prm{N: "G0", V: G},
		&fun.Prm{N: "pr", V: 1.0},
		&fun.Prm{N: "p0", V: 0.0},
		&fun.Prm{N: "ev0", V: 0.0},
		&fun.Prm{N: "le", V: 0},
		&fun.Prm{N: "K0", V: K},
	})
	drv.CheckD = true
	drv.TolD = 1e-6
	drv.VerD = true // verbose
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
	err = pth.SetPQstrain(ndim, nincs, niout, K, G, p0, DP, DQ, noise)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// run
	err = drv.Run(&pth)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}
	io.PfYel("\nresults\n")
	for i, res := range drv.Res {
		io.Pfyel("ε=%v σ=%v α0=%v\n", drv.Eps[i], res.Sig, res.Alp[0])
	}

	// plot
	//if true {
	if false {
		var plr Plotter
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

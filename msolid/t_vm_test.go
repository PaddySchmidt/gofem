// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

func Test_vm01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("vm01")

	// allocate driver
	ndim, pstress := 2, false
	simfnk, modelname := "test", "vm"
	var drv Driver
	drv.Init(simfnk, modelname, ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "K", V: 1.5},
		&fun.Prm{N: "G", V: 1},
		&fun.Prm{N: "qy0", V: 2},
		&fun.Prm{N: "H", V: 0.5},
	})
	drv.CheckD = true

	// vm model
	vm := drv.model.(*VonMises)

	// path
	p0 := 0.0
	Δp := 3.0
	Δq := vm.qy0
	ϵ := 1e-3
	DP := []float64{Δp + ϵ, 3, 2, 1, 0}
	DQ := []float64{Δq + ϵ, 4, 2, 1, 3}
	nincs := 1
	niout := 1
	noise := 0.0
	var pth Path
	pth.SetPQstrain(ndim, nincs, niout, vm.K, vm.G, p0, DP, DQ, noise)

	// run
	err := drv.Run(&pth)
	if err != nil {
		utl.Panic("%v", err.Error())
	}

	// plot
	//if DPsaveFig {
	//}
}
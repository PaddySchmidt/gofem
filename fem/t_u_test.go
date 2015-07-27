// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

func Test_sigini01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sigini01. zero displacements. initial stresses")

	// fem
	fem := NewFEM("data/sigini01.sim", "", true, false, false, chk.Verbose)

	// set stage
	err := fem.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectors
	err = fem.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := fem.Domains[0]

	// check displacements
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
		chk.Vector(tst, "u", tolu, u, nil)
	}

	// analytical solution
	qnV, qnH := -100.0, -50.0
	ν := 0.25
	σx, σy := qnH, qnV
	σz := ν * (σx + σy)
	σref := []float64{σx, σy, σz, 0}

	// check stresses
	e := dom.Elems[0].(*ElemU)
	tols := 1e-13
	for idx, _ := range e.IpsElem {
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		chk.Vector(tst, "σ", tols, σ, σref)
	}
}

func Test_sigini02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sigini02. initial stresses. run simulation")

	// fem
	fem := NewFEM("data/sigini02.sim", "", true, false, false, chk.Verbose)

	// run simulation
	err := fem.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := fem.Domains[0]

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH0", V: -20},
		&fun.Prm{N: "qnV0", V: -20},
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	t := dom.Sol.T
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	e := dom.Elems[0].(*ElemU)
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}

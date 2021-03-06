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
	"github.com/cpmech/gosl/la"
)

func Test_sigini01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sigini01. zero displacements. initial stresses")

	// fem
	analysis := NewFEM("data/sigini01.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectors
	err = analysis.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

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
	analysis := NewFEM("data/sigini02.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

	// external force
	e := dom.Elems[0].(*ElemU)
	io.PfYel("fex = %v\n", e.fex)
	io.PfYel("fey = %v\n", e.fey)
	la.PrintMat("K", e.K, "%10.2f", false)

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
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Cell.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}

func Test_square01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("square01. ini stress free square")

	// fem
	analysis := NewFEM("data/square01.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

	// external force
	e := dom.Elems[0].(*ElemU)
	io.PfYel("fex = %v\n", e.fex)
	io.PfYel("fey = %v\n", e.fey)
	la.PrintMat("K", e.K, "%10.2f", false)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
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
		io.Pfyel("u = %v\n", u)
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Cell.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}

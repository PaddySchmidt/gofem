// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

func Test_nurbs01(tst *testing.T) {

	/*  4 (1,2)             (2,2) 6
	    5  2@o--------------o@3   7
	         |              |
	         |              |      @     -- control point
	         |              |      o     -- node
	         |              |      (a,b) -- span
	         |              |
	         |              |
	    0  0@o--------------o@1   2
	    1 (1,1)             (2,1) 3
	*/

	//verbose()
	chk.PrintTitle("nurb01. square with initial stress")

	// fem
	analysis := NewFEM("data/nurbs01.sim", "", true, false, false, false, chk.Verbose, 0)

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

	// draw NURBS
	if false {
		nurbs := dom.Msh.Cells[0].Shp.Nurbs
		gm.PlotNurbs("/tmp/gofem", "test_nurbs01", nurbs)
	}

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 4)
	chk.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
		chk.StrAssert(nod.Dofs[0].Key, "ux")
		chk.StrAssert(nod.Dofs[1].Key, "uy")
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "eqs", eqs, utl.IntRange(4*2))
	chk.Ints(tst, "nids", nids, []int{0, 1, 2, 3})

	// check Umap
	Umaps := [][]int{
		{0, 1, 2, 3, 4, 5, 6, 7},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemU)
		io.Pfpink("%2d : Umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "Umap", e.Umap, Umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), 4)
	var ct_ux_eqs []int // equations with ux prescribed [sorted]
	var ct_uy_eqs []int // equations with uy prescribed [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pfgrey("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "ux":
			ct_ux_eqs = append(ct_ux_eqs, eq)
		case "uy":
			ct_uy_eqs = append(ct_uy_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	chk.Ints(tst, "equations with ux prescribed", ct_ux_eqs, []int{0, 4})
	chk.Ints(tst, "equations with uy prescribed", ct_uy_eqs, []int{1, 3})

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

func Test_nurbs02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs02. square with initial stress. run")

	// fem
	analysis := NewFEM("data/nurbs02.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

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
		io.Pfyel("u = %v\n", u)
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	if false {
		// check stresses
		tols := 1e-13
		for idx, ip := range e.IpsElem {
			x := e.Cell.Shp.IpRealCoords(e.X, ip)
			σ := e.States[idx].Sig
			io.Pforan("σ = %v\n", σ)
			sol.CheckStress(tst, t, σ, x, tols)
		}
	}
}

func Test_nurbs03(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs03. ini stress free square")

	// fem
	analysis := NewFEM("data/nurbs03.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

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
}

func Test_nurbs04(tst *testing.T) {

	verbose()
	chk.PrintTitle("nurbs04. perforated disk")

	// fem
	analysis := NewFEM("data/nurbs04.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// draw NURBS
	if false {
		dom := analysis.Domains[0]
		nurbs := dom.Msh.Cells[0].Shp.Nurbs
		gm.PlotNurbs("/tmp/gofem", "test_nurbs04", nurbs)
	}
}

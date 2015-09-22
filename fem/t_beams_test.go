// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_beam01a(tst *testing.T) {

	//verbose()
	chk.PrintTitle("beam01a")

	// fem
	analysis := NewFEM("data/beam01.sim", "", true, false, false, false, chk.Verbose, 0)

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

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 2)
	chk.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 3)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 1})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5})

	// check solution arrays
	ny := 6
	nλ := 3
	nyb := ny + nλ
	chk.IntAssert(len(dom.Sol.Y), ny)
	chk.IntAssert(len(dom.Sol.Dydt), 0)
	chk.IntAssert(len(dom.Sol.D2ydt2), 0)
	chk.IntAssert(len(dom.Sol.Psi), 0)
	chk.IntAssert(len(dom.Sol.Zet), 0)
	chk.IntAssert(len(dom.Sol.Chi), 0)
	chk.IntAssert(len(dom.Sol.L), nλ)
	chk.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	chk.IntAssert(len(dom.Fb), nyb)
	chk.IntAssert(len(dom.Wb), nyb)

	// check umap
	e := dom.Elems[0].(*Beam)
	io.Pforan("e = %v\n", e.Umap)
	chk.Ints(tst, "umap", e.Umap, []int{0, 1, 2, 3, 4, 5})

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_ux_eqs []int // constrained ux equations [sorted]
	var ct_uy_eqs []int // constrained uy equations [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pforan("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "ux":
			ct_ux_eqs = append(ct_ux_eqs, eq)
		case "uy":
			ct_uy_eqs = append(ct_uy_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
			return
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{0})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 4})
}

func Test_beam01b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("beam01b. simply supported")

	// start simulation
	analysis := NewFEM("data/beam01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	dom := analysis.Domains[0]
	ele := dom.Elems[0].(*Beam)
	_, M := ele.CalcVandM(dom.Sol, 0.5, 1)
	qn, L := 15.0, 1.0
	Mcentre := qn * L * L / 8.0
	io.Pforan("M = %v (%v)\n", M, Mcentre)
	chk.Scalar(tst, "M @ centre", 1e-17, M[0], Mcentre)

	// check moment using OutIpsData
	idx_centre := 5 // considering 11 stations
	dat := ele.OutIpsData()
	res := dat[idx_centre].Calc(dom.Sol)
	io.Pfcyan("M @ centre (OutIpsData) = %v\n", res["M"])
	chk.Scalar(tst, "M @ centre (OutIpsData)", 1e-17, res["M"], Mcentre)
}

func Test_beam02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("beam02. cantilever")

	// start simulation
	analysis := NewFEM("data/beam02.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	dom := analysis.Domains[0]
	ele := dom.Elems[0].(*Beam)
	_, M := ele.CalcVandM(dom.Sol, 0, 1)
	qn, L := 15.0, 1.0
	Mleft := -qn * L * L / 2.0
	io.Pforan("M = %v (%v)\n", M, Mleft)
	chk.Scalar(tst, "M @ left", 1e-15, M[0], Mleft)
}

func Test_beam03(tst *testing.T) {

	//verbose()
	chk.PrintTitle("beam03. small frame")

	// start simulation
	analysis := NewFEM("data/beam03.sim", "", true, true, false, false, chk.Verbose, 0)

	// set stage and run
	set_and_run := func(stgidx int) {
		err := analysis.SetStage(stgidx)
		if err != nil {
			tst.Errorf("SetStage failed:\n%v", err)
			return
		}
		err = analysis.SolveOneStage(stgidx, true)
		if err != nil {
			tst.Error("SolveOneStage failed:\n%v", err)
			return
		}
	}
	set_and_run(0)

	// domain
	dom := analysis.Domains[0]

	// message
	io.Pf("\nebc = %v\n", dom.EssenBcs.List(0))
	io.Pf("fbc = %v\n\n", dom.PtNatBcs.List(0))
	for _, nod := range dom.Nodes {
		io.Pf("node # %2d ux=%23.15e uy=%23.15e\n", nod.Vert.Id, dom.Sol.Y[nod.GetEq("ux")], dom.Sol.Y[nod.GetEq("uy")])
	}
	io.Pf("\n")

	// define function to check bending moment
	check_M := func(beamId int, s, Mref, tol float64) {
		ele := dom.Cid2elem[beamId].(*Beam)
		_, M := ele.CalcVandM(dom.Sol, s, 1)
		chk.Scalar(tst, io.Sf("Beam %d: M(s=%g) = %.6f", ele.Id(), s, M[0]), tol, M[0], Mref)
	}

	// check
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 34, 1e-13)
	check_M(1, 0, -20.4, 1e-13)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, -54.4, 1e-13)

	nstations, withtext, numfmt, tol, coef := 11, true, "", 1e-10, 0.2
	if chk.Verbose {
		plt.SetForPng(1, 600, 150)
		PlotAllBendingMoments(dom, nstations, withtext, numfmt, tol, coef)
		plt.SaveD("/tmp/gofem", "test_beam03_prob1.png")
	}

	// problem # 2
	set_and_run(1)
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 34, 1e-13)
	check_M(1, 0, -20.4, 1e-13)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		plt.SetForPng(1, 600, 150)
		PlotAllBendingMoments(dom, nstations, withtext, numfmt, tol, coef)
		plt.SaveD("/tmp/gofem", "test_beam03_prob2.png")
	}

	// problem # 3
	set_and_run(2)
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 20, 1e-12)
	check_M(1, 0, 20, 1e-12)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		plt.SetForPng(1, 600, 150)
		PlotAllBendingMoments(dom, nstations, withtext, numfmt, tol, coef)
		plt.SaveD("/tmp/gofem", "test_beam03_prob3.png")
	}

	// problem # 4
	set_and_run(3)
	e0M := func(x float64) float64 { return 5.0752*x - 0.9984*x*x/3.0 - (16.0-x)*0.0624*x*x/6.0 }
	e1M := func(x float64) float64 { return 1.7472*(10.0-x) - 0.6*0.0624*math.Pow(10.0-x, 3.0)/6.0 }
	check_M(0, 0, e0M(0), 1e-13)
	check_M(0, 0.5, e0M(5), 1e-13)
	check_M(0, 1, e0M(10), 1e-12)
	check_M(1, 0, e1M(0), 1e-12)
	check_M(1, 0.5, e1M(5), 1e-12)
	check_M(1, 1, e1M(10), 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		plt.SetForPng(1, 600, 150)
		PlotAllBendingMoments(dom, nstations, withtext, numfmt, tol, coef)
		plt.SaveD("/tmp/gofem", "test_beam03_prob4.png")
	}
}

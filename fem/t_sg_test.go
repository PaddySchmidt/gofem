// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_sg52a(tst *testing.T) {

	/* Smith & Griffths (5th ed) Figure 5.2 p173
	 *
	 *          0.25       0.5      0.25 kN/m
	 *            ↓         ↓         ↓
	 *    ---    ▷0---------1---------2
	 *     |      |       ,'|       ,'|   E = 1e6 kN/m²
	 *     |      |  0  ,'  |  2  ,'  |   ν = 0.3
	 *     |      |   ,'    |   ,'    |
	 *            | ,'   1  | ,'  3   |   connectivity:
	 *    1 m    ▷3'--------4'--------5     0 : 1 0 3
	 *            |       ,'|       ,'|     1 : 3 4 1
	 *     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
	 *     |      |   ,'    |   ,'    |     3 : 4 5 2
	 *     |      | ,'   5  | ,'   7  |     4 : 4 3 6
	 *    ---    ▷6'--------7'--------8     5 : 6 7 4
	 *            △         △         △     6 : 5 4 7
	 *                                      7 : 7 8 5
	 *            |------- 1 m -------|
	 */

	//verbose()
	chk.PrintTitle("sg52a")

	// start simulation
	analysis := NewFEM("data/sg52.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectros
	err = analysis.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 9)
	chk.IntAssert(len(dom.Elems), 8)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{1, 0, 3, 4, 2, 5, 6, 7, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17})

	// check solution arrays
	ny := 9 * 2
	nλ := 6
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
	umaps := [][]int{
		{0, 1, 2, 3, 4, 5},
		{4, 5, 6, 7, 0, 1},
		{8, 9, 0, 1, 6, 7},
		{6, 7, 10, 11, 8, 9},
		{6, 7, 4, 5, 12, 13},
		{12, 13, 14, 15, 6, 7},
		{10, 11, 6, 7, 14, 15},
		{14, 15, 16, 17, 10, 11},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemU)
		io.Pforan("e%d.umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "umap", e.Umap, umaps[i])
	}

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
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{2, 4, 12})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{13, 15, 17})

	// point loads
	chk.IntAssert(len(dom.PtNatBcs.Bcs), 3)
	chk.StrAssert(dom.PtNatBcs.Bcs[0].Key, "fuy")
	chk.StrAssert(dom.PtNatBcs.Bcs[1].Key, "fuy")
	chk.StrAssert(dom.PtNatBcs.Bcs[2].Key, "fuy")
	chk.IntAssert(dom.PtNatBcs.Bcs[0].Eq, 3)
	chk.IntAssert(dom.PtNatBcs.Bcs[1].Eq, 1)
	chk.IntAssert(dom.PtNatBcs.Bcs[2].Eq, 9)
}

func Test_sg52b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg52b")

	// run simulation
	analysis := NewFEM("data/sg52.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-9
	tolu := 1e-17
	tols := 1.56e-15
	TestingCompareResultsU(tst, "data/sg52.sim", "cmp/sg52.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg57(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg57")

	// run simulation
	analysis := NewFEM("data/sg57.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.35
	tolu := 2e-9
	tols := 0.0002
	TestingCompareResultsU(tst, "data/sg57.sim", "cmp/sg57.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg511(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg511")

	// run simulation
	analysis := NewFEM("data/sg511.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.1
	tolu := 3e-14
	tols := 1.56e-7
	TestingCompareResultsU(tst, "data/sg511.sim", "cmp/sg511.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg515(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg515")

	// run simulation
	analysis := NewFEM("data/sg515.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.15
	tolu := 3e-13
	tols := 3e-8
	TestingCompareResultsU(tst, "data/sg515.sim", "cmp/sg515.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg517(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg517")

	// run simulation
	analysis := NewFEM("data/sg517.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.0036
	tolu := 1e-6
	tols := 1e-4
	TestingCompareResultsU(tst, "data/sg517.sim", "cmp/sg517.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg524(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg524")

	// run simulation
	analysis := NewFEM("data/sg524.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-8
	tols := 1e-7
	TestingCompareResultsU(tst, "data/sg524.sim", "cmp/sg524.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_sg530(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sg530")

	// run simulation
	analysis := NewFEM("data/sg530.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-17
	tols := 1e-15
	TestingCompareResultsU(tst, "data/sg530.sim", "cmp/sg530.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

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
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func test_contact01a(tst *testing.T) {

	/*     Nodes                           Equations:
	 *                   7                                  33,34
	 *     6 o-----o-----o-----o-----o 8       35 o-----o-----o-----o-----o 45,46,47
	 *       |    13     |    14     |         36 |   39,40   |   51,52   |
	 *       |           |           |            |         37|           | 48
	 *    18 o    [2]    o19  [3]    o 20      41 o     o   38o     o     o 49
	 *       |    23     |    24     |         42 |   43,44   |   53,54   | 50
	 *       |           |           |            |          4|           |
	 *     3 o-----o-----o-----o-----o 5        6 o-----o-----o-----o-----o 21,22,23
	 *       |    11    4|    12     |          7 |   12,13  5|   29,30   |
	 *       |           |           |            |           |           | 26
	 *    15 o    [0]    o16  [1]    o 17   14,15 o     o   10o     o     o 27
	 *       |    21     |    22     |            |   16,17 11|   31,32   | 28
	 *       |           |           |            |           |           |
	 *     0 o-----o-----o-----o-----o 2        0 o-----o-----o-----o-----o 18
	 *             9     1     10               1       8     2   24,25     19
	 *                                                  9     3             20
	 */

	//verbose()
	chk.PrintTitle("contact01a")

	// start simulation
	analysis := NewFEM("data/contact01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 25)
	chk.IntAssert(len(dom.Elems), 4)

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 1, 4, 3, 9, 16, 11, 15, 21, 2, 5, 10, 17, 12, 22, 7, 6, 19, 13, 18, 23, 8, 20, 14, 24})
	chk.Ints(tst, "eqs", eqs, utl.IntRange(55))

	// vertices with "qbn"
	contactverts := map[int]bool{2: true, 17: true, 5: true, 20: true, 8: true}

	// check dofs
	var contacteqs []int
	for _, nod := range dom.Nodes {
		if contactverts[nod.Vert.Id] {
			chk.IntAssert(len(nod.Dofs), 3)
			contacteqs = append(contacteqs, nod.Dofs[2].Eq)
		} else {
			chk.IntAssert(len(nod.Dofs), 2)
		}
	}
	sort.Ints(contacteqs)
	io.Pforan("contacteqs = %v\n", contacteqs)
	chk.Ints(tst, "contacteqs", contacteqs, []int{20, 23, 28, 47, 50})

	// check Qmap
	e1 := dom.Elems[1].(*ElemU)
	chk.Ints(tst, "e1.Qmap", e1.Qmap, []int{20, 23, 28})
	e3 := dom.Elems[3].(*ElemU)
	chk.Ints(tst, "e3.Qmap", e3.Qmap, []int{23, 47, 50})
}

func Test_contact01b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("contact01b")

	// start simulation
	analysis := NewFEM("data/contact01.sim", "", true, true, false, false, chk.Verbose, 0)

	// for debugging Kb
	//if true {
	if false {
		u_DebugKb(analysis, &testKb{
			tst: tst, eid: 3, tol: 1e-7, verb: chk.Verbose,
			ni: -1, nj: -1, itmin: 1, itmax: 1, tmin: 0.2, tmax: -1,
		})
	}

	// run simulation
	err := analysis.Run()
	if err != nil {
		io.PfRed("Run failed:\n%v", err)
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	//if true {
	if false {

		// domain
		dom := analysis.Domains[0]

		// solution
		var sol ana.CteStressPstrain
		sol.Init(fun.Prms{
			&fun.Prm{N: "qnH", V: 0},
			&fun.Prm{N: "qnV", V: -100},
		})

		// check displacements
		t := dom.Sol.T
		tolu := 1e-16
		for _, n := range dom.Nodes {
			eqx := n.GetEq("ux")
			eqy := n.GetEq("uy")
			u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
			io.Pfgreen("u = %v\n", u)
			sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
		}

		// check stresses
		e := dom.Elems[3].(*ElemU)
		tols := 1e-13
		for idx, ip := range e.IpsElem {
			x := e.Cell.Shp.IpRealCoords(e.X, ip)
			σ := e.States[idx].Sig
			io.Pforan("σ = %v\n", σ)
			sol.CheckStress(tst, t, σ, x, tols)
		}
	}
}

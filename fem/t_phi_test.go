// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_phi01(tst *testing.T) {

	/*     Nodes                           Equations: TODO
	 *                 7                                7
	 *     6 o----o----o----o----o 8        6 o----o----o----o----o 8
	 *       |   13    |   14    |            |   13    |   14    |
	 *       |         |         |            |         |         |
	 *    18 o   [2]   o19 [3]   o 20      18 o    o    o19  o    o 20
	 *       |   23    |   24    |            |   23    |   24    |
	 *       |         |         |            |         |         |
	 *     3 o----o----o----o----o 5        3 o----o----o----o----o 5
	 *       |   11   4|   12    |            |    6   2|   12    |
	 *       |         |         |            |         |         |
	 *    15 o   [0]   o16 [1]   o 17       7 o    o    o5   o    o 17
	 *       |   21    |   22    |            |    8    |   22    |
	 *       |         |         |            |         |         |
	 *     0 o----o----o----o----o 2        0 o----o----o----o----o 2
	 *            9    1    10                     4    1    10
	 */

	//verbose()
	chk.PrintTitle("phi01a")

	// make sure to flush log
	defer End()

	// start simulation
	if !Start("data/phi01.sim", true, chk.Verbose) {
		tst.Errorf("test failed\n")
		return
	}

	// domain
	distr := false
	dom := NewDomain(Global.Sim.Regions[0], distr)
	if dom == nil {
		tst.Errorf("test failed\n")
		return
	}

	// set stage
	if !dom.SetStage(0, Global.Sim.Stages[0], distr) {
		tst.Errorf("test failed\n")
		return
	}

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 27)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 1)
		chk.StrAssert(nod.Dofs[0].Key, "h")
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{
		0, 1, 3, 2,
	})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3})

	/*
		// check pmap
		pmaps := [][]int{
			{0, 1, 2, 3, 4, 5, 6, 7, 8},
			{3, 2, 9, 10, 6, 11, 12, 13, 14},
			{10, 9, 15, 16, 12, 17, 18, 19, 20},
			{16, 15, 21, 22, 18, 23, 24, 25, 26},
		}
		for i, ele := range dom.Elems {
			e := ele.(*ElemP)
			io.Pforan("e%d.pmap = %v\n", e.Id(), e.Pmap)
			chk.Ints(tst, "pmap", e.Pmap, pmaps[i])
		}

		// constraints
		chk.IntAssert(len(dom.EssenBcs.Bcs), 3)
		var ct_pl_eqs []int // equations with pl prescribed [sorted]
		for _, c := range dom.EssenBcs.Bcs {
			chk.IntAssert(len(c.Eqs), 1)
			eq := c.Eqs[0]
			io.Pforan("key=%v eq=%v\n", c.Key, eq)
			switch c.Key {
			case "pl":
				ct_pl_eqs = append(ct_pl_eqs, eq)
			default:
				tst.Errorf("key %s is incorrect", c.Key)
			}
		}
		sort.Ints(ct_pl_eqs)
		chk.Ints(tst, "equations with pl prescribed", ct_pl_eqs, []int{0, 1, 4})

		// initial values @ nodes
		io.Pforan("initial values @ nodes\n")
		for _, nod := range dom.Nodes {
			z := nod.Vert.C[1]
			eq := nod.Dofs[0].Eq
			pl := dom.Sol.Y[eq]
			plC, _, _ := Global.HydroSt.Calc(z)
			chk.Scalar(tst, io.Sf("nod %3d : pl(@ %4g)= %6g", nod.Vert.Id, z, pl), 1e-17, pl, plC)
		}

		// intial values @ integration points
		io.Pforan("initial values @ integration points\n")
		for _, ele := range dom.Elems {
			e := ele.(*ElemP)
			for idx, ip := range e.IpsElem {
				s := e.States[idx]
				z := e.Shp.IpRealCoords(e.X, ip)[1]
				_, ρLC, _ := Global.HydroSt.Calc(z)
				chk.Scalar(tst, io.Sf("sl(@ %18g)= %18g", z, s.A_sl), 1e-17, s.A_sl, 1)
				chk.Scalar(tst, io.Sf("ρL(@ %18g)= %18g", z, s.A_ρL), 1e-13, s.A_ρL, ρLC)
			}
		}
	*/
}

func Test_phi02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("phi01b")

	// make sure to flush log
	defer End()

	// run simulation
	if !Start("data/phi02.sim", true, chk.Verbose) {
		tst.Errorf("test failed\n")
		return
	}

	// run simulation
	if !Run() {
		tst.Errorf("test failed\n")
		return
	}
}

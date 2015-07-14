// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_phi01(tst *testing.T) {

	/*     Nodes                           Equations:
	 *                 7                                15
	 *     6 o----o----o----o----o 8       16 o----o----o----o----o 21
	 *       |   13    |   14    |            |   18    |   23    |
	 *       |         |         |            |         |         |
	 *    18 o   [2]   o19 [3]   o 20      19 o    o    o17  o    o 22
	 *       |   23    |   24    |            |   20    |   24    |
	 *       |         |         |            |         |         |
	 *     3 o----o----o----o----o 5        3 o----o----o----o----o 10
	 *       |   11   4|   12    |            |    6   2|   13    |
	 *       |         |         |            |         |         |
	 *    15 o   [0]   o16 [1]   o 17       7 o    o    o5   o    o 12
	 *       |   21    |   22    |            |    8    |   14    |
	 *       |         |         |            |         |         |
	 *     0 o----o----o----o----o 2        0 o----o----o----o----o 9
	 *            9    1    10                     4    1    11
	 */

	//verbose()
	chk.PrintTitle("phi01")

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
	chk.IntAssert(len(dom.Nodes), 25)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 1)
		chk.StrAssert(nod.Dofs[0].Key, "h")
	}
	// check equations
	nids, eqs := get_nids_eqs(dom)

	chk.Ints(tst, "nids", nids, []int{
		0, 1, 4, 3, 9, 16, 11, 15, 21, 2, 5, 10, 17, 12, 22, 7, 6, 19, 13, 18, 23, 8, 20, 14, 24,
	})
	chk.Ints(tst, "eqs", eqs, []int{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
	})

	// check initial values
	r := 0.25
	chk.Scalar(tst, "h @ nod 4", 1e-17, dom.Sol.Y[2], -r)
	chk.Scalar(tst, "h @ nod 8", 1e-17, dom.Sol.Y[21], math.Sqrt(2)/2.0-r)
	for _, eq := range []int{5, 6, 13, 17} {
		chk.Scalar(tst, io.Sf("h @ eq %d", eq), 1e-17, dom.Sol.Y[eq], 0)
	}
}

func Test_phi02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("phi02")

	// run simulation
	defer End()
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

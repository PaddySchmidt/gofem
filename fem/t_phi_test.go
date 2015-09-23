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

	// start simulation
	analysis := NewFEM("data/phi01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.PrintTitle("phi02 - Transport constant Speed")

	// run simulation
	analysis := NewFEM("data/phi02.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

	eq := []int{
		30, 50, 70, 80, 90,
	}
	node := []int{
		15, 25, 35, 40, 45,
	}
	// check results
	for i, v := range eq {
		x := []float64{dom.Msh.Verts[node[i]].C[0], dom.Msh.Verts[node[i]].C[1]}
		xc := []float64{0.0, 0.025}
		r := 1.0
		d := (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1])
		d = math.Sqrt(d) - r

		// chk.PrintTitle(io.Sf("%25.10e%25.10e", dom.Sol.Y[v], d+1.0))
		chk.Scalar(tst, "h @ nod "+string(i), 2e-2, dom.Sol.Y[v], d+1.0)
	}
}

func Test_phi03(tst *testing.T) {

	// verbose()
	chk.PrintTitle("phi03 - Reinitialisation")

	// run simulation
	analysis := NewFEM("data/phi03.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}
	// domain
	dom := analysis.Domains[0]

	eq := []int{
		0, 20, 54,
	}
	node := []int{
		0, 10, 27,
	}
	// check results
	for i, v := range eq {
		// chk.PrintTitle(io.Sf("%25.10e%25.10e", dom.Sol.Y[v], dom.Msh.Verts[node[i]].C[0]))
		chk.Scalar(tst, "h @ nod "+string(i), 5e-2, dom.Sol.Y[v], dom.Msh.Verts[node[i]].C[0])
	}
}

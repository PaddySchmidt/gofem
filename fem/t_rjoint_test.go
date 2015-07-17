// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_rjoint01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("rjoint01. curved line in 3D")

	// initialisation
	if !Start("data/rjoint01.sim", true, chk.Verbose, false) {
		tst.Errorf("Start failed\n")
		return
	}

	// callback to check consistent tangent operators
	eid := 2 // rjoint element
	if true {
		defer rjoint_DebugKb(&testKb{
			tst: tst, eid: eid, tol: 1e-8, verb: chk.Verbose,
			ni: -1, nj: -1, itmin: 1, itmax: -1, tmin: -1, tmax: -1,
		})()
	}

	// run simulation
	if !RunAll() {
		tst.Errorf("Run failed\n")
		return
	}
}

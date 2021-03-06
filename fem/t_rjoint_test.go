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
	analysis := NewFEM("data/rjoint01.sim", "", true, false, false, false, chk.Verbose, 0)

	// callback to check consistent tangent operators
	eid := 2 // rjoint element
	if true {
		rjoint_DebugKb(analysis, &testKb{
			tst: tst, eid: eid, tol: 1e-8, verb: chk.Verbose,
			ni: -1, nj: -1, itmin: 1, itmax: -1, tmin: -1, tmax: -1,
		})
	}

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}
}

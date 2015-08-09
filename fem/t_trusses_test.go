// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_bridge01a(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bridge01a. simple bridge section")

	// fem
	analysis := NewFEM("data/bridge01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-11
	tolu := 1e-15
	tols := 1e-9
	TestingCompareResultsU(tst, "data/bridge01.sim", "cmp/bridge01.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_bridge01b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bridge01. simple bridge section. ElastRod")

	// fem
	analysis := NewFEM("data/bridge01erod.sim", "", true, true, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// recompute matrices
	for _, elem := range analysis.Domains[0].Elems {
		e := elem.(*ElastRod)
		e.Recompute(true)
	}

	// run
	err = analysis.SolveOneStage(0, true)
	if err != nil {
		tst.Error("SolveOneStage failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-11
	tolu := 1e-15
	tols := 1e-9
	TestingCompareResultsU(tst, "data/bridge01erod.sim", "cmp/bridge01.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

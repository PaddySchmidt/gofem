// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_fourlayers01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("fourlayers01")

	if !Start("data/fourlayers.sim", true, chk.Verbose, false) {
		tst.Errorf("test failed\n")
	}

	distr := false
	dom := NewDomain(Global.Sim.Regions[0], distr)
	if dom == nil {
		tst.Errorf("test failed\n")
		return
	}

	io.Pforan("stage # 0\n")
	if !dom.SetStage(0, Global.Sim.Stages[0], distr) {
		tst.Errorf("test failed\n")
		return
	}
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{1, 2, 14, 12, 0, 10})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 12, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11})

	io.Pforan("stage # 1\n")
	if !dom.SetStage(1, Global.Sim.Stages[1], distr) {
		tst.Errorf("test failed\n")
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 19, 4, 5, 20, 6, 7, 8, 9, 18, 10, 11, 12, 13, 14, 15, 16, 17})

	io.Pforan("stage # 2\n")
	if !dom.SetStage(2, Global.Sim.Stages[2], distr) {
		tst.Errorf("test failed\n")
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 7, 11, 8, 13})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 25, 4, 5, 26, 6, 7, 8, 9, 24, 10, 11, 12, 13, 14, 15, 16, 17, 27, 18, 19, 20, 21, 22, 23})

	io.Pforan("stage # 3\n")
	if !dom.SetStage(3, Global.Sim.Stages[3], distr) {
		tst.Errorf("test failed\n")
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{7, 13, 5, 4, 10, 12, 9, 6, 1, 2, 14, 11, 3, 0, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 33, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 31, 12, 13, 32, 14, 15, 16, 17, 30, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29})
}

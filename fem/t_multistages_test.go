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

	analysis := NewFEM("data/fourlayers.sim", "", true, false, false, false, chk.Verbose, 0)

	doms := NewDomains(analysis.Sim, analysis.DynCfs, analysis.HydSta, 0, 1, false)
	if len(doms) == 0 {
		tst.Errorf("NewDomains failed\n")
		return
	}
	dom := doms[0]

	io.Pforan("stage # 0\n")
	err := dom.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage # 0 failed\n%v", err)
		return
	}
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{1, 2, 14, 12, 0, 10})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 12, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11})

	io.Pforan("stage # 1\n")
	err = dom.SetStage(1)
	if err != nil {
		tst.Errorf("SetStage # 1 failed\n%v", err)
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 19, 4, 5, 20, 6, 7, 8, 9, 18, 10, 11, 12, 13, 14, 15, 16, 17})

	io.Pforan("stage # 2\n")
	err = dom.SetStage(2)
	if err != nil {
		tst.Errorf("SetStage # 2 failed\n%v", err)
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 7, 11, 8, 13})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 25, 4, 5, 26, 6, 7, 8, 9, 24, 10, 11, 12, 13, 14, 15, 16, 17, 27, 18, 19, 20, 21, 22, 23})

	io.Pforan("stage # 3\n")
	err = dom.SetStage(3)
	if err != nil {
		tst.Errorf("SetStage # 3 failed\n%v", err)
		return
	}
	nids, eqs = get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{7, 13, 5, 4, 10, 12, 9, 6, 1, 2, 14, 11, 3, 0, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 33, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 31, 12, 13, 32, 14, 15, 16, 17, 30, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29})
}

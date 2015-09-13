// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_nurbs01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs01")

	msh, err := ReadMsh("data", "nurbs01.msh", 0)
	if err != nil {
		tst.Errorf("test failed:\n%v", err)
		return
	}

	for _, cell := range msh.Cells {
		p0, p1 := cell.Shp.Nurbs.Ord(0), cell.Shp.Nurbs.Ord(1)
		io.Pfcyan("cell # %d : NURBS orders = (%d,%d)\n", cell.Id, p0, p1)
		chk.IntAssert(p0, 2)
		chk.IntAssert(p1, 2)
	}
}

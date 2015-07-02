// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_topol01(tst *testing.T) {

	// test title
	verbose()
	chk.PrintTitle("topol01")

	// constants
	datadir := "data/"
	simfn := "box.sim"

	// run FE simulation
	defer fem.End()
	if !fem.Start(datadir+simfn, true, chk.Verbose) {
		chk.Panic("cannot start FE simulation")
	}
	if !fem.Run() {
		chk.Panic("cannot run FE simulation")
	}

	// start post-processing
	Start(datadir+simfn, 0, 0)

	// vertices on plane
	ftag := -31
	dat := PlaneIterator(ftag)
	vids := []int{48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}
	chk.IntAssert(dat.Plane, 2) // perpendicular to z
	chk.IntAssert(len(dat.Vids), 16)
	for _, vid := range vids {
		if !dat.Vids[vid] {
			chk.Panic("vid=%d is not in map", vid)
		}
	}
	chk.Scalar(tst, "Δx", 1e-15, dat.Dx[0], 4.0/3.0)
	chk.Scalar(tst, "Δy", 1e-15, dat.Dx[1], 4.0/3.0)
	chk.Scalar(tst, "Δz", 1e-15, dat.Dx[2], 0)
	chk.Ints(tst, "Iu", dat.Iu, []int{0, 1})
	chk.Scalar(tst, "Δu", 1e-15, dat.Du[0], 4.0/3.0)
	chk.Scalar(tst, "Δv", 1e-15, dat.Du[1], 4.0/3.0)
	chk.Vector(tst, "Umin", 1e-15, dat.Umin, []float64{0, 0})
	chk.Vector(tst, "Umax", 1e-15, dat.Umax, []float64{4, 4})
	chk.Ints(tst, "Nu", dat.Nu, []int{4, 4})

	dat.Iter()

	return

	// define entities
	lx, ly, lz := 4.0, 4.0, 4.0
	Define("A", At{lx, ly, lz})

	// load results
	LoadResults(nil)
	Plot("t", "uz", "A", plt.Fmt{C: "b", M: "o"}, -1)

	//Draw("", "", true, nil)
	Draw("", "", false, nil)
}

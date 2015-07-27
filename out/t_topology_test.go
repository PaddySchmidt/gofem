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

func Test_topo01(tst *testing.T) {

	// test title
	//verbose()
	chk.PrintTitle("topo01")

	// constants
	datadir := "data/"
	simfn := "box.sim"

	// start simulation
	processing := fem.NewFEM(datadir+simfn, "", true, true, false, false, chk.Verbose)

	// run simulation
	err := processing.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// start post-processing
	Start(datadir+simfn, 0, 0)

	// vertices on plane (indenter/surface)
	ftag := -31
	surf := NodesOnPlane(ftag)
	vids := []int{48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}
	chk.IntAssert(surf.Plane, 2) // perpendicular to z
	chk.IntAssert(len(surf.Ids), 16)
	for _, vid := range vids {
		if !surf.Ids[vid] {
			chk.Panic("vid=%d is not in map", vid)
		}
	}
	chk.Scalar(tst, "Δx", 1e-15, surf.Dx[0], 4.0/3.0)
	chk.Scalar(tst, "Δy", 1e-15, surf.Dx[1], 4.0/3.0)
	chk.Scalar(tst, "Δz", 1e-15, surf.Dx[2], 0)
	chk.Ints(tst, "Iu", surf.Iu, []int{0, 1})
	chk.Scalar(tst, "Δu", 1e-15, surf.Du[0], 4.0/3.0)
	chk.Scalar(tst, "Δv", 1e-15, surf.Du[1], 4.0/3.0)
	chk.Vector(tst, "Umin", 1e-15, surf.Umin, []float64{0, 0})
	chk.Vector(tst, "Umax", 1e-15, surf.Umax, []float64{4, 4})
	chk.Ints(tst, "Nu", surf.Nu, []int{4, 4})

	// define entities
	lx, ly, lz := 4.0, 4.0, 4.0
	Define("A", At{lx, ly, lz})
	Define("surf", surf)

	// load results
	Extrap = []string{"sz"}
	LoadResults(nil)

	// displacements
	Splot("displacements")
	Plot("t", "uz", "A", plt.Fmt{C: "b", M: "o"}, -1)

	// stresses on surface
	V := IntegOnPlane("ex_sz", "surf")
	Splot("integral of stresses")
	Plot(Times, V, "surf", plt.Fmt{C: "r", M: "o"}, -1)
	Csplot.Xlbl = "time"
	Csplot.Ylbl = "integ(sz)"

	//Draw("", "", true, nil)
	Draw("", "", false, nil)
}

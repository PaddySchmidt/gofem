// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
)

func Test_nurbs01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs01")

	verts := [][]float64{
		{5.0, 10, 0, 1}, // global 0
		{8.0, 10, 0, 1}, // global 1
		{8.0, 13, 0, 1}, // global 2
		{5.0, 13, 0, 1}, // global 3
		{6.0, 10, 0, 1}, // global 4
		{6.0, 13, 0, 1}, // global 5
		{7.0, 10, 0, 1}, // global 6
		{7.0, 13, 0, 1}, // global 7
	}
	knots := [][]float64{
		{0, 0, 0, 0.5, 1, 1, 1},
		{0, 0, 1, 1},
	}
	ctrls := []int{
		0, 4, 6, 1, // first level along x
		3, 5, 7, 2, // second level along x
	}

	var nurbs gm.Nurbs
	nurbs.Init(2, []int{2, 1}, knots)
	nurbs.SetControl(verts, ctrls)

	spans := nurbs.Elements()
	ibasis0 := nurbs.IndBasis(spans[0])
	ibasis1 := nurbs.IndBasis(spans[1])
	io.Pforan("spans = %v\n", spans)
	chk.Ints(tst, "span0", spans[0], []int{2, 3, 1, 2})
	chk.Ints(tst, "span1", spans[1], []int{3, 4, 1, 2})
	chk.Ints(tst, "ibasis0", ibasis0, []int{0, 1, 2, 4, 5, 6})
	chk.Ints(tst, "ibasis1", ibasis1, []int{1, 2, 3, 5, 6, 7})

	o := GetShapeNurbs(&nurbs)

	dux := 0.5
	duy := 1.0
	drx := 2.0
	dry := 2.0
	JuCor := (dux / drx) * (duy / dry)

	r := []float64{0.75, 0.75, 0}

	Ju, u, ibasis := NurbsShapeFunc(o.S, o.DSdR, r, true, &nurbs, spans[0])
	io.Pforan("Ju = %v\n", Ju)
	io.Pforan("u = %v\n", u)
	io.Pforan("ibasis = %v\n", ibasis)
	chk.Scalar(tst, "Ju", 1e-17, Ju, JuCor)
	chk.Scalar(tst, "ux", 1e-17, u[0], (1.0+r[0])*dux/drx)
	chk.Scalar(tst, "uy", 1e-17, u[1], (1.0+r[1])*duy/dry)
	chk.Ints(tst, "ibasis", ibasis, []int{0, 1, 2, 4, 5, 6})

	io.Pforan("S(u(r)) = %v\n", o.S)

	Ju, u, ibasis = NurbsShapeFunc(o.S, o.DSdR, r, true, &nurbs, spans[1])
	io.Pfpink("\nJu = %v\n", Ju)
	io.Pfpink("u = %v\n", u)
	io.Pfpink("ibasis = %v\n", ibasis)
	chk.Scalar(tst, "Ju", 1e-17, Ju, JuCor)
	chk.Scalar(tst, "ux", 1e-17, u[0], 0.5+(1.0+r[0])*dux/drx)
	chk.Scalar(tst, "uy", 1e-17, u[1], (1.0+r[1])*duy/dry)
	chk.Ints(tst, "ibasis", ibasis, []int{1, 2, 3, 5, 6, 7})

	if chk.Verbose && false {
		gm.PlotNurbs("/tmp/gofem", "tst_nurbs01", &nurbs)
	}
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

func get_nurbs_A() *gm.Nurbs {
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
	return &nurbs
}

func get_nurbs_B() *gm.Nurbs {
	verts := [][]float64{
		{5.0, 10, 0, 1}, // global 0
		{8.0, 10, 0, 1}, // global 1
		{8.0, 13, 0, 1}, // global 2
		{5.0, 13, 0, 1}, // global 3
		{6.0, 11, 0, 1}, // global 4
		{6.0, 12, 0, 1}, // global 5
		{7.0, 11, 0, 1}, // global 6
		{7.0, 12, 0, 1}, // global 7
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
	return &nurbs
}

func Test_nurbs01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs01")

	nurbs := get_nurbs_A()
	spans := nurbs.Elements()
	ibasis0 := nurbs.IndBasis(spans[0])
	ibasis1 := nurbs.IndBasis(spans[1])
	io.Pforan("spans = %v\n", spans)
	chk.Ints(tst, "span0", spans[0], []int{2, 3, 1, 2})
	chk.Ints(tst, "span1", spans[1], []int{3, 4, 1, 2})
	chk.Ints(tst, "ibasis0", ibasis0, []int{0, 1, 2, 4, 5, 6})
	chk.Ints(tst, "ibasis1", ibasis1, []int{1, 2, 3, 5, 6, 7})

	shape := GetShapeNurbs(nurbs)

	dux := 0.5
	duy := 1.0
	drx := 2.0
	dry := 2.0
	JuCor := (dux / drx) * (duy / dry)

	r := []float64{0.75, 0.75, 0}

	Ju, u, ibasis := shape.NurbsFunc(shape.S, shape.DSdR, r, true, spans[0])
	io.Pforan("Ju = %v\n", Ju)
	io.Pforan("u = %v\n", u)
	io.Pforan("ibasis = %v\n", ibasis)
	chk.Scalar(tst, "Ju", 1e-17, Ju, JuCor)
	chk.Scalar(tst, "ux", 1e-17, u[0], (1.0+r[0])*dux/drx)
	chk.Scalar(tst, "uy", 1e-17, u[1], (1.0+r[1])*duy/dry)
	chk.Ints(tst, "ibasis", ibasis, []int{0, 1, 2, 4, 5, 6})

	io.Pforan("S(u(r)) = %v\n", shape.S)

	Ju, u, ibasis = shape.NurbsFunc(shape.S, shape.DSdR, r, true, spans[1])
	io.Pfpink("\nJu = %v\n", Ju)
	io.Pfpink("u = %v\n", u)
	io.Pfpink("ibasis = %v\n", ibasis)
	chk.Scalar(tst, "Ju", 1e-17, Ju, JuCor)
	chk.Scalar(tst, "ux", 1e-17, u[0], 0.5+(1.0+r[0])*dux/drx)
	chk.Scalar(tst, "uy", 1e-17, u[1], (1.0+r[1])*duy/dry)
	chk.Ints(tst, "ibasis", ibasis, []int{1, 2, 3, 5, 6, 7})

	if chk.Verbose {
		gm.PlotNurbs("/tmp/gofem", "tst_nurbs01", nurbs)
	}
}

func Test_nurbs02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs02")

	shape := GetShapeNurbs(get_nurbs_A())

	C := [][][]float64{
		{{5, 10}, {6.5, 10}, {6.5, 13}, {5, 13}},
		{{6.5, 10}, {8, 10}, {8, 13}, {6.5, 13}},
	}

	tol := 1e-14
	verb := true
	check_nurbs_isoparametric(tst, shape, C)
	check_nurbs_dSdR(tst, shape, []float64{0.75, 0.75, 0}, tol, verb)
}

func Test_nurbs03(tst *testing.T) {

	//verbose()
	chk.PrintTitle("nurbs03")

	shape := GetShapeNurbs(get_nurbs_B())

	C := [][][]float64{
		{{5, 10}, {6.5, 11}, {6.5, 12}, {5, 13}},
		{{6.5, 11}, {8, 10}, {8, 13}, {6.5, 12}},
	}

	tol := 1e-14
	verb := true
	check_nurbs_isoparametric(tst, shape, C)
	check_nurbs_dSdR(tst, shape, []float64{0.75, 0.75, 0}, tol, verb)

	if false {
		gm.PlotNurbs("/tmp/gofem", "tst_nurbs03", shape.Nurbs)
	}
}

// check isoparametric property
//  C -- [nspans_or_elements][4][2] elements coordinates of corners (not control points)
func check_nurbs_isoparametric(tst *testing.T, shape *Shape, C [][][]float64) {

	// auxiliary
	r := []float64{0, 0, 0}
	x := make([]float64, 2)
	qua4_natcoords := [][]float64{
		{-1, 1, 1, -1},
		{-1, -1, 1, 1},
	}

	// loop over elements == spans
	spans := shape.Nurbs.Elements()
	for ie, span := range spans {
		ibasis := shape.Nurbs.IndBasis(span)
		io.Pf("\nelement = %v, ibasis = %v\n", span, ibasis)
		for i := 0; i < 4; i++ {
			for j := 0; j < 2; j++ {
				r[j] = qua4_natcoords[j][i]
			}
			shape.NurbsFunc(shape.S, shape.DSdR, r, false, span)
			for j := 0; j < 2; j++ {
				x[j] = 0
				for k, l := range ibasis {
					q := shape.Nurbs.GetQl(l)
					x[j] += shape.S[k] * q[j]
				}
			}
			io.Pforan("x = %v\n", x)
			chk.Vector(tst, "x", 1e-17, x, C[ie][i])
		}
	}
}

func check_nurbs_dSdR(tst *testing.T, shape *Shape, r []float64, tol float64, verbose bool) {

	// auxiliary
	r_tmp := make([]float64, len(r))
	S_tmp := make([]float64, shape.Nverts)

	// loop over elements == spans
	spans := shape.Nurbs.Elements()
	for _, span := range spans {
		ibasis := shape.Nurbs.IndBasis(span)
		io.Pf("\nelement = %v, ibasis = %v\n", span, ibasis)

		// analytical
		shape.NurbsFunc(shape.S, shape.DSdR, r, true, span)

		// numerical
		for n := 0; n < shape.Nverts; n++ {
			for i := 0; i < shape.Gndim; i++ {
				dSndRi, _ := num.DerivCentral(func(t float64, args ...interface{}) (Sn float64) {
					copy(r_tmp, r)
					r_tmp[i] = t
					shape.NurbsFunc(S_tmp, nil, r_tmp, false, span)
					Sn = S_tmp[n]
					return
				}, r[i], 1e-1)
				if verbose {
					io.Pfgrey2("  dS%ddR%d @ [%5.2f%5.2f%5.2f] = %v (num: %v)\n", n, i, r[0], r[1], r[2], shape.DSdR[n][i], dSndRi)
				}
				if math.Abs(shape.DSdR[n][i]-dSndRi) > tol {
					tst.Errorf("nurbs dS%ddR%d failed with err = %g\n", n, i, math.Abs(shape.DSdR[n][i]-dSndRi))
					return
				}
				//chk.Scalar(tst, fmt.Sprintf("dS%ddR%d", n, i), tol, dSdR[n][i], dSndRi)
			}
		}
	}
}

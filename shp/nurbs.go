// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
)

// GetShapeNurbs returns a shape structure based on NURBS
//  Note: enodes are the local ids of control points in NURBS
func GetShapeNurbs(nurbs *gm.Nurbs) (o *Shape) {
	o = new(Shape)
	o.Type = "nurbs"
	o.Gndim = nurbs.Gnd()
	switch o.Gndim {
	case 1:
		o.BasicType = "lin2"
		o.NatCoords = [][]float64{
			{-1, 1},
		}
	case 2:
		o.BasicType = "qua4"
		o.NatCoords = [][]float64{
			{-1, 1, 1, -1},
			{-1, -1, 1, 1},
		}
	case 3:
		o.BasicType = "hex8"
		o.NatCoords = [][]float64{
			{-1, 1, 1, -1, -1, 1, 1, -1},
			{-1, -1, 1, 1, -1, -1, 1, 1},
			{-1, -1, -1, -1, 1, 1, 1, 1},
		}
	}
	o.Nverts = nurbs.GetElemNumBasis()
	o.VtkCode = VTK_POLY_VERTEX
	if o.Gndim > 1 {
		o.FaceType = "nurbs"
	}
	o.init_scratchpad()
	return
}

func NurbsShapeFunc(S []float64, dSdR [][]float64, r []float64, derivs bool, nurbs *gm.Nurbs, span []int) (Ju float64, u []float64, ibasis []int) {

	// compute mapping to knots space
	nd := nurbs.Gnd()
	u = make([]float64, nd)
	Ju = 1.0 // det(dudr) => du = Ju * dr
	var umin, umax float64
	for i := 0; i < nd; i++ {
		umin = nurbs.U(i, span[i*2])
		umax = nurbs.U(i, span[i*2+1])
		u[i] = ((umax-umin)*r[i] + (umax + umin)) / 2.0
		Ju *= (umax - umin) / 2.0
		if u[i] < umin || u[i] > umax {
			chk.Panic("compute NURBS shape function outide cell range:\nr[%d]=%v, u[%d]=%v, urange=[%v,%v]", i, r[i], i, u[i], umin, umax)
		}
	}

	// local indices of control points
	ibasis = nurbs.IndBasis(span)

	// shape and/or derivatives in knots space
	if derivs {
		nurbs.CalcBasisAndDerivs(u)
	} else {
		nurbs.CalcBasis(u)
	}
	for k, l := range ibasis {
		S[k] = nurbs.GetBasisL(l)
	}

	// derivatives in natural space
	if derivs {
		for k, l := range ibasis {
			nurbs.GetDerivL(dSdR[k], l) // dSdR := dSdU
			for i := 0; i < nd; i++ {
				umin = nurbs.U(i, span[i*2])
				umax = nurbs.U(i, span[i*2+1])
				dSdR[k][i] *= (umax - umin) / 2.0 // dSdR[i] := dSdU[i] * du[i]/dr[i] (no sum on i)
			}
		}
	}
	return
}

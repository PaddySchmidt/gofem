// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
)

// GetShapeNurbs returns a shape structure based on NURBS
//  Note: span are the local ids of control points in NURBS defining elements
func GetShapeNurbs(nurbs *gm.Nurbs, span []int) (o *Shape) {
	o = new(Shape)
	o.Type = "nurbs"
	o.Gndim = nurbs.Gnd()
	switch o.Gndim {
	case 1:
		o.BasicType = "lin2"
	case 2:
		o.BasicType = "qua4"
	case 3:
		o.BasicType = "hex8"
	}
	o.Nverts = nurbs.GetElemNumBasis()
	o.VtkCode = VTK_POLY_VERTEX
	if o.Gndim > 1 {
		o.FaceType = "nurbs"
	}
	o.Func = o.NurbsFunc
	o.Nurbs = nurbs
	o.Span = span
	o.Ibasis = o.Nurbs.IndBasis(o.Span)
	o.U = make([]float64, o.Gndim)
	o.init_scratchpad()
	return
}

// NurbsFunc implements shape/deriv functions for NURBS
func (o *Shape) NurbsFunc(S []float64, dSdR [][]float64, r []float64, derivs bool) {

	// compute mapping to knots space
	o.Ju = 1.0 // det(dudr) => du = Ju * dr
	var umin, umax float64
	for i := 0; i < o.Gndim; i++ {
		umin = o.Nurbs.U(i, o.Span[i*2])
		umax = o.Nurbs.U(i, o.Span[i*2+1])
		o.U[i] = ((umax-umin)*r[i] + (umax + umin)) / 2.0
		o.Ju *= (umax - umin) / 2.0
		if o.U[i] < umin || o.U[i] > umax {
			chk.Panic("compute NURBS shape function outide cell range:\nr[%d]=%v, u[%d]=%v, urange=[%v,%v]", i, r[i], i, o.U[i], umin, umax)
		}
	}

	// shape and/or derivatives in knots space
	if derivs {
		o.Nurbs.CalcBasisAndDerivs(o.U)
	} else {
		o.Nurbs.CalcBasis(o.U)
	}
	for k, l := range o.Ibasis {
		S[k] = o.Nurbs.GetBasisL(l)
	}

	// derivatives in natural space
	if derivs {
		for k, l := range o.Ibasis {
			o.Nurbs.GetDerivL(dSdR[k], l) // dSdR := dSdU
			for i := 0; i < o.Gndim; i++ {
				umin = o.Nurbs.U(i, o.Span[i*2])
				umax = o.Nurbs.U(i, o.Span[i*2+1])
				dSdR[k][i] *= (umax - umin) / 2.0 // dSdR[i] := dSdU[i] * du[i]/dr[i] (no sum on i)
			}
		}
	}
	return
}

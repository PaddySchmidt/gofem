// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/utl"
)

// GetShapeNurbs returns a shape structure based on NURBS
//  Note: span are the local ids of control points in NURBS defining elements
//  Note: FaceLocalVerts does not work for internal surfaces; only those @ boundaries
func GetShapeNurbs(nurbs *gm.Nurbs, nrbfaces []*gm.Nurbs, span []int) (o *Shape) {

	// basic data
	o = new(Shape)
	o.Type = "nurbs"
	o.FaceType = "nurbs"
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
	o.Func = o.NurbsFunc
	o.FaceFunc = o.NurbsFaceFunc
	o.Nurbs = nurbs
	o.Span = span
	o.Ibasis = o.Nurbs.IndBasis(o.Span)
	o.U = make([]float64, o.Gndim)

	// faces basic data
	nfaces := 2 * o.Gndim
	o.FaceLocalVerts = nurbs.ElemBryLocalInds()
	if o.Gndim == 3 {
		o.NurbsFaces = nrbfaces
		o.SpanFace = [][]int{
			span[0:2], span[0:2],
			span[2:4], span[2:4],
			span[4:6], span[4:6],
		}
		o.FaceFlip = []bool{false, true, false, true, false, true} // => point to the inside
	} else {
		o.NurbsFaces = []*gm.Nurbs{nrbfaces[2], nrbfaces[1], nrbfaces[3], nrbfaces[0]}
		o.SpanFace = [][]int{span[0:2], span[2:4], span[0:2], span[2:4]}
		o.FaceFlip = []bool{false, false, true, true} // => point to the inside
	}
	o.IbasisFace = make([][]int, nfaces)
	for idxface, face := range o.NurbsFaces {
		o.IbasisFace[idxface] = face.IndBasis(o.SpanFace[idxface])
		if idxface == 0 {
			o.FaceNvertsMax = len(o.IbasisFace[idxface])
		} else {
			o.FaceNvertsMax = utl.Imax(o.FaceNvertsMax, len(o.IbasisFace[idxface]))
		}
	}

	// allocate stracthpad variables
	o.init_scratchpad()
	return
}

// nurbs_func implements shape/deriv functions for NURBS
func nurbs_func(u, S []float64, dSdR [][]float64, r []float64, derivs bool, nurbs *gm.Nurbs, ibasis, span []int) {

	// compute mapping to knots space
	nd := nurbs.Gnd()
	var umin, umax float64
	for i := 0; i < nd; i++ {
		umin = nurbs.U(i, span[i*2])
		umax = nurbs.U(i, span[i*2+1])
		u[i] = ((umax-umin)*r[i] + (umax + umin)) / 2.0
		if u[i] < umin || u[i] > umax {
			chk.Panic("cannot compute NURBS shape function outide cell range:\nr[%d]=%v, u[%d]=%v, urange=[%v,%v]", i, r[i], i, u[i], umin, umax)
		}
	}

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

// NurbsFunc implements shape/deriv functions for NURBS
func (o *Shape) NurbsFunc(S []float64, dSdR [][]float64, r []float64, derivs bool, idxface int) {
	nurbs_func(o.U, S, dSdR, r, derivs, o.Nurbs, o.Ibasis, o.Span)
}

// NurbsFaceFunc implements shape/deriv functions for faces of NURBS
func (o *Shape) NurbsFaceFunc(S []float64, dSdR [][]float64, r []float64, derivs bool, idxface int) {
	nurbs_func(o.U, S, dSdR, r, derivs, o.NurbsFaces[idxface], o.IbasisFace[idxface], o.SpanFace[idxface])
}

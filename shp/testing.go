// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// Check_S checks that shape functions evaluate to 1.0 @ nodes
func (o *Shape) Check_S(tst *testing.T, tol float64, verbose bool) {

	// loop over all vertices
	errS := 0.0
	r := []float64{0, 0, 0}
	for n := 0; n < o.Nverts; n++ {

		// natural coordinates @ vertex
		for i := 0; i < o.Gndim; i++ {
			r[i] = o.NatCoords[i][n]
		}

		// compute function
		o.Func(o.S, o.DSdR, r, false)

		// check
		if verbose {
			io.Pf("S = %v\n", o.S)
		}
		for m := 0; m < o.Nverts; m++ {
			if n == m {
				errS += math.Abs(o.S[m] - 1.0)
			} else {
				errS += math.Abs(o.S[m])
			}
		}
	}

	// error
	if errS > tol {
		tst.Errorf("%s failed with err = %g\n", o.Type, errS)
		return
	}
}

// Check_Sf checks shape functions @ faces
func (o *Shape) Check_Sf(tst *testing.T, tol float64, verbose bool) {

	// skip 1D shapes
	nfaces := len(o.FaceLocalV)
	if nfaces == 0 {
		return
	}

	// loop over face vertices
	errS := 0.0
	r := []float64{0, 0, 0}
	for k := 0; k < nfaces; k++ {
		for n := range o.FaceLocalV[k] {

			// natural coordinates @ vertex
			for i := 0; i < o.Gndim; i++ {
				r[i] = o.NatCoords[i][n]
			}

			// compute function
			o.Func(o.S, o.DSdR, r, false)

			// check
			if verbose {
				io.Pforan("S = %v\n", o.S)
			}
			for m := range o.FaceLocalV[k] {
				if n == m {
					errS += math.Abs(o.S[m] - 1.0)
				} else {
					errS += math.Abs(o.S[m])
				}
			}
		}
	}

	// error
	if verbose {
		io.Pforan("%g\n", errS)
	}
	if errS > tol {
		tst.Errorf("%s failed with err = %g\n", o.Type, errS)
		return
	}
}

// Check_dSdR compares analytical dSdR with numerical results obtained by finite differences
func (o *Shape) Check_dSdR(tst *testing.T, tol float64, verbose bool) {

	// loop over all vertices
	h := 1e-1
	r := []float64{0, 0, 0}
	r_temp := []float64{0, 0, 0}
	S_temp := make([]float64, o.Nverts)
	for n := 0; n < o.Nverts; n++ {

		// natural coordinates @ vertex
		for i := 0; i < o.Gndim; i++ {
			r[i] = o.NatCoords[i][n]
		}

		// analytical
		o.Func(o.S, o.DSdR, r, true)

		// numerical
		for i := 0; i < o.Gndim; i++ {
			dSndRi, _ := num.DerivCentral(func(x float64, args ...interface{}) (Sn float64) {
				copy(r_temp, r)
				r_temp[i] = x
				o.Func(S_temp, nil, r_temp, false)
				Sn = S_temp[n]
				return
			}, r[i], h)
			if verbose {
				io.Pfgrey2("  dS%ddR%d @ [% 4.1f % 4.1f % 4.1f] = %v (num: %v)\n", n, i, r[0], r[1], r[2], o.DSdR[n][i], dSndRi)
			}
			if math.Abs(o.DSdR[n][i]-dSndRi) > tol {
				tst.Errorf("%s dS%ddR%d failed with err = %g\n", o.Type, n, i, math.Abs(o.DSdR[n][i]-dSndRi))
				return
			}
			//chk.Scalar(tst, fmt.Sprintf("dS%ddR%d", n, i), tol, dSdR[n][i], dSndRi)
		}
	}
}

// CheckNurbsIsop checks isoparametric property with NURBS
//  C -- [4][2] elements coordinates of corners (not control points)
func CheckNurbsIsop(tst *testing.T, shape *Shape, C [][]float64) {

	// auxiliary
	r := []float64{0, 0, 0}
	x := make([]float64, 2)
	qua4_natcoords := [][]float64{
		{-1, 1, 1, -1},
		{-1, -1, 1, 1},
	}

	// check
	io.Pf("\nelement = %v, ibasis = %v\n", shape.Span, shape.Ibasis)
	for i := 0; i < 4; i++ {
		for j := 0; j < 2; j++ {
			r[j] = qua4_natcoords[j][i]
		}
		shape.NurbsFunc(shape.S, shape.DSdR, r, false)
		for j := 0; j < 2; j++ {
			x[j] = 0
			for k, l := range shape.Ibasis {
				q := shape.Nurbs.GetQl(l)
				x[j] += shape.S[k] * q[j]
			}
		}
		io.Pforan("x = %v\n", x)
		chk.Vector(tst, "x", 1e-17, x, C[i])
	}
}

// CheckNurbs_dSdR checks derivatives with NURBS
func CheckNurbs_dSdR(tst *testing.T, shape *Shape, r []float64, tol float64, verbose bool) {

	// auxiliary
	r_tmp := make([]float64, len(r))
	S_tmp := make([]float64, shape.Nverts)

	// loop over elements == spans
	spans := shape.Nurbs.Elements()
	for _, span := range spans {
		ibasis := shape.Nurbs.IndBasis(span)
		io.Pf("\nelement = %v, ibasis = %v\n", span, ibasis)

		// analytical
		shape.NurbsFunc(shape.S, shape.DSdR, r, true)

		// numerical
		for n := 0; n < shape.Nverts; n++ {
			for i := 0; i < shape.Gndim; i++ {
				dSndRi, _ := num.DerivCentral(func(t float64, args ...interface{}) (Sn float64) {
					copy(r_tmp, r)
					r_tmp[i] = t
					shape.NurbsFunc(S_tmp, nil, r_tmp, false)
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

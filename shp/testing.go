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

// CheckShape checks that shape functions evaluate to 1.0 @ nodes
func CheckShape(tst *testing.T, shape *Shape, tol float64, verbose bool) {

	// loop over all vertices
	errS := 0.0
	r := []float64{0, 0, 0}
	for n := 0; n < shape.Nverts; n++ {

		// natural coordinates @ vertex
		for i := 0; i < shape.Gndim; i++ {
			r[i] = shape.NatCoords[i][n]
		}

		// compute function
		shape.Func(shape.S, shape.DSdR, r, false)

		// check
		if verbose {
			io.Pf("S = %v\n", shape.S)
		}
		for m := 0; m < shape.Nverts; m++ {
			if n == m {
				errS += math.Abs(shape.S[m] - 1.0)
			} else {
				errS += math.Abs(shape.S[m])
			}
		}
	}

	// error
	if errS > tol {
		tst.Errorf("%s failed with err = %g\n", shape.Type, errS)
		return
	}
}

// CheckShapeFace checks shape functions @ faces
func CheckShapeFace(tst *testing.T, shape *Shape, tol float64, verbose bool) {

	// skip 1D shapes
	nfaces := len(shape.FaceLocalV)
	if nfaces == 0 {
		return
	}

	// loop over face vertices
	errS := 0.0
	r := []float64{0, 0, 0}
	for k := 0; k < nfaces; k++ {
		for n := range shape.FaceLocalV[k] {

			// natural coordinates @ vertex
			for i := 0; i < shape.Gndim; i++ {
				r[i] = shape.NatCoords[i][n]
			}

			// compute function
			shape.Func(shape.S, shape.DSdR, r, false)

			// check
			if verbose {
				io.Pforan("S = %v\n", shape.S)
			}
			for m := range shape.FaceLocalV[k] {
				if n == m {
					errS += math.Abs(shape.S[m] - 1.0)
				} else {
					errS += math.Abs(shape.S[m])
				}
			}
		}
	}

	// error
	if verbose {
		io.Pforan("%g\n", errS)
	}
	if errS > tol {
		tst.Errorf("%s failed with err = %g\n", shape.Type, errS)
		return
	}
}

// CheckIsop checks isoparametric property
//  C    -- [4][2] elements coordinates of corners (not control points)
//  Cnat -- [2][4] natural coordinates of corners
func CheckIsop(tst *testing.T, shape *Shape, C [][]float64, Cnat [][]float64) {

	// auxiliary
	r := []float64{0, 0, 0}
	x := make([]float64, 2)

	// check
	io.Pf("\nelement = %v, ibasis = %v\n", shape.Span, shape.Ibasis)
	for i := 0; i < 4; i++ {
		for j := 0; j < 2; j++ {
			r[j] = Cnat[j][i]
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

// CheckDSdR checks dSdR derivatives of shape structures
func CheckDSdR(tst *testing.T, shape *Shape, r []float64, tol float64, verbose bool) {

	// auxiliary
	r_tmp := make([]float64, len(r))
	S_tmp := make([]float64, shape.Nverts)

	// analytical
	shape.Func(shape.S, shape.DSdR, r, true)

	// numerical
	for n := 0; n < shape.Nverts; n++ {
		for i := 0; i < shape.Gndim; i++ {
			dSndRi, _ := num.DerivCentral(func(t float64, args ...interface{}) (Sn float64) {
				copy(r_tmp, r)
				r_tmp[i] = t
				shape.Func(S_tmp, nil, r_tmp, false)
				Sn = S_tmp[n]
				return
			}, r[i], 1e-1)
			if verbose {
				io.Pfgrey2("  dS%ddR%d @ %5.2f = %v (num: %v)\n", n, i, r, shape.DSdR[n][i], dSndRi)
			}
			if math.Abs(shape.DSdR[n][i]-dSndRi) > tol {
				tst.Errorf("nurbs dS%ddR%d failed with err = %g\n", n, i, math.Abs(shape.DSdR[n][i]-dSndRi))
				return
			}
		}
	}
}

// CheckDSdx checks G=dSdx derivatives of shape structures
func CheckDSdx(tst *testing.T, shape *Shape, xmat [][]float64, x []float64, tol float64, verbose bool) {

	// find r corresponding to x
	r := make([]float64, 3)
	err := shape.InvMap(r, x, xmat)
	if err != nil {
		tst.Errorf("InvMap failed:\n%v", err)
		return
	}

	// analytical
	err = shape.CalcAtIp(xmat, r, true)
	if err != nil {
		tst.Errorf("CalcAtIp failed:\n%v", err)
		return
	}

	// numerical
	x_tmp := make([]float64, len(x))
	for n := 0; n < shape.Nverts; n++ {
		for i := 0; i < shape.Gndim; i++ {
			dSnDxi, _ := num.DerivCentral(func(t float64, args ...interface{}) (Sn float64) {
				copy(x_tmp, x)
				x_tmp[i] = t
				err = shape.InvMap(r, x_tmp, xmat)
				if err != nil {
					tst.Errorf("InvMap failed:\n%v", err)
					return
				}
				err = shape.CalcAtIp(xmat, r, false)
				if err != nil {
					tst.Errorf("CalcAtIp failed:\n%v", err)
					return
				}
				Sn = shape.S[n]
				return
			}, x[i], 1e-1)
			if verbose {
				io.Pfgrey2("  dS%dDx%d @ %5.2f = %v (num: %v)\n", n, i, x, shape.G[n][i], dSnDxi)
			}
			if math.Abs(shape.G[n][i]-dSnDxi) > tol {
				tst.Errorf("nurbs dS%dDx%d failed with err = %g\n", n, i, math.Abs(shape.G[n][i]-dSnDxi))
				return
			}
		}
	}
}

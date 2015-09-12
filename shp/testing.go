// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math"
	"testing"

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

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/tsr"
)

// Mmatch computes M=q/p and qy0 from c and φ corresponding to the strength that would
// be modelled by the Mohr-Coulomb model matching one of the following cones:
//  typ == 0 : compression cone (outer)
//      == 1 : extension cone (inner)
//      == 2 : plane-strain
func Mmatch(c, φ float64, typ int) (M, qy0 float64, err error) {
	φr := φ * math.Pi / 180.0
	si := math.Sin(φr)
	co := math.Cos(φr)
	var ξ float64
	switch typ {
	case 0: // compression cone (outer)
		M = 6.0 * si / (3.0 - si)
		ξ = 6.0 * co / (3.0 - si)
	case 1: // extension cone (inner)
		M = 6.0 * si / (3.0 + si)
		ξ = 6.0 * co / (3.0 + si)
	case 2: // plane-strain
		t := si / co
		d := math.Sqrt(3.0 + 4.0*t*t)
		M = 3.0 * t / d
		ξ = 3.0 / d
	default:
		return 0, 0, chk.Err("typ=%d is invalid", typ)
	}
	qy0 = ξ * c
	return
}

// SpectralCompose recreates tensor m from its spectral decomposition
// m   -- 2nd order tensor in Mandel basis
// λ   -- eigenvalues
// n   -- eigenvectors [ncp][nvecs]
// tmp -- temporary matrix [3][3]
func SpectralCompose(m, λ []float64, n, tmp [][]float64) {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			tmp[i][j] = λ[0]*n[i][0]*n[j][0] + λ[1]*n[i][1]*n[j][1] + λ[2]*n[i][2]*n[j][2]
		}
	}
	tsr.Ten2Man(m, tmp)
}

// Eigenprojectors computes the Mandel eigenprojectors for given eigenvectors
// n  -- eigenvectors [ncp][nvecs]
func Eigenprojectors(P [][]float64, n [][]float64) {
	P[0][0] = n[0][0] * n[0][0]
	P[0][1] = n[1][0] * n[1][0]
	P[0][2] = n[2][0] * n[2][0]
	P[0][3] = n[0][0] * n[1][0] * tsr.SQ2

	P[1][0] = n[0][1] * n[0][1]
	P[1][1] = n[1][1] * n[1][1]
	P[1][2] = n[2][1] * n[2][1]
	P[1][3] = n[0][1] * n[1][1] * tsr.SQ2

	P[2][0] = n[0][2] * n[0][2]
	P[2][1] = n[1][2] * n[1][2]
	P[2][2] = n[2][2] * n[2][2]
	P[2][3] = n[0][2] * n[1][2] * tsr.SQ2
	if len(P[0]) == 6 {
		P[0][4] = n[1][0] * n[2][0] * tsr.SQ2
		P[0][5] = n[2][0] * n[0][0] * tsr.SQ2

		P[1][4] = n[1][1] * n[2][1] * tsr.SQ2
		P[1][5] = n[2][1] * n[0][1] * tsr.SQ2

		P[2][4] = n[1][2] * n[2][2] * tsr.SQ2
		P[2][5] = n[2][2] * n[0][2] * tsr.SQ2
	}
}

/*
func Eigenprojectors(P0, P1, P2 []float64, n [][]float64) {
	P0[0] = n[0][0] * n[0][0]
	P0[1] = n[1][0] * n[1][0]
	P0[2] = n[2][0] * n[2][0]
	P0[3] = n[0][0] * n[1][0] * tsr.SQ2

	P1[0] = n[0][1] * n[0][1]
	P1[1] = n[1][1] * n[1][1]
	P1[2] = n[2][1] * n[2][1]
	P1[3] = n[0][1] * n[1][1] * tsr.SQ2

	P2[0] = n[0][2] * n[0][2]
	P2[1] = n[1][2] * n[1][2]
	P2[2] = n[2][2] * n[2][2]
	P2[3] = n[0][2] * n[1][2] * tsr.SQ2
	if len(P0) == 6 {
		P0[4] = n[1][0] * n[2][0] * tsr.SQ2
		P0[5] = n[2][0] * n[0][0] * tsr.SQ2

		P1[4] = n[1][1] * n[2][1] * tsr.SQ2
		P1[5] = n[2][1] * n[0][1] * tsr.SQ2

		P2[4] = n[1][2] * n[2][2] * tsr.SQ2
		P2[5] = n[2][2] * n[0][2] * tsr.SQ2
	}
}
*/

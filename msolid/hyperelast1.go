// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// HyperElast1 implements a nonlinear hyperelastic model for powders and porous media
type HyperElast1 struct {

	// parameters
	κ   float64 // κ
	κb  float64 // \bar{κ}
	G0  float64 // G0
	pr  float64 // pr
	pt  float64 // pt
	p0  float64 // p0
	εv0 float64 // εv0

	// derived
	pa float64 // pa = pr + pt
	a  float64 // a = 1 / κ
}

// add model to factory
func init() {
	allocators["hyperelast1"] = func() Model { return new(HyperElast1) }
}

// Init initialises model
func (o *HyperElast1) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// parameters
	for _, p := range prms {
		switch p.N {
		case "kap":
			o.κ = p.V
		case "kapb":
			o.κb = p.V
		case "G0":
			o.G0 = p.V
		case "pr":
			o.pr = p.V
		case "pt":
			o.pt = p.V
		case "p0":
			o.p0 = p.V
		case "ev0":
			o.εv0 = p.V
		}
	}

	// derived
	o.pa = o.pr + o.pt
	o.a = 1.0 / o.κ
	return
}

// GetPrms gets (an example) of parameters
func (o HyperElast1) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "kap", V: 0.05},
		&fun.Prm{N: "kapb", V: 0.001},
		&fun.Prm{N: "G0", V: 10000},
		&fun.Prm{N: "pr", V: 2.0},
		&fun.Prm{N: "pt", V: 1.0},
		&fun.Prm{N: "p0", V: 1.5},
		&fun.Prm{N: "ev0", V: 0.0},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o HyperElast1) InitIntVars() (s *State, err error) {
	chk.Panic("HyperElast1: InitIntVars: not ready yet")
	return
}

// Update updates stresses for given strains
func (o *HyperElast1) Update(s *State, ε, Δε []float64) (err error) {
	chk.Panic("HyperElast1: Update: not ready yet")
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *HyperElast1) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	chk.Panic("HyperElast1: CalcD: not ready yet")
	return
}

// Invs computes p and q for given elastic εv and εd
func (o HyperElast1) Invs(εve, εde float64) (p, q float64) {
	pv := (o.pa + o.p0) * math.Exp(o.a*(o.εv0-εve))
	p = (1.0+1.5*o.a*o.κb*εde*εde)*pv - o.pa
	q = 3.0 * (o.G0 + o.κb*pv) * εde
	return
}

// Moduli computes the following derivatives:
//  Dvv = ∂²ψ/(∂εve ∂εve)
//  Dvd = ∂²ψ/(∂εve ∂εde)
//  Ddd = ∂²ψ/(∂εde ∂εde)
func (o HyperElast1) Moduli(εve, εde float64) (Dvv, Dvd, Ddd float64) {
	pv := (o.pa + o.p0) * math.Exp(o.a*(o.εv0-εve))
	Dvv = o.a * (1.0 + 1.5*o.a*o.κb*εde*εde) * pv
	Dvd = -3.0 * o.a * o.κb * εde * pv
	Ddd = 3.0 * (o.G0 + o.κb*pv)
	return
}

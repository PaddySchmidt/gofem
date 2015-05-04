// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/tsr"
)

// HyperElast1 implements a nonlinear hyperelastic model for powders and porous media
type HyperElast1 struct {

	// constants
	Nsig   int     // number of stress components
	EnoMin float64 // minimum value of ||dev(ε)||

	// parameters
	κ   float64 // κ
	κb  float64 // \bar{κ}
	G0  float64 // G0
	pr  float64 // pr
	pt  float64 // pt
	p0  float64 // p0
	εv0 float64 // εv0
	le  bool    // use linear elastic model
	K0  float64 // K0 (for linear model)

	// derived
	pa float64 // pa = pr + pt
	a  float64 // a = 1 / κ

	// auxiliary
	e []float64 // e = dev(ε)
}

// add model to factory
func init() {
	allocators["hyperelast1"] = func() Model { return new(HyperElast1) }
}

// Init initialises model
func (o *HyperElast1) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// constants
	o.Nsig = 2 * ndim
	o.EnoMin = 1e-8

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
		case "le":
			o.le = p.V > 0
		case "K0":
			o.K0 = p.V
		}
	}

	// derived
	o.pa = o.pr + o.pt
	o.a = 1.0 / o.κ

	// auxiliary
	o.e = make([]float64, 2*ndim)
	return
}

// Set_pt sets pt
func (o *HyperElast1) Set_pt(pt float64) {
	o.pt = pt
	o.pa = o.pr + o.pt
}

// Set_p0_ev0 sets p0 and εv0
func (o *HyperElast1) Set_p0_ev0(p0, εv0 float64) {
	o.p0 = p0
	o.εv0 = εv0
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
func (o HyperElast1) InitIntVars(σ []float64) (s *State, err error) {
	chk.Panic("HyperElast1: InitIntVars: not ready yet")
	return
}

// Update updates stresses for given strains
func (o *HyperElast1) Update(s *State, ε, dummy []float64) (err error) {
	eno, εv, εd := tsr.M_devε(o.e, ε)
	p, q := o.Calc_pq(εv, εd)
	if eno > o.EnoMin {
		for i := 0; i < o.Nsig; i++ {
			s.Sig[i] = -p*tsr.Im[i] + tsr.SQ2by3*q*o.e[i]/eno
		}
		return
	}
	for i := 0; i < o.Nsig; i++ {
		s.Sig[i] = -p * tsr.Im[i]
	}
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *HyperElast1) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	o.L_CalcD(D, s.EpsE)
	return
}

// principal strains /////////////////////////////////////////////////////////////////////////////

// Calc_pq computes p and q for given elastic εv and εd
func (o HyperElast1) Calc_pq(εv, εd float64) (p, q float64) {
	if o.le {
		p = o.p0 - o.K0*(εv-o.εv0)
		q = 3.0 * o.G0 * εd
		return
	}
	pv := (o.pa + o.p0) * math.Exp(o.a*(o.εv0-εv))
	p = (1.0+1.5*o.a*o.κb*εd*εd)*pv - o.pa
	q = 3.0 * (o.G0 + o.κb*pv) * εd
	return
}

// L_update computes principal stresses for given principal strains
func (o *HyperElast1) L_update(σ, ε []float64) (p, q float64) {
	eno, εv, εd := tsr.M_devε(o.e, ε) // using principal values since len(ε)=3
	p, q = o.Calc_pq(εv, εd)
	if eno > o.EnoMin {
		for i := 0; i < 3; i++ {
			σ[i] = -p*tsr.Im[i] + tsr.SQ2by3*q*o.e[i]/eno
		}
		return
	}
	for i := 0; i < 3; i++ {
		σ[i] = -p * tsr.Im[i]
	}
	return
}

// L_CalcD computes De in principal components for given principal elastic strains
//  D -- [ncp][ncp] elastic modulus
//  ε -- [ncp] elastic strains
//
//  Note: this method works also for non-principal components
//
func (o HyperElast1) L_CalcD(D [][]float64, ε []float64) {

	// number of components
	ncp := len(ε)

	// elastic modulus
	I, Psd := tsr.Im, tsr.Psd
	if o.le {
		for i := 0; i < ncp; i++ {
			for j := 0; j < ncp; j++ {
				D[i][j] = o.K0*I[i]*I[j] + 2.0*o.G0*Psd[i][j]
			}
		}
		return
	}

	// invariants of strain and normalised deviatoric direction
	eno, εv, εd := tsr.M_devε(o.e, ε)
	if eno > o.EnoMin {
		for i := 0; i < ncp; i++ {
			o.e[i] /= eno
		}
	} else {
		for i := 0; i < ncp; i++ {
			o.e[i] = 0
		}
	}

	// Dvv = ∂²ψ/(∂εve ∂εve)
	// Dvd = (∂²ψ/(∂εve ∂εde)) * sqrt(2/3)
	// Ddd2 = (∂²ψ/(∂εde ∂εde)) * 2 / 3
	pv := (o.pa + o.p0) * math.Exp(o.a*(o.εv0-εv))
	Dvv := o.a * (1.0 + 1.5*o.a*o.κb*εd*εd) * pv
	DvdS := -3.0 * o.a * o.κb * εd * pv * tsr.SQ2by3
	Ddd2 := 2.0 * (o.G0 + o.κb*pv)
	for i := 0; i < ncp; i++ {
		for j := 0; j < ncp; j++ {
			D[i][j] = Dvv*I[i]*I[j] + Ddd2*Psd[i][j] + DvdS*(I[i]*o.e[j]+o.e[i]*I[j])
		}
	}
}

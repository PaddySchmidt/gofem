// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

// CamClayMod implements the modified CamClay model
type CamClayMod struct {

	// basic data
	Nsig int            // number of σ and ε components
	CS   tsr.NcteM      // slope of cs line
	HE   HyperElast1    // hyper elasticity
	PU   PrincStrainsUp // stress updater

	// parameters
	λ   float64 // slope of isotropic compression model
	ocr float64 // initial over-consolidation ratio

	// auxiliary
	ch float64     // 1/(κ-λ)
	s  []float64   // dev(σ)
	n  []float64   // dM/dσ
	m  [][]float64 // d²M/(dσ dσ)
}

// add model to factory
func init() {
	allocators["ccm"] = func() Model { return new(CamClayMod) }
}

// Init initialises model
func (o *CamClayMod) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// basic data
	o.Nsig = 2 * ndim

	// parameters for CS model
	pp := []string{"φ", "Mfix"}
	vv := []float64{25, 1}
	for _, p := range prms {
		switch p.N {
		case "phi":
			vv[0] = p.V
		case "Mfix":
			vv[1] = p.V
		}
	}
	o.CS.Init(pp, vv)

	// parameters
	var pt float64
	for _, p := range prms {
		switch p.N {
		case "pt":
			pt = p.V
		case "c":
			pt = p.V / o.CS.Tanφ
		case "lam":
			o.λ = p.V
		case "ocr":
			o.ocr = p.V
		}
	}

	// parameters for HE model
	err = o.HE.Init(ndim, pstress, prms)
	if err != nil {
		return
	}
	o.HE.Set_pt(pt)

	// stress updater
	o.PU.Init(ndim, prms, o)

	// auxiliary
	o.ch = 1.0 / (o.HE.κ - o.λ)
	o.s = make([]float64, o.Nsig)
	o.n = make([]float64, o.Nsig)
	o.m = la.MatAlloc(o.Nsig, o.Nsig)
	return
}

// GetPrms gets (an example) of parameters
func (o *CamClayMod) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "phi", V: 25},
		&fun.Prm{N: "Mfix", V: 1},
		&fun.Prm{N: "c", V: 10},
		&fun.Prm{N: "lam", V: 0.1},
		&fun.Prm{N: "ocr", V: 1},
		&fun.Prm{N: "kap", V: 0.05},
		&fun.Prm{N: "kapb", V: 0},
		&fun.Prm{N: "G0", V: 10000},
		&fun.Prm{N: "pr", V: 1.0},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o *CamClayMod) InitIntVars(σ []float64) (s *State, err error) {

	// compute α0
	p, q, w := tsr.M_pqw(σ)
	M := o.CS.M(w)
	pt := o.HE.pt
	var α0 float64
	if math.Abs(p+pt) < 1e-8 {
		α0 = 1e-8
	} else {
		α0 = p + q*q/(M*M*(p+pt))
	}

	// set state
	nalp := 1 // alp[0] = α0 (yield surface size controller)
	s = NewState(o.Nsig, nalp, false, true)
	copy(s.Sig, σ)
	s.Alp[0] = α0 * o.ocr

	// compute initial strains
	o.HE.CalcEps0(s)
	return
}

// Update updates stresses for given strains
func (o *CamClayMod) Update(s *State, ε, Δε []float64, eid, ipid int) (err error) {
	return o.PU.Update(s, ε, Δε, eid, ipid)
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *CamClayMod) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	return o.PU.CalcD(D, s)
}

// ContD computes D = dσ_new/dε_new continuous
func (o *CamClayMod) ContD(D [][]float64, s *State) (err error) {
	chk.Panic("CCM: ContD is not available")
	return
}

// EPmodel ///////////////////////////////////////////////////////////////////////////////////////////

// Info returns some information and data from this model
func (o *CamClayMod) Info() (nalp, nsurf int) {
	return 1, 1
}

// Get_phi gets φ or returns 0
func (o *CamClayMod) Get_phi() float64 { return 0 }

// Get_bsmp gets b coefficient if using SMP invariants
func (o *CamClayMod) Get_bsmp() float64 { return 0 }

// Set_bsmp sets b coefficient if using SMP invariants
func (o *CamClayMod) Set_bsmp(b float64) {}

// L_YieldFunc computes the yield function value for given principal stresses (σ)
func (o *CamClayMod) L_YieldFunc(σ, α []float64) float64 {
	p, q, w := tsr.M_pqw(σ)
	M := o.CS.M(w)
	pt := o.HE.pt
	n0 := (p + pt) * (p - α[0])
	return q*q + M*M*n0
}

// YieldFuncs computes yield function values
func (o *CamClayMod) YieldFuncs(s *State) []float64 {
	p, q, w := tsr.M_pqw(s.Sig)
	M := o.CS.M(w)
	pt := o.HE.pt
	α0 := s.Alp[0]
	n0 := (p + pt) * (p - α0)
	return []float64{q*q + M*M*n0}
}

// ElastUpdate updates state with an elastic response
func (o *CamClayMod) ElastUpdate(s *State, ε []float64) {
	o.HE.Update(s, ε, nil, 0, 0)
}

// ElastD returns continuum elastic D
func (o *CamClayMod) ElastD(D [][]float64, s *State) {
	o.HE.CalcD(D, s, false)
}

// E_CalcSig computes principal stresses for given principal elastic strains
func (o *CamClayMod) E_CalcSig(σ, εe []float64) {
	o.HE.L_update(σ, εe)
}

// E_CalcDe computes elastic modulus in principal components
func (o *CamClayMod) E_CalcDe(De [][]float64, εe []float64) {
	o.HE.L_CalcD(De, εe)
}

// L_FlowHard computes model variabes for given principal values
func (o *CamClayMod) L_FlowHard(Nb, h, σ, α []float64) (f float64, err error) {
	p, q, w := tsr.M_pqws(o.s, σ)
	M := o.CS.M(w)
	pt := o.HE.pt
	n0 := (p + pt) * (p - α[0])
	n1 := 2.0*p + pt - α[0]
	I := tsr.Im
	for i := 0; i < 3; i++ {
		Nb[i] = 3.0*o.s[i] - M*M*n1*I[i]/3.0 + 2.0*M*n0*o.n[i]
	}
	trNb := Nb[0] + Nb[1] + Nb[2]
	h[0] = o.ch * (o.HE.pa + α[0]) * trNb
	f = q*q + M*M*n0
	return
}

// L_SecondDerivs computes second order derivatives
//  N    -- ∂f/∂σ     [nsig]
//  Nb   -- ∂g/∂σ     [nsig]
//  A    -- ∂f/∂α_i   [nalp]
//  h    -- hardening [nalp]
//  Mb   -- ∂Nb/∂εe   [nsig][nsig]
//  a_i  -- ∂Nb/∂α_i  [nalp][nsig]
//  b_i  -- ∂h_i/∂εe  [nalp][nsig]
//  c_ij -- ∂h_i/∂α_j [nalp][nalp]
func (o *CamClayMod) L_SecondDerivs(N, Nb, A, h []float64, Mb, a, b, c [][]float64, σ, α []float64) (err error) {
	p, _, w := tsr.M_pqws(o.s, σ)
	M := o.CS.M(w)
	pt := o.HE.pt
	n0 := (p + pt) * (p - α[0])
	n1 := 2.0*p + pt - α[0]
	I := tsr.Im
	for i := 0; i < 3; i++ {
		Nb[i] = 3.0*o.s[i] - M*M*n1*I[i]/3.0 + 2.0*M*n0*o.n[i]
		N[i] = Nb[i]
	}
	d0 := 2.0 * M * M / 9.0
	d1 := -2.0 * M * n1 / 3.0
	d2 := 2.0 * n0
	d3 := 2.0 * M * n0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			Mb[i][j] = 3.0*tsr.Psd[i][j] + d0*I[i]*I[j] + d1*(I[i]*o.n[j]+o.n[i]*I[j]) + d2*o.n[i]*o.n[j] + d3*o.m[i][j]
		}
		a[0][i] = M*M*I[i]/3.0 - 2.0*M*(p+pt)*o.n[i]
		b[0][i] = o.ch * (o.HE.pa + α[0]) * M * (2.0*M*I[i]/3.0 - 2.0*n1*o.n[i])
	}
	trNb := Nb[0] + Nb[1] + Nb[2]
	A[0] = -M * M * (p + pt)
	h[0] = o.ch * (o.HE.pa + α[0]) * trNb
	c[0][0] = o.ch*trNb + o.ch*(o.HE.pa+α[0])*M*M
	return
}

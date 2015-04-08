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

// CamClayMod implements the modified CamClay model
type CamClayMod struct {

	// basic data
	Nsig int            // number of σ and ε components
	CS   tsr.NcteM      // slope of cs line
	HE   HyperElast1    // hyper elasticity
	PU   PrincStrainsUp // stress updater
	Lσ   []float64      // principal stresses

	// parameters
	λ   float64 // slope of isotropic compression model
	ocr float64 // initial over-consolidation ratio

	// auxiliary
	ch   float64   // 1/(κ-λ)
	devσ []float64 // dev(σ)
	nvec []float64 // dM/dσ
}

// add model to factory
func init() {
	allocators["ccm"] = func() Model { return new(CamClayMod) }
}

// Init initialises model
func (o *CamClayMod) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// basic data
	o.Nsig = 2 * ndim
	o.Lσ = make([]float64, 3)

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
	o.devσ = make([]float64, 2*ndim)
	o.nvec = make([]float64, 2*ndim)
	return
}

// GetPrms gets (an example) of parameters
func (o CamClayMod) GetPrms() fun.Prms {
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
		&fun.Prm{N: "p0", V: 0.0},
		&fun.Prm{N: "ev0", V: 0.0},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o CamClayMod) InitIntVars(σ []float64) (s *State, err error) {

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

	// set HE model
	o.HE.Set_p0_ev0(p, 0)

	// set state
	nalp := 1      // alp[0] = α0 (yield surface size controller)
	nphi := o.Nsig // phi[0] = εe (elastic strains)
	s = NewState(o.Nsig, nalp, nphi, false)
	copy(s.Sig, σ)
	s.Alp[0] = α0 * o.ocr
	return
}

// Update updates stresses for given strains
func (o *CamClayMod) Update(s *State, ε, Δε []float64) (err error) {
	return o.PU.Update(s, ε, Δε)
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *CamClayMod) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	return
}

// ContD computes D = dσ_new/dε_new continuous
func (o *CamClayMod) ContD(D [][]float64, s *State) (err error) {
	chk.Panic("CCM: ContD is not available")
	return
}

// EPmodel ///////////////////////////////////////////////////////////////////////////////////////////

// Info returns some information and data from this model
func (o CamClayMod) Info() (nalp, nsurf int, fcoef, pt, pr float64) {
	return 1, 1, o.HE.pr * o.HE.pr, o.HE.pt, o.HE.pr
}

// IsoF returns the isotropic function, if any
func (o CamClayMod) IsoF() *tsr.IsoFun {
	return nil
}

// YieldFuncs computes yield function values
func (o CamClayMod) YieldFuncs(s *State) []float64 {
	p, q, w := tsr.M_pqw(s.Sig)
	M := o.CS.M(w)
	pt := o.HE.pt
	α0 := s.Alp[0]
	n0 := (p + pt) * (p - α0)
	return []float64{q*q + M*M*n0}
}

// ElastUpdate updates state with an elastic response
func (o CamClayMod) ElastUpdate(s *State, ε, Δε []float64) {
	εe := s.Phi
	for i := 0; i < o.Nsig; i++ {
		εe[i] += Δε[i]
	}
	o.HE.Update(s, εe, Δε)
}

// PVE_CalcSig computes principal stresses for given principal elastic strains
func (o CamClayMod) PVE_CalcSig(σ, εe []float64) {
	o.HE.L_update(σ, εe)
}

// PVE_FlowHard computes model variabes for given elastic strains (principal values)
func (o CamClayMod) PVE_FlowHard(Nb, h, σ, α []float64) (f float64, err error) {
	p, q, w := tsr.M_pqws(o.devσ, σ)
	M := o.CS.M(w)
	pt := o.HE.pt
	n0 := (p + pt) * (p - α[0])
	n1 := 2.0*p + pt - α[0]
	for i := 0; i < 3; i++ {
		Nb[i] = 3.0*o.devσ[i] - M*M*n1*tsr.Im[i]/3.0 + 2.0*M*n0*o.nvec[i]
	}
	trNb := Nb[0] + Nb[1] + Nb[2]
	h[0] = o.ch * (o.HE.pa + α[0]) * trNb
	f = q*q + M*M*n0
	//io.Pforan("Nb = %v\n", Nb)
	//io.Pforan("h0 = %v\n", h[0])
	return
}

func (o CamClayMod) PVE_SecondDerivs(Nb, h, σ, α []float64) (f float64, err error) {
	return
}

/*
func (o CamClayMod) SecondDerivs(dAdσ, dAdα, dNdσ, dNdα, dNbdσ, dNbdα, dhdσ, dhdα [][]float64, α []float64) (err error) {
	n0 := (p + pt) * (p - α[0])
	n1 := 2.0*p + pt - α[0]
	M := o.Calc.M(v.w)
	M2 := M * M
	d0 := 2.0 * M2 / 9.0
	d1 := -2.0 * M * n1 / 3.0
	d2 := 2.0 * n0
	d3 := 2.0 * M * n0
	o.CS.Deriv2(o.d2Mdσdσ, o.dMdσ, v.σ, v.s, p, v.q, v.w)
	for i := 0; i < o.nσ; i++ {
		for j := 0; j < o.nσ; j++ {
			dNdσ[i][j] = 3.0*tsr.Psd[i][j] + d0*tsr.Im[i]*tsr.Im[j] + d1*(tsr.Im[i]*o.dMdσ[j]+o.dMdσ[i]*tsr.Im[j]) + d2*o.dMdσ[i]*o.dMdσ[j] + d3*o.d2Mdσdσ[i][j]
			dNbdσ[i][j] = dNdσ[i][j]
		}
		b1_i := 2.0*M2*tsr.Im[i]/3.0 - 2.0*M*n1*o.dMdσ[i]
		dAdσ[0][i] = M2*tsr.Im[i]/3.0 - 2.0*M*(p+pt)*o.dMdσ[i]
		dNdα[i][0] = dAdσ[0][i]
		dNbdα[i][0] = dNdα[i][0]
		dhdσ[0][i] = o.ψ * (o.pr + α[0]) * b1_i
	}
	trNb := -M2 * n1
	b0 := M2
	dAdα[0][0] = 0
	dhdα[0][0] = o.ψ*trNb + o.ψ*(o.pr+α[0])*b0
}
*/

// PVE_LoadSet sets state with new data (principal strains) from elastoplastic loading
//func (o CamClayMod) PVE_LoadSet(s *State, Δγ float64, εe, α []float64, P [][]float64) (err error) {
//s.Dgam = Δγ
//s.Loading = true
//copy(s.Alp, α)
//o.HE.L_update(o.Lσ, εe)
//return
//}

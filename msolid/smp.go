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

// SmpInvs implements a model with SMP invariants similar to Drucker-Prager model
type SmpInvs struct {

	// basic data
	Nsig int            // number of σ and ε components
	HE   HyperElast1    // hyper elasticity
	PU   PrincStrainsUp // stress updater
	Isof tsr.IsoFun     // isotropic function structure
	RM   fun.RefDecSp1  // reference model for smoothing

	// parameters
	r    float64 // radius
	rtyp int     // rounding type
	pe   float64 // large p value for rounding with Bezier

	// derived parameters
	M  float64 // slope of yield line in a p-q graph
	sα float64 // sin(α) with α=atan(M) (for rounding with circle)
	cα float64 // cos(α) with α=atan(M) (for rounding with circle)
}

// add model to factory
func init() {
	allocators["smp"] = func() Model { return new(SmpInvs) }
}

// Init initialises model
func (o *SmpInvs) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// basic data
	o.Nsig = 2 * ndim

	// parameters
	var c, φ, a, b, β, ϵ float64
	βrm := 20.0 // β_{reference-model}
	for _, p := range prms {
		switch p.N {
		case "c":
			c = p.V
		case "phi":
			φ = p.V
		case "r":
			o.r = p.V
		case "a":
			a = p.V
		case "b":
			b = p.V
		case "bet":
			β = p.V
		case "eps":
			ϵ = p.V
		case "betrm":
			βrm = p.V
		case "pe":
			o.pe = p.V
		case "rtyp":
			o.rtyp = int(p.V)
		}
	}

	// auxiliary
	o.M = tsr.NewSmpCalcμ(φ, a, b, β, ϵ)
	α := math.Atan(o.M)
	o.sα = math.Sin(α)
	o.cα = math.Cos(α)
	shift := c / math.Tan(φ*math.Pi/180.0)

	// parameters for HE model
	err = o.HE.Init(ndim, pstress, prms)
	if err != nil {
		return
	}
	o.HE.Set_pt(shift)

	// stress updater
	o.PU.Init(ndim, prms, o)

	// rounding
	switch o.rtyp {
	case 1: // circle
		m := o.r/o.sα - o.r
		if m > shift {
			shift = m
		}
	case 2: // reference model
		m := o.r/o.sα - o.r
		if m > shift {
			shift = m
		}
		pc := o.r / o.sα
		pb := pc - o.r
		pa := 0.0
		o.RM.Init([]*fun.Prm{
			&fun.Prm{N: "bet", V: βrm},
			&fun.Prm{N: "lam1", V: 1.0 / o.M},
			&fun.Prm{N: "ya", V: -pa},
			&fun.Prm{N: "yb", V: -pb},
		})
	case 3: // o2 Bezier
		pa := 0.0
		pb := o.r
		m := pb - pa
		if m > shift {
			shift = m
		}
		pemin := 2.0*pb - pa + 1e-7
		if o.pe < pemin {
			chk.Panic("pe(%g) must be greater than %g", o.pe, pemin)
		}
	}

	// set isotropic function
	o.Isof.Init(a, b, β, ϵ, shift, o.Nsig, o.ffcn, o.gfcn, o.hfcn)
	return
}

// GetPrms gets (an example) of parameters
func (o SmpInvs) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "phi", V: 25},
		&fun.Prm{N: "Mfix", V: 1},
		&fun.Prm{N: "c", V: 10},
		&fun.Prm{N: "K0", V: 10000},
		&fun.Prm{N: "G0", V: 10000},
		&fun.Prm{N: "pr", V: 1.0},
		&fun.Prm{N: "p0", V: 0.0},
		&fun.Prm{N: "ev0", V: 0.0},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o SmpInvs) InitIntVars(σ []float64) (s *State, err error) {

	// set HE model
	p := tsr.M_p(σ)
	o.HE.Set_p0_ev0(p, 0)

	// set state
	nalp := 1 // alp[0] = εpb (cumulated plastic strain)
	s = NewState(o.Nsig, nalp, false)
	copy(s.Sig, σ)
	s.Alp[0] = 0
	return
}

// Update updates stresses for given strains
func (o *SmpInvs) Update(s *State, ε, Δε []float64) (err error) {
	return o.PU.Update(s, ε, Δε)
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *SmpInvs) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	return o.PU.CalcD(D, s)
}

// ContD computes D = dσ_new/dε_new continuous
func (o *SmpInvs) ContD(D [][]float64, s *State) (err error) {
	chk.Panic("Smp: ContD is not available")
	return
}

// EPmodel ///////////////////////////////////////////////////////////////////////////////////////////

// Info returns some information and data from this model
func (o SmpInvs) Info() (nalp, nsurf int, fcoef, pt, pr float64) {
	return 1, 1, o.HE.pr, o.HE.pt, o.HE.pr
}

// IsoF returns the isotropic function
func (o SmpInvs) IsoF() *tsr.IsoFun {
	return &o.Isof
}

// YieldFuncs computes yield function values
func (o SmpInvs) YieldFuncs(s *State) []float64 {
	res, err := o.Isof.Fa(s.Sig, s.Alp[0])
	if err != nil {
		chk.Panic("cannot compute isotrpic function Fa: %v", err)
	}
	return []float64{res}
}

// ElastUpdate updates state with an elastic response
func (o SmpInvs) ElastUpdate(s *State, ε, Δε []float64) {
	o.HE.Update(s, s.EpsE, nil)
}

// ElastD returns continuum elastic D
func (o SmpInvs) ElastD(D [][]float64, s *State) {
	o.HE.CalcD(D, s, false)
}

// E_CalcSig computes principal stresses for given principal elastic strains
func (o SmpInvs) E_CalcSig(σ, εe []float64) {
	o.HE.L_update(σ, εe)
}

// E_CalcDe computes elastic modulus in principal components
func (o SmpInvs) E_CalcDe(De [][]float64, εe []float64) {
	o.HE.L_CalcD(De, εe)
}

// L_FlowHard computes model variabes for given principal values
func (o SmpInvs) L_FlowHard(Nb, h, σ, α []float64) (f float64, err error) {
	f, err = o.Isof.Gp(σ, α)
	if err != nil {
		return
	}
	for i := 0; i < 3; i++ {
		Nb[i] = o.Isof.Dfdλ[i]
	}
	// no hardening yet
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
func (o SmpInvs) L_SecondDerivs(N, Nb, A, h []float64, Mb, a, b, c [][]float64, σ, α []float64) (err error) {
	_, err = o.Isof.Gp(σ, α)
	if err != nil {
		return
	}
	err = o.Isof.HafterGp(σ)
	if err != nil {
		return
	}
	for i := 0; i < 3; i++ {
		N[i] = o.Isof.Dfdλ[i]
		Nb[i] = N[i]
		for j := 0; j < 3; j++ {
			Mb[i][j] = o.Isof.Dgdλ[i][j]
		}
	}
	// no hardening yet
	return
}

// auxiliary functions ///////////////////////////////////////////////////////////////////////////////

func (o SmpInvs) auxvars(p, q float64) (pc, R, pd float64) {
	pc = o.r / o.sα
	R = math.Sqrt(q*q + (p-pc)*(p-pc))
	pd = pc - R*o.sα
	return
}

func (o SmpInvs) auxvarsB(p, q float64) (pb, m, n, s, t float64) {
	pa := 0.0
	pb = o.r
	m = pb - pa
	n = o.pe - pb
	s = math.Sqrt(o.M * (m*m*o.M + (n-m)*q))
	t = (s - m*o.M) / ((n - m) * o.M)
	return
}

func (o SmpInvs) ffcn(p, q float64, args ...interface{}) (fval float64) {
	switch o.rtyp {
	case 0: // straight line
		fval = q - o.M*p
	case 1: // circle
		pc, R, pd := o.auxvars(p, q)
		if p < pd {
			fval = R - o.r
		} else {
			fval = q*o.cα - (p-pc)*o.sα - o.r
		}
	case 2: // reference model
		fval = -p - o.RM.F(q, nil)
	case 3: // o2 Bezier
		pb, _, n, _, t := o.auxvarsB(p, q)
		fval = n*t*t + pb - p
	}
	return
}

func (o SmpInvs) gfcn(p, q float64, args ...interface{}) (dfdp, dfdq float64) {
	switch o.rtyp {
	case 0: // straight line
		dfdp = -o.M
		dfdq = 1.0
	case 1: // circle
		pc, R, pd := o.auxvars(p, q)
		if p < pd {
			dfdp = (p - pc) / R
			dfdq = q / R
		} else {
			dfdp = -o.sα
			dfdq = o.cα
		}
	case 2: // reference model
		dfdp = -1.0
		dfdq = -o.RM.G(q, nil)
	case 3: // o2 Bezier
		_, _, n, s, t := o.auxvarsB(p, q)
		dfdp = -1.0
		dfdq = n * t / s
	}
	return
}

func (o SmpInvs) hfcn(p, q float64, args ...interface{}) (d2fdp2, d2fdq2, d2fdpdq float64) {
	switch o.rtyp {
	case 0: // straight line
	case 1: // circle
		pc, R, pd := o.auxvars(p, q)
		if p < pd {
			R3 := R * R * R
			d2fdp2 = 1.0/R - (p-pc)*(p-pc)/R3
			d2fdq2 = 1.0/R - q*q/R3
			d2fdpdq = -(p - pc) * q / R3
		}
	case 2: // reference model
		d2fdq2 = -o.RM.H(q, nil)
	case 3: // o2 Bezier
		_, m, n, s, t := o.auxvarsB(p, q)
		s3 := s * s * s
		d2fdq2 = (n*s - n*t*(n-m)*o.M) / (2.0 * s3)
	}
	return
}

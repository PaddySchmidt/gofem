// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/tsr"
)

// DruckerPrager implements Drucker-Prager plasticity model
type DruckerPrager struct {
	SmallElasticity
	M   float64   // slope of fc line
	Mb  float64   // slope of fc line of plastic potential
	qy0 float64   // initial qy
	H   float64   // hardening variable
	ten []float64 // auxiliary tensor
}

// add model to factory
func init() {
	allocators["dp"] = func() Model { return new(DruckerPrager) }
}

// Init initialises model
func (o *DruckerPrager) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// parse parameters
	err = o.SmallElasticity.Init(ndim, pstress, prms)
	if err != nil {
		return
	}
	var c, φ float64
	var typ int
	for _, p := range prms {
		switch p.N {
		case "M":
			o.M = p.V
		case "Mb":
			o.Mb = p.V
		case "qy0":
			o.qy0 = p.V
		case "H":
			o.H = p.V
		case "c":
			c = p.V
		case "phi":
			φ = p.V
		case "typ":
			typ = int(p.V)
		case "E", "nu", "l", "G", "K", "rho":
		default:
			return chk.Err("dp: parameter named %q is incorrect\n", p.N)
		}
	}

	// compute M from φ
	//  typ == 0 : compression cone (outer)
	//      == 1 : extension cone (inner)
	//      == 2 : plane-strain
	if φ > 0 {
		o.M, o.qy0, err = Mmatch(c, φ, typ)
		if err != nil {
			return
		}
		o.Mb = o.M
	}
	//io.Pforan("E=%v nu=%v\n", o.E, o.Nu)
	//io.Pforan("c=%v phi=%v M=%v qy0=%v\n", c, φ, o.M, o.qy0)

	// auxiliary structures
	o.ten = make([]float64, o.Nsig)
	return
}

// GetPrms gets (an example) of parameters
func (o DruckerPrager) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "M", V: 1},
		&fun.Prm{N: "Mb", V: 1},
		&fun.Prm{N: "qy0", V: 0.5},
		&fun.Prm{N: "H", V: 0},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o DruckerPrager) InitIntVars(σ []float64) (s *State, err error) {
	s = NewState(o.Nsig, 1, false, false)
	copy(s.Sig, σ)
	return
}

// Update updates stresses for given strains
func (o *DruckerPrager) Update(s *State, ε, Δε []float64, eid, ipid int) (err error) {

	// set flags
	s.Loading = false    // => not elastoplastic
	s.ApexReturn = false // => not return-to-apex
	s.Dgam = 0           // Δγ := 0

	// accessors
	σ := s.Sig
	α0 := &s.Alp[0]

	// copy of α0 at beginning of step
	α0ini := *α0

	// trial stress
	var devΔε_i float64
	trΔε := Δε[0] + Δε[1] + Δε[2]
	for i := 0; i < o.Nsig; i++ {
		devΔε_i = Δε[i] - trΔε*tsr.Im[i]/3.0
		o.ten[i] = σ[i] + o.K*trΔε*tsr.Im[i] + 2.0*o.G*devΔε_i // ten := σtr
	}
	ptr, qtr := tsr.M_p(o.ten), tsr.M_q(o.ten)

	// trial yield function
	ftr := qtr - o.M*ptr - o.qy0 - o.H*(*α0)

	// elastic update
	if ftr <= 0.0 {
		copy(σ, o.ten) // σ := ten = σtr
		return
	}

	// elastoplastic update
	var str_i float64
	hp := 3.0*o.G + o.K*o.M*o.Mb + o.H
	s.Dgam = ftr / hp
	*α0 += s.Dgam
	pnew := ptr + s.Dgam*o.K*o.Mb
	m := 1.0 - s.Dgam*3.0*o.G/qtr
	for i := 0; i < o.Nsig; i++ {
		str_i = o.ten[i] + ptr*tsr.Im[i]
		σ[i] = m*str_i - pnew*tsr.Im[i]
	}
	s.Loading = true

	// check for apex singularity
	acone := qtr - s.Dgam*3.0*o.G
	if acone < 0 {
		s.Dgam = (-o.M*ptr - o.qy0 - o.H*α0ini) / (3.0*o.K*o.M + o.H)
		*α0 = α0ini + s.Dgam
		pnew = ptr + s.Dgam*3.0*o.K
		for i := 0; i < o.Nsig; i++ {
			σ[i] = -pnew * tsr.Im[i]
		}
		s.ApexReturn = true
	}
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *DruckerPrager) CalcD(D [][]float64, s *State, firstIt bool) (err error) {

	// set first Δγ
	if firstIt {
		s.Dgam = 0
	}

	// elastic
	if !s.Loading {
		return o.SmallElasticity.CalcD(D, s)
	}

	// return to apex
	if s.ApexReturn {
		a1 := o.K * o.H / (3.0*o.K*o.M + o.H)
		for i := 0; i < o.Nsig; i++ {
			for j := 0; j < o.Nsig; j++ {
				D[i][j] = a1 * tsr.Im[i] * tsr.Im[j]
			}
		}
		return
	}

	// elastoplastic => consistent stiffness
	σ := s.Sig
	Δγ := s.Dgam
	p, q := tsr.M_p(σ), tsr.M_q(σ)
	qtr := q + Δγ*3.0*o.G
	m := 1.0 - Δγ*3.0*o.G/qtr
	nstr := tsr.SQ2by3 * qtr // norm(str)
	for i := 0; i < o.Nsig; i++ {
		o.ten[i] = (σ[i] + p*tsr.Im[i]) / (m * nstr) // ten := unit(str) = snew / (m * nstr)
	}
	hp := 3.0*o.G + o.K*o.M*o.Mb + o.H
	a1 := o.K - o.K*o.K*o.Mb*o.M/hp
	a2 := -2.0 * o.G * o.K * o.Mb * tsr.SQ3by2 / hp
	b1 := -tsr.SQ6 * o.G * o.M * o.K / hp
	b2 := 6.0 * o.G * o.G * (Δγ/qtr - 1.0/hp)
	for i := 0; i < o.Nsig; i++ {
		for j := 0; j < o.Nsig; j++ {
			D[i][j] = 2.0*o.G*m*tsr.Psd[i][j] +
				a1*tsr.Im[i]*tsr.Im[j] +
				a2*tsr.Im[i]*o.ten[j] +
				b1*o.ten[i]*tsr.Im[j] +
				b2*o.ten[i]*o.ten[j]
		}
	}
	return
}

// ContD computes D = dσ_new/dε_new continuous
func (o *DruckerPrager) ContD(D [][]float64, s *State) (err error) {

	// elastic part
	err = o.SmallElasticity.CalcD(D, s)
	if err != nil {
		return
	}

	// only elastic
	if !s.Loading {
		return
	}

	// elastoplastic
	σ := s.Sig
	d1 := o.K*o.Mb*o.M + 3.0*o.G + o.H
	a1 := o.K * o.K * o.Mb * o.M / d1
	a2 := tsr.SQ6 * o.K * o.G * o.Mb / d1
	a3 := tsr.SQ6 * o.K * o.G * o.M / d1
	a4 := 6.0 * o.G * o.G / d1
	sno, _, _ := tsr.M_devσ(o.ten, σ) // ten := dev(σ)
	for i := 0; i < o.Nsig; i++ {
		for j := 0; j < o.Nsig; j++ {
			D[i][j] -= a1*tsr.Im[i]*tsr.Im[j] +
				a2*tsr.Im[i]*o.ten[j]/sno +
				a3*o.ten[i]*tsr.Im[j]/sno +
				a4*o.ten[i]*o.ten[j]/(sno*sno)
		}
	}
	return
}

// EPmodel ///////////////////////////////////////////////////////////////////////////////////////////

// Info returns some information and data from this model
func (o DruckerPrager) Info() (nalp, nsurf int) {
	return 1, 1
}

// Get_phi gets φ or returns 0
func (o DruckerPrager) Get_phi() float64 { return 0 }

// Get_bsmp gets b coefficient if using SMP invariants
func (o DruckerPrager) Get_bsmp() float64 { return 0 }

// Set_bsmp sets b coefficient if using SMP invariants
func (o *DruckerPrager) Set_bsmp(b float64) {}

// L_YieldFunc computes the yield function value for given principal stresses (σ)
func (o *DruckerPrager) L_YieldFunc(σ, α []float64) float64 {
	chk.Panic("DruckerPrager: L_YieldFunc is not implemented yet")
	return 0
}

// YieldFs computes the yield functions
func (o DruckerPrager) YieldFuncs(s *State) []float64 {
	p, q := tsr.M_p(s.Sig), tsr.M_q(s.Sig)
	α0 := s.Alp[0]
	return []float64{q - o.M*p - o.qy0 - o.H*α0}
}

// ElastUpdate updates state with an elastic response
func (o DruckerPrager) ElastUpdate(s *State, ε []float64) {
	var devε_i float64
	trε := ε[0] + ε[1] + ε[2]
	for i := 0; i < o.Nsig; i++ {
		devε_i = ε[i] - trε*tsr.Im[i]/3.0
		s.Sig[i] = o.K*trε*tsr.Im[i] + 2.0*o.G*devε_i
	}
}

// ElastD returns continuum elastic D
func (o DruckerPrager) ElastD(D [][]float64, s *State) {
}

// E_CalcSig computes principal stresses for given principal elastic strains
func (o DruckerPrager) E_CalcSig(σ, εe []float64) {
}

// E_CalcDe computes elastic modulus in principal components
func (o DruckerPrager) E_CalcDe(De [][]float64, εe []float64) {
}

// L_FlowHard computes model variabes for given principal values
func (o DruckerPrager) L_FlowHard(Nb, h, σ, α []float64) (f float64, err error) {
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
func (o DruckerPrager) L_SecondDerivs(N, Nb, A, h []float64, Mb, a, b, c [][]float64, σ, α []float64) (err error) {
	return
}

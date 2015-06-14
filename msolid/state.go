// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/la"

// State holds all continuum mechanics data, including for updating the state
type State struct {

	// essential
	Sig  []float64 // σ: current Cauchy stress tensor (effective) [nsig]
	Eps0 []float64 // ε0: initial strains

	// for plasticity (if len(α) > 0)
	EpsE       []float64 // elastic strain
	EpsTr      []float64 // trial elastic strain
	Alp        []float64 // α: internal variables of rate type [nalp]
	Dgam       float64   // Δγ: increment of Lagrange multiplier (for plasticity only)
	Loading    bool      // unloading flag (for plasticity only)
	ApexReturn bool      // return-to-apex (for plasticity only)

	// for large deformations
	F [][]float64 // deformation gradient [3][3]
}

// NewState allocates state structure for small or large deformation analyses
//  large -- large deformation analyses; otherwise small strains
//  nle   -- non-linear elastic
func NewState(nsig, nalp int, large, nle bool) *State {

	// essential
	var state State
	state.Sig = make([]float64, nsig)
	state.Eps0 = make([]float64, nsig)

	// for plasticity
	if nalp > 0 {
		state.EpsTr = make([]float64, nsig)
		state.Alp = make([]float64, nalp)
	}

	// non-linear elasticity
	if nalp > 0 || nle {
		state.EpsE = make([]float64, nsig)
	}

	// large deformations
	if large {
		state.F = la.MatAlloc(3, 3)
	}
	return &state
}

// Set copies states
//  Note: 1) this and other states must have been pre-allocated with the same sizes
//        2) this method does not check for errors
func (o *State) Set(other *State) {

	// essential
	copy(o.Sig, other.Sig)
	copy(o.Eps0, other.Eps0)

	// for plasticity
	if len(o.Alp) > 0 {
		copy(o.EpsTr, other.EpsTr)
		copy(o.Alp, other.Alp)
		o.Dgam = other.Dgam
		o.Loading = other.Loading
		o.ApexReturn = other.ApexReturn
	}

	// non-linear elasticity
	if len(o.EpsE) > 0 {
		copy(o.EpsE, other.EpsE)
	}

	// large deformations
	if len(o.F) > 0 {
		la.MatCopy(o.F, 1, other.F)
	}
}

// GetCopy returns a copy of this state
func (o *State) GetCopy() *State {
	large := len(o.F) > 0
	other := NewState(len(o.Sig), len(o.Alp), large, len(o.EpsE) > 0)
	other.Set(o)
	return other
}

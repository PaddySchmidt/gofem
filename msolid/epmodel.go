// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/tsr"

// EPmodel implements an elasto-plastic model
//  PVS -- principal values formulation with given principal stresses
//  PVE -- principal values formulation with given principal elastic strains
type EPmodel interface {
	Model
	Small

	Info() (nalp, nsurf int)           // Info returns some information and data from this model
	IsoF() *tsr.IsoFun                 // IsoF returns the isotropic function, if any
	YieldFuncs(s *State) []float64     // YieldFs computes the yield functions
	ElastUpdate(s *State, ε []float64) // ElastUpdate updates state with an elastic response
	ElastD(D [][]float64, s *State)    // ElastD returns continuum elastic D

	// E_CalcSig computes principal stresses for given principal elastic strains
	E_CalcSig(σ, εe []float64)

	// E_CalcDe computes elastic modulus in principal components
	E_CalcDe(De [][]float64, εe []float64)

	// L_FlowHard computes model variabes for given principal values
	L_FlowHard(Nb, h, εe, α []float64) (f float64, err error)

	// L_SecondDerivs computes second order derivatives in principal values
	//  N    -- ∂f/∂σ     [nsig]
	//  Nb   -- ∂g/∂σ     [nsig]
	//  A    -- ∂f/∂α_i   [nalp]
	//  h    -- hardening [nalp]
	//  Mb   -- ∂Nb/∂εe   [nsig][nsig]
	//  a_i  -- ∂Nb/∂α_i  [nalp][nsig]
	//  b_i  -- ∂h_i/∂εe  [nalp][nsig]
	//  c_ij -- ∂h_i/∂α_j [nalp][nalp]
	L_SecondDerivs(N, Nb, A, h []float64, Mb, a, b, c [][]float64, σ, α []float64) (err error)
}

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

	Info() (nalp, nsurf int, fcoef, pt, pr float64) // Info returns some information and data from this model
	IsoF() *tsr.IsoFun                              // IsoF returns the isotropic function, if any
	YieldFuncs(s *State) []float64                  // YieldFs computes the yield functions
	ElastUpdate(s *State, ε, Δε []float64)          // ElastUpdate updates state with an elastic response

	// PVE_CalcSig computes principal stresses for given principal elastic strains
	PVE_CalcSig(σ, εe []float64)

	// PVE_FlowHard computes model variabes for given elastic strains (principal values)
	PVE_FlowHard(Nb, h, εe, α []float64) (f float64, err error)

	// PVE_LoadSet sets state with new data (principal strains) from elastoplastic loading
	//PVE_LoadSet(s *State, Δγ float64, εe, α []float64, P [][]float64) (err error)
}

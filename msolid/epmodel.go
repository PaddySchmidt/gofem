// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/tsr"

type EPmodel interface {
	Model
	Small
	Info() (nsurf int, fcoef, pt, pr float64) // Info returns some information and data from this model
	IsoF() *tsr.IsoFun                        // IsoF returns the isotropic function, if any
	YieldFuncs(sta *State) []float64          // YieldFs computes the yield functions
	ElastUpdate(sta *State, Δε []float64)     // ElastUpdate updates state with an elastic response
}

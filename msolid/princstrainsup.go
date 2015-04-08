// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/tsr"
)

// PrincStrainsUp implements stress-update in principal strains space
type PrincStrainsUp struct {

	// constants
	Pert  float64 // perturbation values
	EvTol float64 // tolerance to detect repeated eigenvalues
	Zero  float64 // minimum λ to be considered zero
	Fzero float64 // zero yield function value
	Nsig  int     // number of stress components

	// model
	mdl   EPmodel // elastoplastic model
	nalp  int     // number of α
	nsurf int     // number of yield functions
	fcoef float64 // coefficient to normalise yield function

	// variables
	Lσ   []float64   // eigenvalues of stresses
	Lε   []float64   // eigenvalues of strains
	P    [][]float64 // eigenprojectors of strains and stresses
	εetr []float64   // trial elastic state
	αn   []float64   // α at beginning of update

	// gradients
	Nb []float64 // principal values: plastic flow direction
	h  []float64 // principal values: hardening

	// for Jacobian
	N  []float64
	Mb [][]float64
	a  [][]float64
	b  [][]float64
	c  [][]float64

	// nonlinear solver
	x   []float64    // {εe0, εe1, εe2, α0, α1, ..., Δγ}
	nls num.NlSolver // nonlinear solver
}

// Init initialises this structure
func (o *PrincStrainsUp) Init(ndim int, prms fun.Prms, mdl EPmodel) (err error) {

	// constants
	o.Pert = 1e-7
	o.EvTol = tsr.EV_EVTOL
	o.Zero = tsr.EV_ZERO
	o.Fzero = 1e-9
	o.Nsig = 2 * ndim

	// model
	o.mdl = mdl
	o.nalp, o.nsurf, o.fcoef, _, _ = o.mdl.Info()

	// variables
	o.Lσ = make([]float64, 3)
	o.Lε = make([]float64, 3)
	o.P = la.MatAlloc(3, o.Nsig)
	o.εetr = make([]float64, o.Nsig)
	o.αn = make([]float64, o.nalp)
	o.Nb = make([]float64, 3)
	o.h = make([]float64, 3)
	o.x = make([]float64, 4+o.nalp)

	// nonlinear solver function
	ffcn := func(fx, x []float64) error {
		εe, α, Δγ := x[:3], x[3:3+o.nalp], x[3+o.nalp]
		εetr := o.Lε
		o.mdl.PVE_CalcSig(o.Lσ, εe)
		f, err := o.mdl.PVE_FlowHard(o.Nb, o.h, o.Lσ, α)
		if err != nil {
			return err
		}
		for i := 0; i < 3; i++ {
			fx[i] = εe[i] - εetr[i] + Δγ*o.Nb[i]
		}
		for i := 0; i < o.nalp; i++ {
			fx[3+i] = α[i] - o.αn[i] - Δγ*o.h[i]
		}
		fx[3+o.nalp] = f / o.fcoef
		return nil
	}

	// nonlinear solver Jacobian
	JfcnD := func(dfdx [][]float64, x []float64) error {
		εe, α, Δγ := x[:3], x[3:3+o.nalp], x[3+o.nalp]
		_ = α
		o.mdl.PVE_CalcSig(o.Lσ, εe)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				dfdx[i][j] = tsr.IIm[i][j] + Δγ*o.Mb[i][j]
			}
			for j := 0; j < o.nalp; j++ {
				dfdx[i][3+j] = Δγ * o.a[j][i]
				dfdx[3+j][i] = -Δγ * o.b[j][i]
			}
			dfdx[i][3+o.nalp] = o.Nb[i]
			dfdx[3+o.nalp][i] = o.N[i] / o.fcoef
		}
		// TODO
		//for i := 0; i < o.nalp; i++ {
		//dfdx[3+i]
		//}
		return nil
	}

	// nonlinear solver
	o.nls.Init(4+o.nalp, ffcn, nil, JfcnD, false, true, map[string]float64{})
	return
}

// Update updates state
func (o PrincStrainsUp) Update(s *State, ε, Δε []float64) (err error) {

	// trial strains and stresses
	o.mdl.ElastUpdate(s, ε, Δε)

	// check loading condition
	ftr := o.mdl.YieldFuncs(s)[0]
	io.Pforan("ftr = %v\n", ftr)
	if ftr <= o.Fzero {
		s.Loading = false
		return
	}

	// eigenvalues/projectors
	copy(o.εetr, s.Phi)
	//io.Pfred("εetr = %v\n", o.εetr)
	_, err = tsr.M_FixZeroOrRepeated(o.Lε, o.εetr, o.Pert, o.EvTol, o.Zero)
	if err != nil {
		return
	}
	err = tsr.M_EigenValsProjsNum(o.P, o.Lε, o.εetr)
	if err != nil {
		return
	}
	//io.Pforan("Δε = %v\n", Δε)
	//io.Pforan("Lε = %v\n", o.Lε)
	//io.Pfblue2("P[0] = %v\n", o.P[0])
	//io.Pfcyan("P[1] = %v\n", o.P[1])
	//io.Pfblue2("P[2] = %v\n", o.P[2])

	// trial values
	for i := 0; i < 3; i++ {
		o.x[i] = o.Lε[i]
	}
	for i := 0; i < o.nalp; i++ {
		o.αn[i] = s.Alp[i]
		o.x[3+i] = s.Alp[i]
	}
	o.x[3+o.nalp] = 0 // Δγ

	// solve
	silent := false
	err = o.nls.Solve(o.x, silent)
	if err != nil {
		return
	}

	// set new state
	εe, α, Δγ := o.x[:3], o.x[3:3+o.nalp], o.x[3+o.nalp]
	o.mdl.PVE_CalcSig(o.Lσ, εe)
	for i := 0; i < o.Nsig; i++ {
		s.Sig[i] = o.Lσ[0]*o.P[0][i] + o.Lσ[1]*o.P[1][i] + o.Lσ[2]*o.P[2][i]
	}
	copy(s.Alp, α)
	copy(s.Phi, εe)
	s.Dgam = Δγ
	s.Loading = true
	return
}

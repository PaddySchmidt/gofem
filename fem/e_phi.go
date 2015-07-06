// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// ElemPhi implementes a general element to solve the following equation
//     dφ       ∂φ
//     -- + v . -- = s(x)
//     dt       ∂x
// Notes: v is a constant vector
type ElemPhi struct {

	// basic data
	Cid int         // cell/element id
	X   [][]float64 // [ndim][nnode] matrix of nodal coordinates
	Shp *shp.Shape  // shape structure
	Nu  int         // total number of unknowns == number of vertices

	// integration points
	IpsElem []*shp.Ipoint // [nip] integration points of element

	// local starred variables
	ψs []float64 // [nip] ψ* = β1.φ + β2.dφdt

	// problem variables
	Umap []int // assembly map (location array/element equations)
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["phi"] = func(cellType string, faceConds []*FaceCond) *Info {

		// new info
		var info Info

		// return information
		return &info
	}

	// element allocator
	eallocators["phi"] = func(cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemPhi
		o.Cid = cid
		o.X = x
		o.Shp = shp.Get(cellType) // cellType: e.g. "tri6", "qua8"
		o.Nu = o.Shp.Nverts

		// integration points
		o.IpsElem, _ = GetIntegrationPoints(edat.Nip, edat.Nipf, cellType)
		if o.IpsElem == nil {
			return nil // => failed
		}

		// local starred variables
		nip := len(o.IpsElem)
		o.ψs = make([]float64, nip)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemPhi) Id() int {
	return o.Cid
}

// SetEqs set equations
func (o *ElemPhi) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Nu; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return true
}

// SetEleConds set element conditions
func (o *ElemPhi) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	return true
}

// InterpStarVars interpolate star variables to integration points
func (o *ElemPhi) InterpStarVars(sol *Solution) (ok bool) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		//interpolate starred variables
		o.ψs[idx] = 0
		for m := 0; m < o.Nu; m++ {
			o.ψs[idx] += o.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return true
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemPhi) AddToRhs(fb []float64, sol *Solution) (ok bool) {
	return true
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemPhi) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {
	return true
}

// Update perform (tangent) update
func (o *ElemPhi) Update(sol *Solution) (ok bool) {
	return true
}

// Encode encodes internal variables
func (o *ElemPhi) Encode(enc Encoder) (ok bool) {
	return true
}

// Decode decodes internal variables
func (o *ElemPhi) Decode(dec Decoder) (ok bool) {
	return true
}

// OutIpsData returns data from all integration points for output
func (o *ElemPhi) OutIpsData() (data []*OutIpData) {
	return
}

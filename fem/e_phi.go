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

	// scratchpad. computed @ each ip
	K [][]float64 // [nu][nu] consistent tangent matrix

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

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)

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
		if LogErr(o.Shp.CalcAtIp(o.X, ip, false), "InterpStarVars") {
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

	// auxiliary
	β1 := Global.DynCoefs.β1
	ndim := Global.Ndim
	nverts := o.Shp.Nverts

	v := []float64{1, 0, 0} // TODO: find a way to input the velocity

	// for each integration point
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// auxiliary variables
		coef := o.Shp.J * ip.W
		S := o.Shp.S
		G := o.Shp.G

		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			r := o.Umap[m] // row in the global vector
			for n := 0; n < nverts; n++ {
				fb[r] -= coef * S[m] * S[n] * β1 * sol.Y[o.Umap[n]]
				for j := 0; j < ndim; j++ {
					fb[r] += coef * v[j] * G[m][j] * S[n]
				}
			}
		}
	}
	return true
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemPhi) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// auxiliary
	β1 := Global.DynCoefs.β1
	nverts := o.Shp.Nverts

	// zero K matrix
	la.MatFill(o.K, 0)

	// for each integration point
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// auxiliary variables
		coef := o.Shp.J * ip.W
		S := o.Shp.S

		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				o.K[m][n] += coef * S[m] * S[n] * β1
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return true
}

// Update perform (tangent) update
func (o *ElemPhi) Update(sol *Solution) (ok bool) {
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

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

	// for each integration point
	ndim := Global.Ndim
	Gphi := make([]float64, ndim)
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}
		G := o.Shp.G

		// calculate function
		calc := func(sol *Solution) (vals map[string]float64) {
			for i := 0; i < ndim; i++ {
				Gphi[i] = 0
				for m := 0; m < o.Nu; m++ {
					Gphi[i] += G[m][i] * sol.Y[o.Umap[m]]
				}
			}
			vals = map[string]float64{
				"Gphix": Gphi[0],
				"Gphiy": Gphi[1],
			}
			if ndim == 3 {
				vals["Gphiz"] = Gphi[2]
			}
			return
		}

		// results
		x := o.Shp.IpRealCoords(o.X, ip)
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}

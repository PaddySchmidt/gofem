// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"

	// "fmt"
	"math"
)

// ElemPhi implementes a general element to solve the following equation
//     dφ       ∂φ
//     -- + v . -- = s(x)
//     dt       ∂x
// Notes: v is a constant vector
type ElemPhi struct {
	reinit bool

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // [ndim][nnode] matrix of nodal coordinates

	Nu   int // total number of unknowns == number of vertices
	Ndim int // space dimension

	// integration points
	IpsElem []shp.Ipoint // integration points of element
	IpsFace []shp.Ipoint // integration points corresponding to faces

	// local starred variables
	ψs      []float64 // [nip] ψ* = β1.φ + β2.dφdt
	PhiSign []float64
	// scratchpad. computed @ each ip
	K [][]float64 // [nu][nu] consistent tangent matrix
	M [][]float64

	// problem variables
	Umap []int // assembly map (location array/element equations)

	v_0 [][]float64 // TODO: Should be input parameter
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["phi"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		nverts := cell.Shp.Nverts
		ykeys := []string{"h"}

		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		info.T1vars = ykeys

		// return information
		return &info
	}

	// element allocator
	eallocators["phi"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemPhi
		o.Cell = cell
		o.X = x
		o.Nu = o.Cell.Shp.Nverts
		o.Ndim = sim.Ndim

		// integration points
		var err error
		o.IpsElem, o.IpsFace, err = o.Cell.Shp.GetIps(edat.Nip, edat.Nipf)
		if err != nil {
			chk.Panic("cannot allocate integration points of solid element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}

		// local starred variables
		nip := len(o.IpsElem)
		o.ψs = make([]float64, nip)
		o.PhiSign = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.v_0 = la.MatAlloc(nip, o.Ndim)
		o.reinit = false
		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemPhi) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *ElemPhi) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Nu; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return
}

// SetEleConds set element conditions
func (o *ElemPhi) SetEleConds(key string, f fun.Func, extra string) (err error) {
	return
}

// InterpStarVars interpolate star variables to integration points
func (o *ElemPhi) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		//interpolate starred variables
		o.ψs[idx] = 0
		for m := 0; m < o.Nu; m++ {
			o.ψs[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemPhi) AddToRhs(fb []float64, sol *Solution) (err error) {

	// auxiliary
	ndim := o.Ndim
	nverts := o.Cell.Shp.Nverts
	dt := sol.Dt * 2
	// for each integration point
	for iip, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// auxiliary variables
		coef := o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		o.Speed(S, G, ndim, iip, sol)
		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			r := o.Umap[m] // row in the global vector

			// additional term for reinitialisation process
			if o.reinit {
				fb[r] -= coef * S[m] * o.PhiSign[iip]
			}
			for n := 0; n < nverts; n++ {
				fb[r] += coef * (S[m]*S[n] + dt*dt/6.0*(o.v_0[iip][0]*G[m][0]+o.v_0[iip][1]*G[m][1])*(o.v_0[iip][0]*G[n][0]+o.v_0[iip][1]*G[n][1])) * (sol.Y[o.Umap[n]] - sol.Psi[o.Umap[n]]) / dt
				fb[r] += coef * (S[m] + o.v_0[iip][0]*dt/2.0*G[m][0] + o.v_0[iip][1]*dt/2.0*G[m][1]) * (o.v_0[iip][0]*G[n][0]*sol.Psi[o.Umap[n]] + o.v_0[iip][1]*G[n][1]*sol.Psi[o.Umap[n]])
			}
		}
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemPhi) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// auxiliary
	nverts := o.Cell.Shp.Nverts

	// zero K matrix
	la.MatFill(o.K, 0)
	dt := sol.Dt * 2

	// for each integration point
	for iip, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// auxiliary variables
		coef := o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				o.K[m][n] -= coef * (S[m]*S[n] + dt*dt/6.0*(o.v_0[iip][0]*G[m][0]+o.v_0[iip][1]*G[m][1])*(o.v_0[iip][0]*G[n][0]+o.v_0[iip][1]*G[n][1])) / dt
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// TODO: Remove Speed function and treat speed as input parameter
func (o *ElemPhi) Speed(S []float64, G [][]float64, ndim, iip int, sol *Solution) {
	// Speed orthogonal if reinit function
	if o.reinit {
		var GradPhi []float64
		var Phi float64

		for inode := 0; inode < 4; inode++ {
			// phi value at integration point
			Phi = 0.0
			// gradient of phi allocation
			if ndim == 2 {
				GradPhi = []float64{0., 0.}
			} else if ndim == 3 {
				GradPhi = []float64{0., 0., 0.}
			}

			// Get phi and gradient of phi at integration point level
			for m := 0; m < 4; m++ {
				Phi += S[m] * sol.Psi[o.Umap[m]]
				for j := 0; j < ndim; j++ {
					GradPhi[j] += G[m][j] * sol.Psi[o.Umap[m]]
				}
			}
		}

		// calculate length of gradient
		lenGrad := math.Sqrt(math.Pow(GradPhi[0], 2) + math.Pow(GradPhi[1], 2))

		// avoiding division by zero
		if math.Abs(lenGrad) <= 1.e-11 {
			lenGrad = 1.e-10
		}

		// smooth approximation of sign function
		o.PhiSign[iip] = Phi / math.Sqrt(Phi*Phi+lenGrad/10.0)

		// calculate orthogonal speed at integration point
		o.v_0[iip][0] = o.PhiSign[iip] * GradPhi[0] / lenGrad
		o.v_0[iip][1] = o.PhiSign[iip] * GradPhi[1] / lenGrad
	} else {
		// Arbitrary speed
		o.v_0[iip][0] = 1.0
		o.v_0[iip][1] = 0.0
		if ndim == 3 {
			o.v_0[iip][2] = 1.0
		}
	}
}

// Update perform (tangent) update
func (o *ElemPhi) Update(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *ElemPhi) Encode(enc Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *ElemPhi) Decode(dec Decoder) (err error) {
	return
}

// OutIpsData returns data from all integration points for output
func (o *ElemPhi) OutIpsData() (data []*OutIpData) {

	// for each integration point
	Gphi := make([]float64, o.Ndim)
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		err := o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		G := o.Cell.Shp.G

		// calculate function
		calc := func(sol *Solution) (vals map[string]float64) {
			for i := 0; i < o.Ndim; i++ {
				Gphi[i] = 0
				for m := 0; m < o.Nu; m++ {
					Gphi[i] += G[m][i] * sol.Y[o.Umap[m]]
				}
			}
			vals = map[string]float64{
				"Gphix": Gphi[0],
				"Gphiy": Gphi[1],
			}
			if o.Ndim == 3 {
				vals["Gphiz"] = Gphi[2]
			}
			return
		}

		// results
		x := o.Cell.Shp.IpRealCoords(o.X, ip)
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}

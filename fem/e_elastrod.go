// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// ElasticRod represents a structural rod element (for axial loads only) with 2 nodes only and
// simply implemented with constant stiffness matrix; i.e. no numerical integration is needed
type ElasticRod struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Nu   int         // total number of unknowns == 2 * nsn
	Ndim int         // space dimension

	// parameters and properties
	E float64 // Young's modulus
	A float64 // cross-sectional area
	L float64 // length of rod

	// variables for dynamics
	Rho  float64  // density of solids
	Gfcn fun.Func // gravity function

	// vectors and matrices
	T [][]float64 // [ndim][nu] transformation matrix: system aligned to rod => element system
	K [][]float64 // [nu][nu] element K matrix
	M [][]float64 // [nu][nu] element M matrix

	// problem variables
	Umap []int // assembly map (location array/element equations)

	// scratchpad. computed @ each ip
	grav []float64 // [ndim] gravity vector
	us   []float64 // [ndim] displacements @ ip
	fi   []float64 // [nu] internal forces
	ue   []float64 // local u vector
}

// register element
func init() {

	// information allocator
	infogetters["elastrod"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		// solution variables
		ykeys := []string{"ux", "uy"}
		if sim.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz"}
		}
		info.Dofs = make([][]string, 2)
		for m := 0; m < 2; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz"}

		// t1 and t2 variables
		info.T2vars = ykeys
		return &info
	}

	// element allocator
	eallocators["elastrod"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// check
		ndim := len(x)
		if ndim == 3 {
			chk.Panic("elastrod is not implemented for 3D yet")
		}

		// basic data
		var o ElasticRod
		o.Cell = cell
		o.X = x
		o.Ndim = sim.Ndim
		o.Nu = o.Ndim * 2

		// parameters
		matdata := sim.MatParams.Get(edat.Mat)
		if matdata == nil {
			chk.Panic("cannot get materials data for elastic rod element {tag=%d id=%d material=%q}", cell.Tag, cell.Id, edat.Mat)
		}

		// parameters
		for _, p := range matdata.Prms {
			switch p.N {
			case "E":
				o.E = p.V
			case "A":
				o.A = p.V
			case "rho":
				o.Rho = p.V
			}
		}

		// vectors and matrices
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.M = la.MatAlloc(o.Nu, o.Nu)
		o.ue = make([]float64, o.Nu)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, o.Ndim)
		o.us = make([]float64, o.Ndim)
		o.fi = make([]float64, o.Nu)

		// geometry
		x0 := o.X[0][0]
		y0 := o.X[1][0]
		x1 := o.X[0][1]
		y1 := o.X[1][1]
		dx := x1 - x0
		dy := y1 - y0
		o.L = math.Sqrt(dx*dx + dy*dy)

		// global-to-local transformation matrix
		c := dx / o.L
		s := dy / o.L
		o.T = [][]float64{
			{c, s, 0, 0},
			{0, 0, c, s},
		}

		// K and M matrices
		α := o.E * o.A / o.L
		β := o.Rho * o.A * o.L / 6.0
		o.K = [][]float64{
			{+α * c * c, +α * c * s, -α * c * c, -α * c * s},
			{+α * c * s, +α * s * s, -α * c * s, -α * s * s},
			{-α * c * c, -α * c * s, +α * c * c, +α * c * s},
			{-α * c * s, -α * s * s, +α * c * s, +α * s * s},
		}
		o.M = [][]float64{
			{2.0 * β, 0.0, 1.0 * β, 0.0},
			{0.0, 2.0 * β, 0.0, 1.0 * β},
			{1.0 * β, 0.0, 2.0 * β, 0.0},
			{0.0, 1.0 * β, 0.0, 2.0 * β},
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElasticRod) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *ElasticRod) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < 2; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Umap[r] = eqs[m][i]
		}
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *ElasticRod) InterpStarVars(sol *Solution) (err error) {
	// TODO: dynamics
	chk.Panic("ElasticRod cannot handle dynamics yet")
	return
}

// SetEleConds set element conditions
func (o *ElasticRod) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" {
		chk.Panic("ElasticRod cannot handle gravity yet")
		o.Gfcn = f
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o ElasticRod) AddToRhs(fb []float64, sol *Solution) (err error) {
	/*
		for m := 0; m < 2; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.Umap[i+m*o.Ndim]
				for n := 0; n < 2; n++ {
					for j := 0; j < o.Ndim; j++ {
						c := o.Umap[j+n*o.Ndim]
						fb[r] -= o.K[r][c] * sol.Y[c]
					}
				}
			}
		}
	*/
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			fb[I] -= o.K[i][j] * sol.Y[J] // -fi
		}
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o ElasticRod) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// Update perform (tangent) update
func (o *ElasticRod) Update(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElasticRod) Encode(enc Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o ElasticRod) Decode(dec Decoder) (err error) {
	return
}

// OutIpsData returns data from all integration points for output
func (o ElasticRod) OutIpsData() (data []*OutIpData) {
	x := make([]float64, o.Ndim)
	for i := 0; i < o.Ndim; i++ {
		x[i] = (o.X[i][0] + o.X[i][1]) / 2.0 // centroid
	}
	calc := func(sol *Solution) (vals map[string]float64) {
		ua := make([]float64, 2) // axial displacements
		for i := 0; i < 2; i++ {
			for j, J := range o.Umap {
				ua[i] += o.T[i][j] * sol.Y[J]
			}
		}
		ea := (ua[1] - ua[0]) / o.L // axial strain
		vals = make(map[string]float64)
		vals["sig"] = o.E * ea // axial stress
		return
	}
	data = append(data, &OutIpData{o.Id(), x, calc})
	return
}

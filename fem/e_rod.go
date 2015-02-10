// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Rod represents a structural rod element (for only axial loads)
type Rod struct {

	// basic data
	Cell *inp.Cell   // cell
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Nu   int         // total number of unknowns == 2 * nsn

	// parameters
	E float64 // Young's modulus
	A float64 // cross-sectional area

	// variables for dynamics
	Rho  float64  // density of solids
	Gfcn fun.Func // gravity function

	// integration points
	IpsElem []*shp.Ipoint // integration points of element

	// vectors and matrices
	K   [][]float64 // global K matrix
	M   [][]float64 // global M matrices
	Rus []float64   // residual: Rus = fi - fx

	// problem variables
	Umap []int // assembly map (location array/element equations)

	// internal variables
	Sig []float64
	//States    []*msolid.State
	//StatesBkp []*msolid.State

	// scratchpad. computed @ each ip
	grav []float64 // [ndim] gravity vector
	us   []float64 // [ndim] displacements @ ip
	fi   []float64 // [nu] internal forces
	ue   []float64 // local u vector
}

// register element
func init() {

	// information allocator
	iallocators["rod"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) (info Info) {

		// solution variables
		ykeys := []string{"ux", "uy"}
		if msh.Ndim == 3 {
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
		return
	}

	// element allocator
	eallocators["rod"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o Rod
		o.Cell = msh.Cells[cid]
		o.X = BuildCoordsMatrix(o.Cell, msh)
		o.Ndim = msh.Ndim
		o.Nu = o.Ndim * o.Cell.Shp.Nverts

		// parameters
		mat := global.Mdb.GetOrPanic(edat.Mat)
		for _, p := range mat.Prms {
			switch p.N {
			case "E":
				o.E = p.V
			case "A":
				o.A = p.V
			}
		}

		// integration points
		var nip int
		if s_nip, found := utl.Keycode(edat.Extra, "nip"); found {
			nip = utl.Atoi(s_nip)
		}
		o.IpsElem = shp.GetIps(o.Cell.Shp.Type, nip)
		nip = len(o.IpsElem)

		// state
		o.Sig = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.M = la.MatAlloc(o.Nu, o.Nu)
		o.ue = make([]float64, o.Nu)
		o.Rus = make([]float64, o.Nu)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, o.Ndim)
		o.us = make([]float64, o.Ndim)
		o.fi = make([]float64, o.Nu)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *Rod) SetEqs(eqs [][]int, mixedform_eqs []int) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Umap[r] = eqs[m][i]
		}
	}
}

// InterpStarVars interpolates star variables to integration points
func (o *Rod) InterpStarVars(sol *Solution) (err error) {

	// skip steady cases
	if global.Sim.Data.Steady {
		return
	}

	/*
		// for each integration point
		for idx, ip := range o.IpsElem {

			// interpolation functions and gradients
			err = o.log(o.Cell.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars")
			if err != nil {
				return
			}

			// interpolate starred variables
			o.divχs[idx] = 0
			for i := 0; i < o.Ndim; i++ {
				o.ζs[idx][i] = 0
				o.χs[idx][i] = 0
				for m := 0; m < o.Cell.Shp.Nverts; m++ {
					r := o.Umap[i+m*o.Ndim]
					o.ζs[idx][i] += o.Cell.Shp.S[m] * sol.Zet[r]
					o.χs[idx][i] += o.Cell.Shp.S[m] * sol.Chi[r]
					o.divχs[idx] += o.Cell.Shp.G[m][i] * sol.Chi[r]
				}
			}
		}
	*/

	// success
	return
}

// SetEleConds set element conditions
func (o *Rod) SetEleConds(key string, f fun.Func, extra string) {

	// gravity
	if key == "g" {
		o.Gfcn = f
		return
	}
}

// adds -R to global residual vector fb
func (o Rod) AddToRhs(fb []float64, sol *Solution) (err error) {

	// for each integration point
	dc := global.DynCoefs
	_ = dc
	nverts := o.Cell.Shp.Nverts
	ndim := o.Ndim
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.log(o.ipvars(idx, sol), "AddToRhs")
		if err != nil {
			return
		}
		coef := ip.W
		G := o.Cell.Shp.Gvec
		Jvec := o.Cell.Shp.Jvec3d

		//σ := o.States[ix].Sig[0]
		for m := 0; m < nverts; m++ {
			for i := 0; i < ndim; i++ {
				r := o.Umap[i+m*ndim]
				fb[r] -= coef * o.A * o.Sig[idx] * G[m] * Jvec[i] // +fi
			}
			//if o.hasg {
			//o.Rus[o.nd-1 + m*o.nd] -= o.coef * gcmp * o.gu.S[m] // -fx
			//}
		}
	}

	return

}

// adds element K to global Jacobian matrix Kb
func (o Rod) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// zero K matrix
	la.MatFill(o.K, 0)
	la.MatFill(o.M, 0)

	// for each integration point
	dc := global.DynCoefs
	_ = dc
	nverts := o.Cell.Shp.Nverts
	ndim := o.Ndim
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.log(o.ipvars(idx, sol), "AddToKb")
		if err != nil {
			return
		}

		coef := ip.W
		G := o.Cell.Shp.Gvec
		J := o.Cell.Shp.J
		Jvec := o.Cell.Shp.Jvec3d

		// add contribution to consistent tangent matrix
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				for i := 0; i < ndim; i++ {
					for j := 0; j < ndim; j++ {
						r := i + m*ndim
						c := j + n*ndim
						o.K[r][c] += coef * o.A * o.E * G[m] * G[n] * Jvec[i] * Jvec[j] / J
						//if !steady {
						//o.M[r][c] += o.ipe[idx][3] * o.ρ * o.A * o.gu.S[m] * o.gu.S[n] * o.gu.Jvec[i] * o.gu.Jvec[j] / o.gu.J
						//}
					}
				}
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}

	// success
	return
}

// Update perform (tangent) update
func (o *Rod) Update(sol *Solution) (err error) {

	// for each integration point
	nverts := o.Cell.Shp.Nverts
	ndim := o.Ndim
	for idx, _ := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.log(o.ipvars(idx, sol), "Update")
		if err != nil {
			return
		}

		G := o.Cell.Shp.Gvec
		J := o.Cell.Shp.J
		Jvec := o.Cell.Shp.Jvec3d

		// compute strains
		Δε := 0.0
		for m := 0; m < nverts; m++ {
			for i := 0; i < ndim; i++ {
				r := o.Umap[i+m*ndim]
				Δε += G[m] * Jvec[i] * sol.Y[r] / J
			}
		}

		// call model update => update stresses
		o.Sig[idx] = Δε * o.E
	}

	// success
	return
}

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *Rod) SetSurfLoads(key string, idxface int, f fun.Func, extra string) {
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// log logs errors
func (o *Rod) log(err error, msg string) error {
	if err != nil {
		log.Printf("Rod: eid=%d %s failed with %v\n", o.Cell.Id, msg, err)
	}
	return err
}

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Rod) ipvars(idx int, sol *Solution) (err error) {

	// interpolation functions and gradients
	err = o.log(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars")
	if err != nil {
		return
	}

	// skip if steady (this must be after CalcAtIp, because callers will need S and G)
	if global.Sim.Data.Steady {
		return
	}

	// clear variables
	for i := 0; i < o.Ndim; i++ {
		o.us[i] = 0
	}

	// recover u-variables @ ip
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := o.Umap[i+m*o.Ndim]
			o.us[i] += o.Cell.Shp.S[m] * sol.Y[r]
		}
	}

	// success
	return
}
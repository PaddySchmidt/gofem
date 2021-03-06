// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"
	"math"
	"sort"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// EssentialBc holds information about essential bounday conditions such as constrained nodes.
// Lagrange multipliers are used to implement both single- and multi-point constraints.
//  In general, essential bcs / constraints are defined by means of:
//
//      A * y = c
//
//  The resulting Kb matrix will then have the following form:
//      _       _
//     |  K  At  | / δy \   / -R - At*λ \
//     |         | |    | = |           |
//     |_ A   0 _| \ δλ /   \  c - A*y  /
//         Kb       δyb          fb
//
type EssentialBc struct {
	Key   string    // ux, uy, rigid, incsup
	Eqs   []int     // equations
	ValsA []float64 // values for matrix A
	Fcn   fun.Func  // function that implements the "c" in A * y = c
	Inact bool      // inactive
}

// EssentialBcs implements a structure to record the definition of essential bcs / constraints.
// Each constraint will have a unique Lagrange multiplier index.
type EssentialBcs struct {
	HydFcn *HydroStatic   // for computing hydrostatic conditions
	Eq2idx map[int][]int  // maps eq number to indices in BcsTmp
	Bcs    []*EssentialBc // active essential bcs / constraints
	A      la.Triplet     // matrix of coefficients 'A'
	Am     *la.CCMatrix   // compressed form of A matrix

	// temporary
	BcsTmp eqbcpairs // temporary essential bcs / constraints, including inactive ones. maps the first equation number to bcs
}

// Reset initialises this structure. It also performs a reset of internal structures.
func (o *EssentialBcs) Init(hydfcn *HydroStatic) {
	o.BcsTmp = make([]eqbcpair, 0)
	o.Eq2idx = make(map[int][]int)
	o.Bcs = make([]*EssentialBc, 0)
	o.HydFcn = hydfcn
}

// Build builds this structure and its iternal data
//  nλ -- is the number of essential bcs / constraints == number of Lagrange multipliers
//  nnzA -- is the number of non-zeros in matrix 'A'
func (o *EssentialBcs) Build(ny int) (nλ, nnzA int) {

	// sort bcs to make sure all processors will number Lagrange multipliers in the same order
	sort.Sort(o.BcsTmp)

	// count number of active constraints and non-zeros in matrix A
	for _, pair := range o.BcsTmp {
		if !pair.bc.Inact {
			o.Bcs = append(o.Bcs, pair.bc)
			nλ += 1
			nnzA += len(pair.bc.ValsA)
		}
	}

	// skip if there are no constraints
	if nλ == 0 {
		return
	}

	// set matrix A
	o.A.Init(nλ, ny, nnzA)
	for i, c := range o.Bcs {
		for j, eq := range c.Eqs {
			o.A.Put(i, eq, c.ValsA[j])
		}
	}
	o.Am = o.A.ToMatrix(nil)

	// debug
	if false {
		log.Printf("\n\nAm=%v\n", o.Am)
	}
	return
}

// AddtoRhs adds the essential bcs / constraints terms to the augmented fb vector
func (o EssentialBcs) AddToRhs(fb []float64, sol *Solution) {

	// skip if there are no constraints
	if len(o.Bcs) == 0 {
		return
	}

	// add -At*λ to fb
	la.SpMatTrVecMulAdd(fb, -1, o.Am, sol.L) // fb += -1 * At * λ

	// assemble -rc = c - A*y into fb
	ny := len(sol.Y)
	for i, c := range o.Bcs {
		fb[ny+i] = c.Fcn.F(sol.T, nil)
	}
	la.SpMatVecMulAdd(fb[ny:], -1, o.Am, sol.Y) // fb += -1 * A * y
}

// add adds new essential bcs / constraint and sets map eq2idx
func (o *EssentialBcs) add(key string, eqs []int, valsA []float64, fcn fun.Func) {
	idx := len(o.BcsTmp)
	o.BcsTmp = append(o.BcsTmp, eqbcpair{eqs[0], &EssentialBc{key, eqs, valsA, fcn, false}})
	for _, eq := range eqs {
		utl.IntIntsMapAppend(&o.Eq2idx, eq, idx)
	}
}

// add_single adds single-point constraint
func (o *EssentialBcs) add_single(key string, eq int, fcn fun.Func) {
	for _, idx := range o.Eq2idx[eq] {
		pair := o.BcsTmp[idx]
		if pair.bc.Key == "rigid" || pair.bc.Key == "incsup" {
			return
		}
		pair.bc.Inact = true
	}
	o.add(key, []int{eq}, []float64{1}, fcn)
}

// GetFirstYandCmap returns the initial "yandc" map with additional keys that EssentialBcs can handle
//  rigid  -- define rigid element constraints
//  incsup -- inclined support constraints
//  hst    -- set hydrostatic pressures
func GetIsEssenKeyMap() map[string]bool {
	return map[string]bool{"rigid": true, "incsup": true, "hst": true}
}

// Set sets a constraint if it does NOT exist yet.
//  key   -- can be Dof key such as "ux", "uy" or constraint type such as "mpc" or "rigid"
//  extra -- is a keycode-style data. e.g. "!type:incsup2d !alp:30"
//  Notes: 1) the default for key is single point constraint; e.g. "ux", "uy", ...
//         2) hydraulic head can be set with key == "H"
func (o *EssentialBcs) Set(key string, nodes []*Node, fcn fun.Func, extra string) (err error) {

	// len(nod) must be greater than 0
	chk.IntAssertLessThan(0, len(nodes)) // 0 < len(nod)

	// skip nil node
	if nodes[0] == nil {
		return
	}

	// space dimension
	ndim := len(nodes[0].Vert.C)

	// rigid element
	if key == "rigid" {
		a := nodes[0].Dofs
		for i := 1; i < len(nodes); i++ {
			for j, b := range nodes[i].Dofs {
				o.add(key, []int{a[j].Eq, b.Eq}, []float64{1, -1}, &fun.Zero)
			}
		}
		return // success
	}

	// inclined support
	if key == "incsup" {

		// check
		if ndim != 2 {
			return chk.Err("inclined support works only in 2D for now")
		}

		// get data
		var α float64
		if val, found := io.Keycode(extra, "alp"); found {
			α = io.Atof(val) * math.Pi / 180.0
		}
		co, si := math.Cos(α), math.Sin(α)

		// set for all nodes
		for _, nod := range nodes {

			// find existent constraints and deactivate them
			eqx := nod.Dofs[0].Eq
			eqy := nod.Dofs[1].Eq
			for _, eq := range []int{eqx, eqy} {
				for _, idx := range o.Eq2idx[eq] {
					pair := o.BcsTmp[idx]
					if pair.bc.Key != "rigid" {
						pair.bc.Inact = true
					}
				}
			}

			// set constraint
			o.add(key, []int{eqx, eqy}, []float64{co, si}, &fun.Zero)
		}
		return // success
	}

	// hydraulic head
	if key == "hst" {

		// set for all nodes
		for _, nod := range nodes {

			// create function
			// Note: fcn is a shift such that  pl = pl(z) - shift(t)
			d := nod.GetDof("pl")
			if d == nil {
				continue // node doesn't have key. ex: pl in qua8/qua4 elements
			}
			z := nod.Vert.C[1] // 2D
			if ndim == 3 {
				z = nod.Vert.C[2] // 3D
			}
			plVal, _, err := o.HydFcn.Calc(z)
			if err != nil {
				return chk.Err("cannot set hst (hydrostatic) essential boundary condition")
			}
			pl := fun.Add{
				B: 1, Fb: &fun.Cte{C: plVal},
				A: -1, Fa: fcn,
			}

			// set constraint
			o.add_single("pl", d.Eq, &pl)
		}
		return // success
	}

	// single-point constraint
	for _, nod := range nodes {
		d := nod.GetDof(key)
		if d == nil {
			return // success
		}
		o.add_single(key, d.Eq, fcn)
	}

	// success
	return
}

// auxiliary /////////////////////////////////////////////////////////////////////////////////////////

type eqbcpair struct {
	eq int
	bc *EssentialBc
}

type eqbcpairs []eqbcpair

func (o eqbcpairs) Len() int           { return len(o) }
func (o eqbcpairs) Swap(i, j int)      { o[i], o[j] = o[j], o[i] }
func (o eqbcpairs) Less(i, j int) bool { return o[i].eq < o[j].eq }

// List returns a simple list logging bcs at time t
func (o *EssentialBcs) List(t float64) (l string) {
	var pairs eqbcpairs
	for _, bc := range o.Bcs {
		for _, eq := range bc.Eqs {
			pairs = append(pairs, eqbcpair{eq, bc})
		}
	}
	sort.Sort(pairs)
	l = "\n  ====================================================================================\n"
	l += io.Sf("  %8s%8s%23s%23s\n", "eq", "key", "value @ t=0", io.Sf("value @ t=%g", t))
	l += "  ------------------------------------------------------------------------------------\n"
	for _, p := range pairs {
		l += io.Sf("  %8d%8s%23.13f%23.13f\n", p.eq, p.bc.Key, p.bc.Fcn.F(0, nil), p.bc.Fcn.F(t, nil))
	}
	l += "  ====================================================================================\n"
	return
}

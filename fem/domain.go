// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Solution holds the solution data @ nodes.
//        / u \         / u \
//        |   | => y =  |   |
//  yb =  | p |         \ p / (ny x 1)
//        |   |
//        \ λ / (nyb x 1)
//
type Solution struct {

	// state
	T      float64   // current time
	Y      []float64 // DOFs (solution variables); e.g. y = {u, p}
	Dydt   []float64 // dy/dt
	D2ydt2 []float64 // d²y/dt²

	// auxiliary
	ΔY  []float64 // total increment (for nonlinear solver)
	Psi []float64 // t1 star vars; e.g. ψ* = β1.p + β2.dpdt
	Zet []float64 // t2 star vars; e.g. ζ* = α1.u + α2.v + α3.a
	Chi []float64 // t2 star vars; e.g. χ* = α4.u + α5.v + α6.a
	L   []float64 // Lagrange multipliers
}

// Domain holds all Nodes and Elements active during a stage in addition to the Solution at nodes.
// Only elements in this processor are recorded here; however information from
// all cells might be recorded as well.
type Domain struct {

	// init: region, mesh, linear solver
	Reg    *inp.Region // region data
	Msh    *inp.Mesh   // mesh data
	LinSol la.LinSol   // linear solver

	// stage: nodes (active) and elements (active AND in this processor)
	Nodes []*Node // active nodes (for each stage)
	Elems []Elem  // only active elements in this processor (for each stage)

	// stage: auxiliary maps for dofs and equation types
	F2Y      map[string]string // converts f-keys to y-keys; e.g.: "ux" => "fx"
	YandC    map[string]bool   // y and constraints keys; e.g. "ux", "pl", "H", "incsup", "rigid"
	Dof2Tnum map[string]int    // {t1,t2}-types: dof => t_number; e.g. "ux" => 2, "pl" => 1

	// stage: auxiliary maps for nodes and elements
	Vid2node   []*Node // [nverts] VertexId => index in Nodes. Inactive vertices are 'nil'
	Cid2elem   []Elem  // [ncells] CellId => index in Elems. Cells in other processors or inactive are 'nil'
	Cid2active []bool  // [ncells] CellId => whether cell is active or not in ANY processor

	// stage: subsets of elements
	ElemIntvars []ElemIntvars   // elements with internal vars in this processor
	ElemConnect []ElemConnector // connector elements in this processor
	ElemWriters []ElemWriter    // writer elements in this processor

	// stage: coefficients and prescribed forces
	EssenBcs EssentialBcs // constraints (Lagrange multipliers)
	PtNatBcs PtNaturalBcs // point loads such as prescribed forces at nodes

	// stage: t1 and t2 variables
	T1eqs []int // first t-derivative variables; e.g.:  dp/dt vars (subset of ykeys)
	T2eqs []int // second t-derivative variables; e.g.: d²u/dt² vars (subset of ykeys)

	// stage: dimensions
	NnzKb int // number of nonzeros in Kb matrix
	Ny    int // total number of dofs, except λ
	Nlam  int // total number of Lagrange multipliers
	NnzA  int // number of nonzeros in A (constraints) matrix
	Nyb   int // total number of equations: ny + nλ

	// stage: solution and linear solver
	Sol      *Solution   // solution state
	Kb       *la.Triplet // Jacobian == dRdy
	Fb       []float64   // residual == -fb
	Wb       []float64   // workspace
	InitLSol bool        // flag telling that linear solver needs to be initialised prior to any further call
}

// NewDomain returns a new domain
func NewDomain(reg *inp.Region) (o *Domain) {
	o = new(Domain)
	o.Reg = reg
	o.Msh = inp.ReadMsh(reg.Mshfile, global.Root)
	if global.Distr {
		PanicOrNot(global.Nproc != len(o.Msh.Part2cells), "number of processors must be equal to number of partitions defined in mesh file. %d != %d", global.Nproc, len(o.Msh.Part2cells))
	}
	o.LinSol = la.GetSolver(global.Sim.LinSol.Name)
	return
}

// SetStage set nodes, equation numbers and auxiliary data for given stage
func (o *Domain) SetStage(idxstg int, stg *inp.Stage) {

	// backup state
	if idxstg > 0 {
		o.create_stage_copy()
		o.fix_inact_flags(stg.Activate, false)
		o.fix_inact_flags(stg.Deactivate, true)
	}

	// nodes (active) and elements (active AND in this processor)
	o.Nodes = make([]*Node, 0)
	o.Elems = make([]Elem, 0)

	// auxiliary maps for dofs and equation types
	o.F2Y = make(map[string]string)
	o.YandC = GetIsEssenKeyMap()
	o.Dof2Tnum = make(map[string]int)

	// auxiliary maps for nodes and elements
	o.Vid2node = make([]*Node, len(o.Msh.Verts))
	o.Cid2elem = make([]Elem, len(o.Msh.Cells))
	o.Cid2active = make([]bool, len(o.Msh.Cells))

	// subsets of elements
	o.ElemConnect = make([]ElemConnector, 0)
	o.ElemWriters = make([]ElemWriter, 0)
	o.ElemIntvars = make([]ElemIntvars, 0)

	// allocate nodes and cells (active only) -------------------------------------------------------

	// for each cell
	var eq int // current equation number => total number of equations @ end of loop
	o.NnzKb = 0
	for _, c := range o.Msh.Cells {

		// get element data and information structure
		edat := o.Reg.Etag2data(c.Tag)
		if edat.Inact {
			continue
		}
		o.Cid2active[c.Id] = true
		info := GetElemInfo(edat, c.Id, o.Msh)
		utl.IntAssert(len(info.Dofs), len(c.Verts))

		// store y and f information
		for ykey, fkey := range info.Y2F {
			o.F2Y[fkey] = ykey
			o.YandC[ykey] = true
		}

		// t1 and t2 equations
		for _, ykey := range info.T1vars {
			o.Dof2Tnum[ykey] = 1
		}
		for _, ykey := range info.T2vars {
			o.Dof2Tnum[ykey] = 2
		}

		// loop over nodes of this element
		var eNdof int // number of DOFs of this elmeent
		for j, v := range c.Verts {

			// new or existent node
			var nod *Node
			if o.Vid2node[v] == nil {
				nod = NewNode(o.Msh.Verts[v])
				o.Vid2node[v] = nod
				o.Nodes = append(o.Nodes, nod)
			} else {
				nod = o.Vid2node[v]
			}

			// set DOFs and equation numbers
			for _, ukey := range info.Dofs[j] {
				eq = nod.AddDofAndEq(ukey, eq)
				eNdof += 1
			}
		}

		// allocate element
		mycell := c.Part == global.Rank // cell belongs to this processor
		if !global.Distr {
			mycell = true // not distributed simulation => this processor has all cells
		}
		if mycell {

			// new element
			ele := NewElem(edat, c.Id, o.Msh)
			o.Cid2elem[c.Id] = ele
			o.Elems = append(o.Elems, ele)

			// give equation numbers to new element
			eqs := make([][]int, len(c.Verts))
			for j, v := range c.Verts {
				for _, dof := range o.Vid2node[v].dofs {
					eqs[j] = append(eqs[j], dof.Eq)
				}
			}
			ele.SetEqs(eqs, nil)

			// subsets of elements
			o.add_element_to_subsets(ele)
		}

		// number of non-zeros
		o.NnzKb += eNdof * eNdof
	}

	// logging
	if global.Root {
		log.Printf("dom: stage # %d %s\n", idxstg, stg.Desc)
		log.Printf("dom: nnodes=%d nelems=%d\n", len(o.Nodes), len(o.Elems))
	}

	// element conditions, essential and natural boundary conditions --------------------------------

	// (re)set constraints and prescribed forces structures
	o.EssenBcs.Reset(o.Msh.Ndim)
	o.PtNatBcs.Reset()

	// element conditions
	for _, ec := range stg.EleConds {
		cells, ok := o.Msh.CellTag2cells[ec.Tag]
		PanicOrNot(!ok, "cannot find cells with tag = %d to assign conditions", ec.Tag)
		for _, c := range cells {
			e := o.Cid2elem[c.Id]
			if e != nil { // set conditions only for this processor's / active element
				for j, key := range ec.Keys {
					e.SetEleConds(key, global.Sim.Functions.GetOrPanic(ec.Funcs[j]), ec.Extra)
				}
			}
		}
	}

	// face boundary conditions
	for _, fc := range stg.FaceBcs {
		pairs, ok := o.Msh.FaceTag2cells[fc.Tag]
		PanicOrNot(!ok, "cannot find cells with face tag = %d to assign face boundary conditions", fc.Tag)
		for _, p := range pairs {
			if !o.Cid2active[p.C.Id] { // skip inactive element
				continue
			}
			localverts := p.C.Shp.FaceLocalV[p.Fid]
			var enodes []*Node
			for _, l := range localverts {
				v := p.C.Verts[l]
				enodes = append(enodes, o.Vid2node[v])
			}
			for j, key := range fc.Keys {
				if o.YandC[key] {
					o.EssenBcs.Set(key, enodes, global.Sim.Functions.GetOrPanic(fc.Funcs[j]), fc.Extra)
				} else {
					e := o.Cid2elem[p.C.Id]
					if e != nil { // set natural BCs only for this processor's / active element
						e.SetSurfLoads(key, p.Fid, global.Sim.Functions.GetOrPanic(fc.Funcs[j]), fc.Extra)
					}
				}
			}
		}
	}

	// vertex bounday conditions
	for _, nc := range stg.NodeBcs {
		verts, ok := o.Msh.VertTag2verts[nc.Tag]
		PanicOrNot(!ok, "cannot find vertices with tag = %d to assign node boundary conditions", nc.Tag)
		for _, v := range verts {
			if o.Vid2node[v.Id] != nil { // set BCs only for active nodes
				n := o.Vid2node[v.Id]
				for j, key := range nc.Keys {
					if o.YandC[key] {
						o.EssenBcs.Set(key, []*Node{n}, global.Sim.Functions.GetOrPanic(nc.Funcs[j]), nc.Extra)
					} else {
						o.PtNatBcs.Set(o.F2Y[key], n, global.Sim.Functions.GetOrPanic(nc.Funcs[j]), nc.Extra)
					}
				}
			}
		}
	}

	// resize slices --------------------------------------------------------------------------------

	// t1 and t2 equations
	o.T1eqs = make([]int, 0)
	o.T2eqs = make([]int, 0)
	for _, nod := range o.Nodes {
		for _, dof := range nod.dofs {
			switch o.Dof2Tnum[dof.Ukey] {
			case 1:
				o.T1eqs = append(o.T1eqs, dof.Eq)
			case 2:
				o.T2eqs = append(o.T2eqs, dof.Eq)
			default:
				PanicOrNot(true, "t1 and t2 equations are incorrectly set")
			}
		}
	}

	// size of arrays
	o.Ny = eq
	o.Nlam, o.NnzA = o.EssenBcs.Build(o.Ny)
	o.Nyb = o.Ny + o.Nlam

	// solution structure and linear solver
	o.Sol = new(Solution)
	o.Kb = new(la.Triplet)
	o.Fb = make([]float64, o.Nyb)
	o.Wb = make([]float64, o.Nyb)
	o.Kb.Init(o.Nyb, o.Nyb, o.NnzKb+2*o.NnzA)
	o.InitLSol = true // tell solver that lis has to be initialised before use

	// allocate arrays
	o.Sol.Y = make([]float64, o.Ny)
	o.Sol.ΔY = make([]float64, o.Ny)
	o.Sol.L = make([]float64, o.Nlam)
	if !global.Sim.Data.Steady {
		o.Sol.Dydt = make([]float64, o.Ny)
		o.Sol.D2ydt2 = make([]float64, o.Ny)
		o.Sol.Psi = make([]float64, o.Ny)
		o.Sol.Zet = make([]float64, o.Ny)
		o.Sol.Chi = make([]float64, o.Ny)
	}

	// initialise internal variables
	for _, e := range o.ElemIntvars {
		e.InitIvs(o.Sol)
	}

	// connect elements
	for _, e := range o.ElemConnect {
		e.Connect(o.Elems, o.Cid2elem)
	}

	// logging
	if global.Root {
		log.Printf("dom: essenbcs=%v\n", o.EssenBcs.List(stg.Control.Tf))
		log.Printf("dom: ptnatbcs=%v\n", o.PtNatBcs.List(stg.Control.Tf))
		log.Printf("dom: ny=%d nlam=%d nnzKb=%d nnzA=%d nt1eqs=%d nt2eqs=%d\n", o.Ny, o.Nlam, o.NnzKb, o.NnzA, len(o.T1eqs), len(o.T2eqs))
	}
}

// auxiliary functions //////////////////////////////////////////////////////////////////////////////

// add_element_to_subsets adds an Elem to many subsets as it fits
func (o *Domain) add_element_to_subsets(ele Elem) {
	if e, ok := ele.(ElemIntvars); ok {
		o.ElemIntvars = append(o.ElemIntvars, e)
	}
	if e, ok := ele.(ElemConnector); ok {
		o.ElemConnect = append(o.ElemConnect, e)
	}
	if e, ok := ele.(ElemWriter); ok {
		o.ElemWriters = append(o.ElemWriters, e)
	}
}

// star_vars computes starred variables
func (o *Domain) star_vars(Δt float64) (err error) {

	// skip if steady simulation
	if global.Sim.Data.Steady {
		return
	}

	// recompute coefficients
	dc := global.DynCoefs
	err = dc.CalcBoth(Δt)
	if err != nil {
		return
	}

	// compute starred vectors
	for _, I := range o.T1eqs {
		o.Sol.Psi[I] = dc.β1*o.Sol.Y[I] + dc.β2*o.Sol.Dydt[I]
	}
	for _, I := range o.T2eqs {
		o.Sol.Zet[I] = dc.α1*o.Sol.Y[I] + dc.α2*o.Sol.Dydt[I] + dc.α3*o.Sol.D2ydt2[I]
		o.Sol.Chi[I] = dc.α4*o.Sol.Y[I] + dc.α5*o.Sol.Dydt[I] + dc.α6*o.Sol.D2ydt2[I]
	}

	// set internal starred variables
	for _, e := range o.Elems {
		e.InterpStarVars(o.Sol)
	}
	return
}

// create_stage_copy creates a copy of current stage => to be used later when activating/deactivating elements
func (o *Domain) create_stage_copy() {
}

// set_act_deact_flags sets inactive flags for new active/inactive elements
func (o *Domain) fix_inact_flags(eids_or_tags []int, deactivate bool) {
	for _, tag := range eids_or_tags {
		if tag >= 0 { // this meahs that tag == cell.Id
			cell := o.Msh.Cells[tag]
			tag = cell.Tag
		}
		edat := o.Reg.Etag2data(tag)
		edat.Inact = deactivate
	}
}
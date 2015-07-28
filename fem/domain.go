// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
)

// Solution holds the solution data @ nodes.
//
//        / u \         / u \
//        |   | => y =  |   |
//  yb =  | p |         \ p / (ny x 1)
//        |   |
//        \ λ / (nyb x 1)
//
type Solution struct {

	// current state
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

	// problem definition and constants
	Steady  bool      // [from Sim] steady simulation
	Axisym  bool      // [from Sim] axisymmetric
	Pstress bool      // [from Sim] plane-stress
	DynCfs  *DynCoefs // [from FEM] coefficients for dynamics/transient simulations
}

// Domain holds all Nodes and Elements active during a stage in addition to the Solution at nodes.
// Only elements in this processor are recorded here; however information from
// all cells might be recorded as well.
type Domain struct {

	// init: auxiliary variables
	Distr  bool            // distributed/parallel run
	Proc   int             // this processor number
	Sim    *inp.Simulation // [from FEM] input data
	Reg    *inp.Region     // region data
	Msh    *inp.Mesh       // mesh data
	LinSol la.LinSol       // linear solver
	DynCfs *DynCoefs       // [from FEM] coefficients for dynamics/transient simulations
	HydSta *HydroStatic    // [from FEM] function to compute hydrostatic state

	// stage: nodes (active) and elements (active AND in this processor)
	Nodes  []*Node // active nodes (for each stage)
	Elems  []Elem  // [procNcells] only active elements in this processor (for each stage)
	MyCids []int   // [procNcells] the ids of cells in this processor

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

	// for divergence control
	bkpSol *Solution // backup solution
}

// Clean cleans memory allocated by domain
func (o *Domain) Clean() {
	o.LinSol.Clean()
}

// NewDomains returns domains
func NewDomains(sim *inp.Simulation, dyncfs *DynCoefs, hydsta *HydroStatic, proc, nproc int, distr bool) (doms []*Domain) {
	doms = make([]*Domain, len(sim.Regions))
	for i, reg := range sim.Regions {
		doms[i] = new(Domain)
		doms[i].Distr = distr
		doms[i].Proc = proc
		doms[i].Sim = sim
		doms[i].Reg = reg
		doms[i].Msh = reg.Msh
		if distr {
			if nproc != len(reg.Msh.Part2cells) {
				chk.Panic("number of processors must be equal to the number of partitions defined in mesh file. %d != %d", nproc, len(reg.Msh.Part2cells))
			}
		}
		doms[i].LinSol = la.GetSolver(sim.LinSol.Name)
		doms[i].DynCfs = dyncfs
		doms[i].HydSta = hydsta
	}
	return
}

// SetStage set nodes, equation numbers and auxiliary data for given stage
func (o *Domain) SetStage(stgidx int) (err error) {

	// pointer to stage structure
	stg := o.Sim.Stages[stgidx]

	// backup state
	if stgidx > 0 {
		o.create_stage_copy()
		err = o.fix_inact_flags(stg.Activate, false)
		if err != nil {
			return
		}
		err = o.fix_inact_flags(stg.Deactivate, true)
		if err != nil {
		}
	}

	// nodes (active) and elements (active AND in this processor)
	o.Nodes = make([]*Node, 0)
	o.Elems = make([]Elem, 0)
	o.MyCids = make([]int, 0)

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
	o.ElemIntvars = make([]ElemIntvars, 0)

	// allocate nodes and cells (active only) -------------------------------------------------------

	// for each cell
	var eq int // current equation number => total number of equations @ end of loop
	o.NnzKb = 0
	for _, cell := range o.Msh.Cells {

		// set cell's face boundary conditions
		cell.SetFaceConds(stg, o.Sim.Functions)

		// get element info
		info, inactive, err := GetElemInfo(cell, o.Reg, o.Sim)
		if err != nil {
			return chk.Err("get element information failed:\n%v", err)
		}

		// skip inactive element
		if inactive {
			continue
		}
		o.Cid2active[cell.Id] = true

		// for non-joint elements, add new DOFs
		if !cell.IsJoint {
			chk.IntAssert(len(info.Dofs), len(cell.Verts))

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
			for j, v := range cell.Verts {

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

			// number of non-zeros
			o.NnzKb += eNdof * eNdof
		}

		// allocate element
		mycell := cell.Part == o.Proc // cell belongs to this processor
		if mycell || !o.Distr {

			// new element
			ele, err := NewElem(cell, o.Reg, o.Sim)
			if err != nil {
				return chk.Err("new element failed:%\v", err)
			}
			o.Cid2elem[cell.Id] = ele
			o.Elems = append(o.Elems, ele)
			o.MyCids = append(o.MyCids, ele.Id())

			// give equation numbers to new element
			eqs := make([][]int, len(cell.Verts))
			for j, v := range cell.Verts {
				for _, dof := range o.Vid2node[v].Dofs {
					eqs[j] = append(eqs[j], dof.Eq)
				}
			}
			err = ele.SetEqs(eqs, nil)
			if err != nil {
				return chk.Err("cannot set element equations:\n%v", err)
			}

			// subsets of elements
			o.add_element_to_subsets(ele)
		}
	}

	// connect elements (e.g. Joints)
	for _, e := range o.ElemConnect {
		nnz, err := e.Connect(o.Cid2elem, o.Msh.Cells[e.Id()])
		if err != nil {
			return chk.Err("cannot connect rod-joint elements with solids and rods:\n%v", err)
		}
		o.NnzKb += nnz
	}

	// element conditions, essential and natural boundary conditions --------------------------------

	// (re)set constraints and prescribed forces structures
	o.EssenBcs.Init(o.HydSta)
	o.PtNatBcs.Reset()

	// element conditions
	for _, ec := range stg.EleConds {
		cells, ok := o.Msh.CellTag2cells[ec.Tag]
		if !ok {
			return chk.Err("cannot find cells with tag = %d to assign conditions", ec.Tag)
		}
		for _, cell := range cells {
			e := o.Cid2elem[cell.Id]
			if e != nil { // set conditions only for this processor's / active element
				for j, key := range ec.Keys {
					fcn := o.Sim.Functions.Get(ec.Funcs[j])
					if fcn == nil {
						return chk.Err("cannot find function named %q\n", ec.Funcs[j])
					}
					e.SetEleConds(key, fcn, ec.Extra)
				}
			}
		}
	}

	// face essential boundary conditions
	for _, cellsAndFaces := range o.Msh.FaceTag2cells {
		for _, pair := range cellsAndFaces {
			cell := pair.C
			for _, fc := range cell.FaceBcs {
				var enodes []*Node
				for _, v := range fc.GlobalVerts {
					enodes = append(enodes, o.Vid2node[v])
				}
				if o.YandC[fc.Cond] {
					err = o.EssenBcs.Set(fc.Cond, enodes, fc.Func, fc.Extra)
					if err != nil {
						return chk.Err("setting of essential boundary conditions failed:\n%v", err)
					}
				}
			}
		}
	}

	// vertex bounday conditions
	for _, nc := range stg.NodeBcs {
		verts, ok := o.Msh.VertTag2verts[nc.Tag]
		if !ok {
			return chk.Err("cannot find vertices with tag = %d to assign node boundary conditions", nc.Tag)
		}
		for _, v := range verts {
			if o.Vid2node[v.Id] != nil { // set BCs only for active nodes
				n := o.Vid2node[v.Id]
				for j, key := range nc.Keys {
					fcn := o.Sim.Functions.Get(nc.Funcs[j])
					if fcn == nil {
						return chk.Err("cannot find function named %q\n", nc.Funcs[j])
					}
					if o.YandC[key] {
						o.EssenBcs.Set(key, []*Node{n}, fcn, nc.Extra)
					} else {
						o.PtNatBcs.Set(o.F2Y[key], n, fcn, nc.Extra)
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
		for _, dof := range nod.Dofs {
			switch o.Dof2Tnum[dof.Key] {
			case 1:
				o.T1eqs = append(o.T1eqs, dof.Eq)
			case 2:
				o.T2eqs = append(o.T2eqs, dof.Eq)
			default:
				chk.Panic("t1 and t2 equations are incorrectly set")
			}
		}
	}

	// size of arrays
	o.Ny = eq
	o.Nlam, o.NnzA = o.EssenBcs.Build(o.Ny)
	o.Nyb = o.Ny + o.Nlam

	// solution structure and linear solver
	o.Sol = new(Solution)
	o.Sol.Steady = o.Sim.Data.Steady
	o.Sol.Axisym = o.Sim.Data.Axisym
	o.Sol.Pstress = o.Sim.Data.Pstress
	o.Sol.DynCfs = o.DynCfs

	// linear system and linear solver
	o.Kb = new(la.Triplet)
	o.Fb = make([]float64, o.Nyb)
	o.Wb = make([]float64, o.Nyb)
	o.Kb.Init(o.Nyb, o.Nyb, o.NnzKb+2*o.NnzA)
	o.InitLSol = true // tell solver that lis has to be initialised before use

	// allocate arrays
	o.Sol.Y = make([]float64, o.Ny)
	o.Sol.ΔY = make([]float64, o.Ny)
	o.Sol.L = make([]float64, o.Nlam)
	if !o.Sim.Data.Steady {
		o.Sol.Dydt = make([]float64, o.Ny)
		o.Sol.D2ydt2 = make([]float64, o.Ny)
		o.Sol.Psi = make([]float64, o.Ny)
		o.Sol.Zet = make([]float64, o.Ny)
		o.Sol.Chi = make([]float64, o.Ny)
	}

	// success
	return
}

// SetIniVals sets/resets initial values (nodes and integration points)
func (o *Domain) SetIniVals(stgidx int, zeroSol bool) (err error) {

	// pointer to stage structure
	stg := o.Sim.Stages[stgidx]

	// clear solution vectors
	if zeroSol {
		o.Sol.Reset(o.Sim.Data.Steady)
	}

	// initialise internal variables
	if stg.HydroSt {
		err = o.SetHydroSt(stg)
		if err != nil {
			return
		}
	} else if stg.GeoSt != nil {
		err = o.SetGeoSt(stg)
		if err != nil {
			return
		}
	} else if stg.IniStress != nil {
		err = o.SetIniStress(stg)
		if err != nil {
			return
		}
	} else if stg.Initial != nil {
		err = o.SetInitial(stg)
		if err != nil {
			return
		}
	} else {
		for _, e := range o.ElemIntvars {
			e.SetIniIvs(o.Sol, nil)
		}
	}

	// import results from another set of files
	if stg.Import != nil {
		sum := new(Summary)
		err = sum.Read(stg.Import.Dir, stg.Import.Fnk, o.Sim.EncType)
		if err != nil {
			return chk.Err("cannot import state from %s/%s.sim:\n%v", stg.Import.Dir, stg.Import.Fnk, err)
		}
		err = o.Read(sum, len(sum.OutTimes)-1, o.Proc, false)
		if err != nil {
			return chk.Err("cannot load results into domain:\n%v", err)
		}
		if o.Ny != len(o.Sol.Y) {
			return chk.Err("length of primary variables vector imported is not equal to the one allocated. make sure the number of DOFs of the imported simulation matches this one. %d != %d", o.Ny, len(o.Sol.Y))
		}
		if stg.Import.ResetU {
			for _, ele := range o.ElemIntvars {
				err = ele.Ureset(o.Sol)
				if err != nil {
					return chk.Err("cannot run reset function of element after displacements are zeroed:\n%v", err)
				}
			}
			for _, nod := range o.Nodes {
				for _, ukey := range []string{"ux", "uy", "uz"} {
					eq := nod.GetEq(ukey)
					if eq >= 0 {
						o.Sol.Y[eq] = 0
						if len(o.Sol.Dydt) > 0 {
							o.Sol.Dydt[eq] = 0
							o.Sol.D2ydt2[eq] = 0
						}
					}
				}
			}
		}
	}

	// make sure time is zero at the beginning of simulation
	o.Sol.T = 0
	return
}

// auxiliary functions //////////////////////////////////////////////////////////////////////////////

// Reset clear values
func (o *Solution) Reset(steady bool) {
	o.T = 0
	for i := 0; i < len(o.Y); i++ {
		o.Y[i] = 0
		o.ΔY[i] = 0
	}
	if !steady {
		for i := 0; i < len(o.Y); i++ {
			o.Psi[i] = 0
			o.Zet[i] = 0
			o.Chi[i] = 0
			o.Dydt[i] = 0
			o.D2ydt2[i] = 0
		}
	}
	for i := 0; i < len(o.L); i++ {
		o.L[i] = 0
	}
}

// add_element_to_subsets adds an Elem to many subsets as it fits
func (o *Domain) add_element_to_subsets(ele Elem) {
	if e, ok := ele.(ElemIntvars); ok {
		o.ElemIntvars = append(o.ElemIntvars, e)
	}
	if e, ok := ele.(ElemConnector); ok {
		o.ElemConnect = append(o.ElemConnect, e)
	}
}

// create_stage_copy creates a copy of current stage => to be used later when activating/deactivating elements
func (o *Domain) create_stage_copy() {
}

// set_act_deact_flags sets inactive flags for new active/inactive elements
func (o *Domain) fix_inact_flags(eids_or_tags []int, deactivate bool) (err error) {
	for _, tag := range eids_or_tags {
		if tag >= 0 { // this meahs that tag == cell.Id
			cell := o.Msh.Cells[tag]
			tag = cell.Tag
		}
		edat := o.Reg.Etag2data(tag)
		if edat == nil {
			return chk.Err("cannot get element's data with etag=%d", tag)
		}
		edat.Inact = deactivate
	}
	return
}

// backup saves a copy of solution
func (o *Domain) backup() {
	if o.bkpSol == nil {
		o.bkpSol = new(Solution)
		o.bkpSol.Y = make([]float64, o.Ny)
		o.bkpSol.ΔY = make([]float64, o.Ny)
		o.bkpSol.L = make([]float64, o.Nlam)
		if !o.Sim.Data.Steady {
			o.bkpSol.Dydt = make([]float64, o.Ny)
			o.bkpSol.D2ydt2 = make([]float64, o.Ny)
			o.bkpSol.Psi = make([]float64, o.Ny)
			o.bkpSol.Zet = make([]float64, o.Ny)
			o.bkpSol.Chi = make([]float64, o.Ny)
		}
	}
	o.bkpSol.T = o.Sol.T
	copy(o.bkpSol.Y, o.Sol.Y)
	copy(o.bkpSol.ΔY, o.Sol.ΔY)
	copy(o.bkpSol.L, o.Sol.L)
	if !o.Sim.Data.Steady {
		copy(o.bkpSol.Dydt, o.Sol.Dydt)
		copy(o.bkpSol.D2ydt2, o.Sol.D2ydt2)
		copy(o.bkpSol.Psi, o.Sol.Psi)
		copy(o.bkpSol.Zet, o.Sol.Zet)
		copy(o.bkpSol.Chi, o.Sol.Chi)
	}
	for _, e := range o.ElemIntvars {
		e.BackupIvs(true)
	}
}

// restore restores solution
func (o *Domain) restore() {
	o.Sol.T = o.bkpSol.T
	copy(o.Sol.Y, o.bkpSol.Y)
	copy(o.Sol.ΔY, o.bkpSol.ΔY)
	copy(o.Sol.L, o.bkpSol.L)
	if !o.Sim.Data.Steady {
		copy(o.Sol.Dydt, o.bkpSol.Dydt)
		copy(o.Sol.D2ydt2, o.bkpSol.D2ydt2)
		copy(o.Sol.Psi, o.bkpSol.Psi)
		copy(o.Sol.Zet, o.bkpSol.Zet)
		copy(o.Sol.Chi, o.bkpSol.Chi)
	}
	for _, e := range o.ElemIntvars {
		e.RestoreIvs(true)
	}
}

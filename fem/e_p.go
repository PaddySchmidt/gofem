// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// ElemP implements an element for transient seepage analyses [1]
//  References:
//   [1] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816,
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type ElemP struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Shp  *shp.Shape  // shape structure
	Np   int         // total number of unknowns == number of vertices
	Ndim int         // space dimension

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model
	Mdl *mporous.Model // model

	// problem variables
	Pmap []int // assembly map (location array/element equations)

	// internal variables
	States    []*mporous.State
	StatesBkp []*mporous.State
	StatesAux []*mporous.State

	// gravity
	Gfcn fun.Func // gravity function

	// natural boundary conditions
	NatBcs []*NaturalBc // natural boundary conditions

	// flux boundary conditions (qb == \bar{q})
	ρl_ex     []float64   // [nverts] ρl extrapolted to nodes => if has qb (flux)
	dρldpl_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Emat      [][]float64 // [nverts][nips] extrapolator matrix
	DoExtrap  bool        // do extrapolation of ρl and Cpl => for use with flux and seepage conditions

	// seepage face
	Nf         int         // number of fl variables
	HasSeep    bool        // indicates if this element has seepage faces
	Vid2seepId []int       // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int       // [nf] maps seepage face variable id to local vertex id
	Fmap       []int       // [nf] map of "fl" variables (seepage face)
	Macaulay   bool        // use discrete ramp function instead of smooth ramp
	βrmp       float64     // coefficient for Sramp
	κ          float64     // κ coefficient to normalise equation for seepage face modelling
	Hst        []bool      // [nf] set hydrostatic plmax
	Plmax      [][]float64 // [nf][nipsFace] specified plmax (not corrected by multiplier)

	// local starred variables
	ψl []float64 // [nip] ψl* = β1.p + β2.dpdt

	// scratchpad. computed @ each ip
	g   []float64       // [ndim] gravity vector
	pl  float64         // pl: liquid pressure
	gpl []float64       // [ndim] ∇pl: gradient of liquid pressure
	ρwl []float64       // [ndim] ρl*wl: weighted liquid relative velocity
	tmp []float64       // [ndim] temporary (auxiliary) vector
	Kpp [][]float64     // [np][np] Kpp := dRpl/dpl consistent tangent matrix
	Kpf [][]float64     // [np][nf] Kpf := dRpl/dfl consistent tangent matrix
	Kfp [][]float64     // [nf][np] Kfp := dRfl/dpl consistent tangent matrix
	Kff [][]float64     // [nf][nf] Kff := dRfl/dfl consistent tangent matrix
	res *mporous.LsVars // variable to hold results from CalcLs
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["p"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		// number of nodes in element
		nverts := shp.GetNverts(cell.Type)

		// solution variables
		ykeys := []string{"pl"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"pl": "ql"}

		// vertices on seepage faces
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("seep")
			for _, m := range lverts {
				if m < nverts { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					info.Dofs[m] = append(info.Dofs[m], "fl")
				}
			}
			if len(lverts) > 0 {
				ykeys = append(ykeys, "fl")
				info.Y2F["fl"] = "nil"
			}
		}

		// t1 and t2 variables
		info.T1vars = ykeys
		return &info
	}

	// element allocator
	eallocators["p"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemP
		o.Cell = cell
		o.X = x
		o.Shp = shp.Get(cell.Type, sim.GoroutineId)
		o.Np = o.Shp.Nverts
		o.Ndim = sim.Ndim

		// integration points
		var err error
		o.IpsElem, o.IpsFace, err = GetIntegrationPoints(edat.Nip, edat.Nipf, cell.Type)
		if err != nil {
			chk.Panic("cannot allocate integration points of p-element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}
		nip := len(o.IpsElem)

		// models
		o.Mdl, err = GetAndInitPorousModel(sim.MatParams, edat.Mat, sim.Key)
		if err != nil {
			chk.Panic("cannot get model for p-element {tag=%d id=%d material=%q}:\n%v", cell.Tag, cell.Id, edat.Mat, err)
		}

		// local starred variables
		o.ψl = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.g = make([]float64, o.Ndim)
		o.gpl = make([]float64, o.Ndim)
		o.ρwl = make([]float64, o.Ndim)
		o.tmp = make([]float64, o.Ndim)
		o.Kpp = la.MatAlloc(o.Np, o.Np)
		o.res = new(mporous.LsVars)

		// vertices on seepage faces
		var seepverts []int
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("seep")
			for _, m := range lverts {
				if m < o.Np { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					seepverts = append(seepverts, m)
				}
			}
		}

		o.Nf = len(seepverts)
		o.HasSeep = o.Nf > 0
		if o.HasSeep {

			// vertices on seepage face; numbering
			o.SeepId2vid = seepverts
			o.Vid2seepId = utl.IntVals(o.Np, -1)
			o.Fmap = make([]int, o.Nf)
			for μ, m := range o.SeepId2vid {
				o.Vid2seepId[m] = μ
			}

			// flags
			o.Macaulay, o.βrmp, o.κ = GetSeepFaceFlags(edat.Extra)

			// allocate coupling matrices
			o.Kpf = la.MatAlloc(o.Np, o.Nf)
			o.Kfp = la.MatAlloc(o.Nf, o.Np)
			o.Kff = la.MatAlloc(o.Nf, o.Nf)
		}

		// set natural boundary conditions
		for idx, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})

			// allocate extrapolation structures
			if fc.Cond == "ql" || fc.Cond == "seep" {
				nv := o.Shp.Nverts
				nip := len(o.IpsElem)
				o.ρl_ex = make([]float64, nv)
				o.dρldpl_ex = la.MatAlloc(nv, nv)
				o.Emat = la.MatAlloc(nv, nip)
				o.DoExtrap = true
				err = o.Shp.Extrapolator(o.Emat, o.IpsElem)
				if err != nil {
					chk.Panic("cannot build extrapolator matrix for p-element:\n%v", err)
				}
			}

			// additional seepage condition structures: hydrostatic flags
			if fc.Cond == "seep" {
				if len(o.Hst) == 0 {
					o.Hst = make([]bool, len(cell.FaceBcs))
				}
				if s_val, found := io.Keycode(fc.Extra, "plmax"); found {
					o.Hst[idx] = (s_val == "hst")
				}
			}
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElemP) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *ElemP) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Pmap = make([]int, o.Np)
	for m := 0; m < o.Shp.Nverts; m++ {
		o.Pmap[m] = eqs[m][0]
	}
	if o.HasSeep {
		for i, m := range o.SeepId2vid {
			o.Fmap[i] = eqs[m][1]
		}
	}
	return
}

// SetEleConds sets element conditions
func (o *ElemP) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemP) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.ψl[idx] = 0
		for m := 0; m < o.Shp.Nverts; m++ {
			o.ψl[idx] += o.Shp.S[m] * sol.Psi[o.Pmap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o ElemP) AddToRhs(fb []float64, sol *Solution) (err error) {

	// clear variables
	if o.DoExtrap {
		la.VecFill(o.ρl_ex, 0)
	}

	// for each integration point
	β1 := sol.DynCfs.β1
	nverts := o.Shp.Nverts
	var coef, plt, klr, ρL, ρl, Cpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Shp.J * ip.W
		S := o.Shp.S
		G := o.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		err = o.Mdl.CalcLs(o.res, o.States[idx], o.pl, 0, false)
		if err != nil {
			return
		}
		ρl = o.res.A_ρl
		Cpl = o.res.Cpl

		// compute ρwl. see Eq. (6) of [1]
		for i := 0; i < o.Ndim; i++ {
			o.ρwl[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.ρwl[i] += klr * o.Mdl.Klsat[i][j] * (ρL*o.g[j] - o.gpl[j])
			}
		}

		// add negative of residual term to fb. see Eqs. (12) and (17) of [1]
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			fb[r] -= coef * S[m] * Cpl * plt
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.ρwl[i] // += coef * div(ρl*wl)
			}
			if o.DoExtrap { // Eq. (19)
				o.ρl_ex[m] += o.Emat[m][idx] * ρl
			}
		}
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		return o.add_natbcs_to_rhs(fb, sol)
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o ElemP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// clear matrices
	la.MatFill(o.Kpp, 0)
	nverts := o.Shp.Nverts
	if o.DoExtrap {
		for i := 0; i < nverts; i++ {
			o.ρl_ex[i] = 0
			for j := 0; j < nverts; j++ {
				o.dρldpl_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	Cl := o.Mdl.Cl
	β1 := sol.DynCfs.β1
	var coef, plt, klr, ρL, ρl, Cpl, dCpldpl, dklrdpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Shp.J * ip.W
		S := o.Shp.S
		G := o.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		err = o.Mdl.CalcLs(o.res, o.States[idx], o.pl, 0, true)
		if err != nil {
			return
		}
		ρl = o.res.A_ρl
		Cpl = o.res.Cpl
		dCpldpl = o.res.DCpldpl
		dklrdpl = o.res.Dklrdpl

		// Kpp := dRpl/dpl. see Eqs. (18), (A.2) and (A.3) of [1]
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.tmp[j] = S[n]*dklrdpl*(ρL*o.g[j]-o.gpl[j]) + klr*(S[n]*Cl*o.g[j]-G[n][j])
			}
			for m := 0; m < nverts; m++ {
				o.Kpp[m][n] += coef * S[m] * S[n] * (dCpldpl*plt + β1*Cpl)
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.Kpp[m][n] -= coef * G[m][i] * o.Mdl.Klsat[i][j] * o.tmp[j]
					}
				}
				if o.DoExtrap { // inner summation term in Eq. (22)
					o.dρldpl_ex[m][n] += o.Emat[m][idx] * Cpl * S[n]
				}
			}
			if o.DoExtrap { // Eq. (19)
				o.ρl_ex[n] += o.Emat[n][idx] * ρl
			}
		}
	}

	// add to Kb
	if o.HasSeep {

		// contribution from natural boundary conditions
		err = o.add_natbcs_to_jac(sol)
		if err != nil {
			return
		}

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kpf[i][j])
				Kb.Put(J, I, o.Kfp[j][i])
			}
		}
		for i, I := range o.Fmap {
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kff[i][j])
			}
		}

	} else {

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
		}
	}
	return
}

// Update performs (tangent) update
func (o *ElemP) Update(sol *Solution) (err error) {

	// for each integration point
	var pl, Δpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		// compute pl and Δpl @ ip by means of interpolating from nodes
		pl, Δpl = 0, 0
		for m := 0; m < o.Shp.Nverts; m++ {
			r := o.Pmap[m]
			pl += o.Shp.S[m] * sol.Y[r]
			Δpl += o.Shp.S[m] * sol.ΔY[r]
		}

		// update state
		err = o.Mdl.Update(o.States[idx], Δpl, 0, pl, 0)
		if err != nil {
			return
		}
	}
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o ElemP) Ipoints() (coords [][]float64) {
	coords = la.MatAlloc(len(o.IpsElem), o.Ndim)
	for idx, ip := range o.IpsElem {
		coords[idx] = o.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *ElemP) SetIniIvs(sol *Solution, ignored map[string][]float64) (err error) {

	// auxiliary
	nip := len(o.IpsElem)
	nverts := o.Shp.Nverts
	var ρL, ρG, pl, pg float64

	// allocate slices of states
	o.States = make([]*mporous.State, nip)
	o.StatesBkp = make([]*mporous.State, nip)
	o.StatesAux = make([]*mporous.State, nip)

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		S := o.Shp.S
		G := o.Shp.G

		// interpolate pl variables
		pl = 0
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] = 0
		}
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			pl += S[m] * sol.Y[r]
			for i := 0; i < o.Ndim; i++ {
				o.gpl[i] += G[m][i] * sol.Y[r]
			}
		}

		// compute density from hydrostatic condition => enforce initial ρwl = 0
		ρL = o.Mdl.RhoL0
		o.compute_gvec(sol.T)
		if math.Abs(o.g[o.Ndim-1]) > 0 {
			ρL = o.gpl[o.Ndim-1] / o.g[o.Ndim-1]
		}

		// state initialisation
		o.States[idx], err = o.Mdl.NewState(ρL, ρG, pl, pg)
		if err != nil {
			return
		}

		// backup copy
		o.StatesBkp[idx] = o.States[idx].GetCopy()
		o.StatesAux[idx] = o.States[idx].GetCopy()
	}

	// seepage face structures
	if o.HasSeep {
		o.Plmax = la.MatAlloc(len(o.NatBcs), len(o.IpsFace))
		for idx, nbc := range o.NatBcs {
			iface := nbc.IdxFace
			for jdx, ipf := range o.IpsFace {
				err = o.Shp.CalcAtFaceIp(o.X, ipf, iface)
				if err != nil {
					return
				}
				Sf := o.Shp.Sf
				switch nbc.Key {
				case "seep":
					pl = 0
					for i, m := range o.Shp.FaceLocalV[iface] {
						pl += Sf[i] * sol.Y[o.Pmap[m]]
					}
					o.Plmax[idx][jdx] = pl
				}
			}
		}
	}
	return
}

// BackupIvs creates copy of internal variables
func (o *ElemP) BackupIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.StatesAux {
			s.Set(o.States[i])
		}
		return
	}
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return
}

// RestoreIvs restores internal variables from copies
func (o *ElemP) RestoreIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.States {
			s.Set(o.StatesAux[i])
		}
		return
	}
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *ElemP) Ureset(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemP) Encode(enc Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o ElemP) Decode(dec Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpsData returns data from all integration points for output
func (o ElemP) OutIpsData() (data []*OutIpData) {
	flow := FlowKeys(o.Ndim)
	for idx, ip := range o.IpsElem {
		s := o.States[idx]
		x := o.Shp.IpRealCoords(o.X, ip)
		calc := func(sol *Solution) (vals map[string]float64) {
			err := o.ipvars(idx, sol)
			if err != nil {
				return
			}
			ρL := s.A_ρL
			klr := o.Mdl.Cnd.Klr(s.A_sl)
			vals = map[string]float64{
				"sl": s.A_sl,
				"pl": o.pl,
				"nf": 1.0 - s.A_ns0,
			}
			for i := 0; i < o.Ndim; i++ {
				for j := 0; j < o.Ndim; j++ {
					vals[flow[i]] += klr * o.Mdl.Klsat[i][j] * (o.g[j] - o.gpl[j]/ρL)
				}
			}
			return
		}
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemP) ipvars(idx int, sol *Solution) (err error) {

	// interpolation functions and gradients
	err = o.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// auxiliary
	o.compute_gvec(sol.T)

	// clear pl and its gradient @ ip
	o.pl = 0
	for i := 0; i < o.Ndim; i++ {
		o.gpl[i] = 0
	}

	// compute pl and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.Shp.Nverts; m++ {
		r := o.Pmap[m]
		o.pl += o.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] += o.Shp.G[m][i] * sol.Y[r]
		}
	}
	return
}

// fipvars computes current values @ face integration points
func (o *ElemP) fipvars(fidx int, sol *Solution) (ρl, pl, fl float64) {
	Sf := o.Shp.Sf
	for i, m := range o.Shp.FaceLocalV[fidx] {
		μ := o.Vid2seepId[m]
		ρl += Sf[i] * o.ρl_ex[m]
		pl += Sf[i] * sol.Y[o.Pmap[m]]
		fl += Sf[i] * sol.Y[o.Fmap[μ]]
	}
	return
}

// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o ElemP) add_natbcs_to_rhs(fb []float64, sol *Solution) (err error) {

	// compute surface integral
	var tmp float64
	var ρl, pl, fl, plmax, g, rmp, rx, rf float64
	for idx, nbc := range o.NatBcs {

		// tmp := plmax shift or qlb
		tmp = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for jdx, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Shp.Sf
			Jf := la.VecNorm(o.Shp.Fnvec)
			coef := ipf.W * Jf

			// select natural boundary condition type
			switch nbc.Key {
			case "ql":
				// flux prescribed
				ρl = 0
				for i, m := range o.Shp.FaceLocalV[iface] {
					ρl += Sf[i] * o.ρl_ex[m]
				}
				for i, m := range o.Shp.FaceLocalV[iface] {
					fb[o.Pmap[m]] -= coef * ρl * tmp * Sf[i]
				}
			case "seep":

				// variables extrapolated to face
				ρl, pl, fl = o.fipvars(iface, sol)
				plmax = o.Plmax[idx][jdx] - tmp
				if plmax < 0 {
					plmax = 0
				}

				// compute residuals
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rx = ρl * rmp // Eq. (30)
				rf = fl - rmp // Eq. (26)
				for i, m := range o.Shp.FaceLocalV[iface] {
					μ := o.Vid2seepId[m]
					fb[o.Pmap[m]] -= coef * Sf[i] * rx
					fb[o.Fmap[μ]] -= coef * Sf[i] * rf
				}
			}
		}
	}
	return
}

// add_natbcs_to_jac adds contribution from natural boundary conditions to Jacobian
func (o ElemP) add_natbcs_to_jac(sol *Solution) (err error) {

	// clear matrices
	if o.HasSeep {
		for i := 0; i < o.Np; i++ {
			for j := 0; j < o.Nf; j++ {
				o.Kpf[i][j] = 0
				o.Kfp[j][i] = 0
			}
		}
		la.MatFill(o.Kff, 0)
	}

	// compute surface integral
	nverts := o.Shp.Nverts
	var shift float64
	var ρl, pl, fl, plmax, g, rmp, rmpD float64
	var drxdpl, drxdfl, drfdpl, drfdfl float64
	for idx, nbc := range o.NatBcs {

		// plmax shift
		shift = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for jdx, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Shp.Sf
			Jf := la.VecNorm(o.Shp.Fnvec)
			coef := ipf.W * Jf

			// select natural boundary condition type
			switch nbc.Key {
			case "seep":

				// variables extrapolated to face
				ρl, pl, fl = o.fipvars(iface, sol)
				plmax = o.Plmax[idx][jdx] - shift
				if plmax < 0 {
					plmax = 0
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rmpD = o.rampD1(fl + o.κ*g)
				drxdpl = ρl * o.κ * rmpD // first term in Eq. (A.4) (without Sn)
				drxdfl = ρl * rmpD       // Eq. (A.5) (without Sn)
				drfdpl = -o.κ * rmpD     // Eq. (A.6) (corrected with κ and without Sn)
				drfdfl = 1.0 - rmpD      // Eq. (A.7) (without Sn)
				for i, m := range o.Shp.FaceLocalV[iface] {
					μ := o.Vid2seepId[m]
					for j, n := range o.Shp.FaceLocalV[iface] {
						ν := o.Vid2seepId[n]
						o.Kpp[m][n] += coef * Sf[i] * Sf[j] * drxdpl
						o.Kpf[m][ν] += coef * Sf[i] * Sf[j] * drxdfl
						o.Kfp[μ][n] += coef * Sf[i] * Sf[j] * drfdpl
						o.Kff[μ][ν] += coef * Sf[i] * Sf[j] * drfdfl
					}
					for n := 0; n < nverts; n++ { // Eqs. (18) and (22)
						for l, r := range o.Shp.FaceLocalV[iface] {
							o.Kpp[m][n] += coef * Sf[i] * Sf[l] * o.dρldpl_ex[r][n] * rmp
						}
					}
				}
			}
		}
	}
	return
}

// ramp implements the ramp function
func (o *ElemP) ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.βrmp)
}

// rampderiv returns the ramp function first derivative
func (o *ElemP) rampD1(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.βrmp)
}

// compute_gvec computes gravity vector @ time t
func (o ElemP) compute_gvec(t float64) {
	o.g[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.g[o.Ndim-1] = -o.Gfcn.F(t, nil)
	}
}

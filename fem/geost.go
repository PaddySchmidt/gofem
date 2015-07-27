// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/ode"
)

// geostate holds state @ top of layer
type geostate struct{ pl, ρL, ρ, σV float64 }

// GeoLayer holds information of one soil layer. It computes pressures (σVabs, pl) and
// densities (ρL, ρ) based on the following model (fully liquid saturated)
//
//    ρL  = ρL0 + Cl・pl   thus   dρL/dpl = Cl
//    sl  = 1
//    ρ   = nf・sl・ρL + (1 - nf)・ρS
//    ns  = 1 - nf
//
//    Z(z) = zmax + T・(z - zmax)   with 0 ≤ T ≤ 1
//    dZ   = (z - zmax)・dT
//    dpl  = ρL(pl)・g・(-dZ)
//    dpl  = ρL(pl)・g・(zmax - z)・dT
//    dσV  = ρ(pl)・g・(zmax - z)・dT
//    Δz   = zmax - z
//
//            / dpl/dT \   / ρL(pl)・g・Δz  \
//    dY/dT = | dρL/dT | = | Cl・dpl/dT     |
//            | dρ/dT  |   | nf・sl・dρL/dT |
//            \ dσV/dT /   \ ρ(pl)・g・Δz   /
//
type GeoLayer struct {
	Tags  []int      // tags of cells within this layer
	Zmin  float64    // coordinate (elevation) at bottom of layer
	Zmax  float64    // coordinate (elevation) at top of layer
	Nodes []*Node    // nodes in layer
	Elems []Elem     // elements in layer
	Cl    float64    // liquid compressibility
	RhoS0 float64    // initial density of solids
	nf0   float64    // initial (constant) porosity
	K0    float64    // earth-pressure at rest
	Dpl   float64    // liquid pressure added by this layer
	DsigV float64    // absolute value of vertical stress increment added by this layer
	top   *geostate  // state @ top of layer
	fcn   ode.Cb_fcn // function for ode solver
	Jac   ode.Cb_jac // Jacobian for ode solver
	sol   ode.ODE    // ode solver
}

// Start starts ODE solver for computing state variables in Calc
//  prev -- previous state @ top of this layer
func (o *GeoLayer) Start(prev *geostate, g float64) {

	// set state @ top
	o.top = prev

	// y := {pl, ρL, ρ, σV} == geostate
	nf := o.nf0
	sl := 1.0
	o.fcn = func(f []float64, x float64, y []float64, args ...interface{}) error {
		Δz := args[0].(float64)
		ρL := y[1]
		ρ := y[2]
		f[0] = ρL * g * Δz    // dpl/dT
		f[1] = o.Cl * f[0]    // dρL/dT
		f[2] = nf * sl * f[1] // dρ/dT
		f[3] = ρ * g * Δz     // dσV/dT
		return nil
	}

	// set ODE (using numerical Jacobian)
	silent := true
	o.sol.Init("Radau5", 4, o.fcn, nil, nil, nil, silent)
	o.sol.Distr = false // must be sure to disable this; otherwise it causes problems in parallel runs
}

// Calc computes state @ level z
func (o GeoLayer) Calc(z float64) (*geostate, error) {
	y := []float64{o.top.pl, o.top.ρL, o.top.ρ, o.top.σV}
	Δz := o.Zmax - z
	err := o.sol.Solve(y, 0, 1, 1, false, Δz)
	if err != nil {
		err = chk.Err("geost: failed when calculating state sing ODE solver: %v", err)
		return nil, err
	}
	return &geostate{y[0], y[1], y[2], y[3]}, nil
}

// GeoLayers is a set of Layer
type GeoLayers []*GeoLayer

// Len the length of Layers
func (o GeoLayers) Len() int {
	return len(o)
}

// Swap swaps two Layers
func (o GeoLayers) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Layers: sort from top to bottom
func (o GeoLayers) Less(i, j int) bool {
	return o[i].Zmin > o[j].Zmin
}

// SetGeoSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (err error) {

	// check layers definition
	geo := stg.GeoSt
	if len(geo.Layers) < 1 {
		return chk.Err("geost: layers must be defined by specifying what tags belong to which layer")
	}

	// get region
	if len(o.Sim.Regions) != 1 {
		return chk.Err("geost: can only handle one domain for now")
	}
	reg := o.Sim.Regions[0]

	// gravity
	grav := o.Sim.Gravity.F(0, nil)

	// fix UseK0
	nlayers := len(geo.Layers)
	if len(geo.UseK0) != nlayers {
		geo.UseK0 = make([]bool, nlayers)
	}

	// initialise layers
	var L GeoLayers
	L = make([]*GeoLayer, nlayers)
	ndim := o.Sim.Ndim
	nodehandled := make(map[int]bool)
	ctaghandled := make(map[int]bool) // required to make sure all elements were initialised
	for i, tags := range geo.Layers {

		// new layer
		L[i] = new(GeoLayer)
		L[i].Tags = tags
		L[i].Zmin = o.Sim.MaxElev
		L[i].Zmax = 0
		L[i].Cl = o.Sim.WaterRho0 / o.Sim.WaterBulk

		// get porous parameters
		L[i].RhoS0, L[i].nf0, err = get_porous_parameters(o.Sim.MatParams, reg, tags[0])
		if err != nil {
			return
		}

		// parameters
		if geo.UseK0[i] {
			L[i].K0 = geo.K0[i]
		} else {
			L[i].K0 = geo.Nu[i] / (1.0 - geo.Nu[i])
		}
		if L[i].K0 < 1e-7 {
			return chk.Err("geost: K0 or Nu is incorect: K0=%g, Nu=%g", L[i].K0, geo.Nu)
		}

		// for each tag of cells in this layer
		for _, tag := range tags {

			// check tags
			cells := o.Msh.CellTag2cells[tag]
			if len(cells) < 1 {
				return chk.Err("geost: there are no cells with tag = %d", tag)
			}

			// set nodes and elements and find min and max z-coordinates
			for _, c := range cells {
				L[i].Elems = append(L[i].Elems, o.Cid2elem[c.Id])
				for _, v := range c.Verts {
					if !nodehandled[v] {
						L[i].Nodes = append(L[i].Nodes, o.Vid2node[v])
					}
					L[i].Zmin = min(L[i].Zmin, o.Msh.Verts[v].C[ndim-1])
					L[i].Zmax = max(L[i].Zmax, o.Msh.Verts[v].C[ndim-1])
					nodehandled[v] = true
				}
				ctaghandled[c.Tag] = true
			}
		}
	}

	// make sure all elements tags were handled
	for tag, _ := range o.Msh.CellTag2cells {
		if !ctaghandled[tag] {
			return chk.Err("geost: there are cells not included in any layer: ctag=%d", tag)
		}
	}

	// sort layers from top to bottom
	sort.Sort(L)

	// set previous/top states in layers and compute Sol.Y
	for i, lay := range L {

		// previous state
		var top *geostate
		ρS := lay.RhoS0
		nf := lay.nf0
		sl := 1.0
		if i == 0 {
			pl := (o.Sim.WaterLevel - o.Sim.MaxElev) * o.Sim.WaterRho0 * grav
			ρL := o.Sim.WaterRho0
			ρ := nf*sl*ρL + (1.0-nf)*ρS
			σV := -o.Sim.Data.Surch
			top = &geostate{pl, ρL, ρ, σV}
		} else {
			top, err = L[i-1].Calc(L[i-1].Zmin)
			if err != nil {
				return chk.Err("cannot compute state @ bottom of layer:\n%v", err)
			}
			ρL := top.ρL
			top.ρ = nf*sl*ρL + (1.0-nf)*ρS
		}

		// start layer
		lay.Start(top, grav)

		// set nodes
		for _, nod := range lay.Nodes {
			z := nod.Vert.C[ndim-1]
			s, err := lay.Calc(z)
			if err != nil {
				return chk.Err("cannot compute state @ node z = %g\n%v", z, err)
			}
			dof := nod.GetDof("pl")
			if dof != nil {
				o.Sol.Y[dof.Eq] = s.pl
			}
		}

		// set elements
		for _, elem := range lay.Elems {
			if ele, okk := elem.(ElemIntvars); okk {

				// build slices
				coords := ele.Ipoints()
				nip := len(coords)
				svT := make([]float64, nip) // total vertical stresses
				for i := 0; i < nip; i++ {
					z := coords[i][ndim-1]
					s, err := lay.Calc(z)
					if err != nil {
						return chk.Err("cannot compute state @ ip z = %g\n%v", z, err)
					}
					svT[i] = -s.σV
				}
				ivs := map[string][]float64{"svT": svT, "K0": []float64{lay.K0}}

				// set element's states
				err = ele.SetIniIvs(o.Sol, ivs)
				if err != nil {
					return chk.Err("geost: element's internal values setting failed:\n%v", err)
				}
			}
		}
	}
	return
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////////

// get_porous_parameters extracts parameters based on region data
func get_porous_parameters(mdb *inp.MatDb, reg *inp.Region, ctag int) (RhoS0, nf0 float64, err error) {
	edat := reg.Etag2data(ctag)
	mat := mdb.Get(edat.Mat)
	if mat.Model != "group" {
		err = chk.Err("geost: material type describing layer must be 'group' with porous data")
		return
	}
	if matname, found := io.Keycode(mat.Extra, "p"); found {
		m := mdb.Get(matname)
		for _, p := range m.Prms {
			switch p.N {
			case "RhoS0":
				RhoS0 = p.V
			case "nf0":
				nf0 = p.V
			}
		}
	}
	if RhoS0 < 1e-7 {
		err = chk.Err("geost: initial density of solids RhoS0=%g is incorrect", RhoS0)
		return
	}
	if nf0 < 1e-7 {
		err = chk.Err("geost: initial porosity nf0=%g is incorrect", nf0)
	}
	return
}

// String prints geostate
func (o *geostate) String() string {
	return io.Sf("pl=%g ρL=%g ρ=%g σV=%g\n", o.pl, o.ρL, o.ρ, o.σV)
}

// String prints a json formatted string with GeoLayers' content
func (o GeoLayers) String() string {
	if len(o) == 0 {
		return "[]"
	}
	l := "[\n"
	for i, lay := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  { \"Tags\":%v, \"Zmin\":%g, \"Zmax\":%g, \"nf0\":%g, \"RhoS0\":%g, \"Cl\":%g\n", lay.Tags, lay.Zmin, lay.Zmax, lay.nf0, lay.RhoS0, lay.Cl)
		l += "    \"Nodes\":["
		for j, nod := range lay.Nodes {
			if j > 0 {
				l += ","
			}
			l += io.Sf("%d", nod.Vert.Id)
		}
		l += "],\n    \"Elems\":["
		for j, ele := range lay.Elems {
			if j > 0 {
				l += ","
			}
			l += io.Sf("%d", ele.Id())
		}
		l += "] }"
	}
	l += "\n]"
	return l
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"math"

	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

const DIST_TOL = 1e-6 // tolerance to compare distances or proximity between points

// PlaneData holds data for handling planes described by {u,v}-coordinates
type PlaneData struct {
	Plane int          // Plane indicator: {0,1,2} == {x-pane, y-plane, z-plane}. plane perpendicular to {x,y,z}
	Ids   map[int]bool // Ids in Dom.Msh.Verts or Ipoints of notes/ips on plane
	Dx    []float64    // Increments between cells in global coordinates
	Iu    []int        // Index of global coordinates corresponding to {u,v} plane
	Du    []float64    // Increments between cells in {u,v} coordinates
	Umin  []float64    // min {u,v} coordinates
	Umax  []float64    // max {u,v} coordinates
	Nu    []int        // number of divisions along {u,v}
	Ubins gm.Bins      // bins to search points in {u,v} grid
	F     [][]float64  // [nu][nv] values of a function evaluated at the uv coords; i.e. f(u,v)

	// derived (set by ConnectResults)
	id2pt map[int]*Point // results: connects Ids to Points in Results structure
}

// Locate finds points on plane => implement the Locator interface
func (o PlaneData) Locate() (res Points) {
	for id, _ := range o.Ids {
		q := get_nod_point(id, nil) // nodes
		if q != nil {
			res = append(res, q)
		} else {
			q := get_ip_point(id, nil) // integration points
			if q != nil {
				res = append(res, q)
			}
		}
	}
	return
}

// ConnectResults connects Ids to Points in Results structure
func (o *PlaneData) ConnectResults(alias string) {
	pts, ok := Results[alias]
	if !ok {
		chk.Panic("cannot get points/results with alias=%q", alias)
	}
	o.id2pt = make(map[int]*Point)
	for _, p := range pts {
		if p.Vid >= 0 {
			if !o.Ids[p.Vid] {
				chk.Panic("Vertex with vid=%d in Results[%q] is not on plane", p.Vid, alias)
			}
			o.id2pt[p.Vid] = p
		}
		if p.IpId >= 0 {
			if !o.Ids[p.IpId] {
				chk.Panic("Integration point with ipid=%d in Results[%q] is not on plane", p.IpId, alias)
			}
			o.id2pt[p.IpId] = p
		}
	}
}

// NodesOnPlane finds vertices located on {x,y,z} plane and returns an iterator
// that can be used to extract results on the plane. An {u,v}-coordinates system is built
// on the plane.
//  Input:
//    ftag -- cells' face tag; e.g. -31
//  Output:
//    dat -- plane data holding vertices on plane and other information; see above
//  Note:
//    1) this function works with 3D cells only
//    2) faces on cells must be either "qua4" or "qua8" types
//    3) middle nodes in "qua8" are disregarded in order to buidl an {u,v} grid
//    4) the resulting mesh on plane must be equally spaced; i.e. Δx and Δy are constant;
//       but not necessarily equal to each other
func NodesOnPlane(ftag int) *PlaneData {

	// check
	ndim := Dom.Msh.Ndim
	if ndim != 3 {
		chk.Panic("this function works in 3D only")
	}

	// loop over cells on face with given face tag
	var dat PlaneData
	dat.Plane = -1
	dat.Ids = make(map[int]bool)
	dat.Dx = make([]float64, ndim)
	Δx := make([]float64, ndim)
	first := true
	cells := Dom.Msh.FaceTag2cells[ftag]
	for _, cell := range cells {
		ctype := cell.C.Type
		for fidx, cftag := range cell.C.FTags {
			if cftag == ftag {

				// check face type
				ftype := shp.GetFaceType(ctype)
				if !(ftype == "qua4" || ftype == "qua8") {
					chk.Panic("can only handle qua4 or qua8 faces for now. ftype=%q", ftype)
				}

				// vertices on face
				flvids := shp.GetFaceLocalVerts(ctype, fidx)
				nv := len(flvids)
				if nv == 8 {
					nv = 4 // avoid middle nodes
				}
				for i := 0; i < nv; i++ {
					vid := cell.C.Verts[flvids[i]]
					dat.Ids[vid] = true
				}

				// compute and check increments in global coordinates
				for i := 1; i < nv; i++ {
					a := cell.C.Verts[flvids[i]]
					b := cell.C.Verts[flvids[i-1]]
					xa := Dom.Msh.Verts[a].C
					xb := Dom.Msh.Verts[b].C
					for j := 0; j < ndim; j++ {
						Δx[j] = utl.Max(Δx[j], math.Abs(xa[j]-xb[j]))
					}
				}
				if first {
					for j := 0; j < ndim; j++ {
						dat.Dx[j] = Δx[j]
					}
					first = false
				} else {
					for j := 0; j < ndim; j++ {
						if math.Abs(dat.Dx[j]-Δx[j]) > 1e-10 {
							chk.Panic("all faces must have the same Δx,Δy,Δz")
						}
					}
				}

				// find which plane vertices are located on => which plane this face is parallel to
				perpto := -1 // "perpendicular to" indicator
				switch {
				case Δx[0] < DIST_TOL: // plane perpendicular to x-axis
					perpto = 0
				case Δx[1] < DIST_TOL: // plane perpendicular to y-axis
					perpto = 1
				case Δx[2] < DIST_TOL: // plane perpendicular to z-axis
					perpto = 2
				default:
					chk.Panic("planes must be perpendicular to one of the x-y-z axes")
				}
				if dat.Plane < 0 {
					dat.Plane = perpto
				} else {
					if perpto != dat.Plane {
						chk.Panic("all planes must be perperdicular to the same axis")
					}
				}
			}
		}
	}

	// uv: indices and increments
	dat.Iu = make([]int, 2)
	dat.Du = make([]float64, 2)
	switch dat.Plane {
	case 0: // plane perpendicular to x-axis
		dat.Iu = []int{1, 2}
	case 1: // plane perpendicular to y-axis
		dat.Iu = []int{0, 2}
	case 2: // plane perpendicular to z-axis
		dat.Iu = []int{0, 1}
	}
	for j := 0; j < 2; j++ {
		dat.Du[j] = dat.Dx[dat.Iu[j]]
	}

	// uv: limits
	dat.Umin = make([]float64, 2)
	dat.Umax = make([]float64, 2)
	first = true
	for vid, _ := range dat.Ids {
		x := Dom.Msh.Verts[vid].C
		if first {
			for j := 0; j < 2; j++ {
				dat.Umin[j] = x[dat.Iu[j]]
				dat.Umax[j] = x[dat.Iu[j]]
			}
			first = false
		} else {
			for j := 0; j < 2; j++ {
				dat.Umin[j] = utl.Min(dat.Umin[j], x[dat.Iu[j]])
				dat.Umax[j] = utl.Max(dat.Umax[j], x[dat.Iu[j]])
			}
		}
	}

	// uv: size of grid
	dat.Nu = make([]int, 2)
	for j := 0; j < 2; j++ {
		dat.Nu[j] = int((dat.Umax[j]-dat.Umin[j])/dat.Du[j]) + 1
	}

	// uv: bins
	dd := DIST_TOL * 2
	nb := 20
	dat.Ubins.Init([]float64{dat.Umin[0] - dd, dat.Umin[1] - dd}, []float64{dat.Umax[0] + dd, dat.Umax[1] + dd}, nb)
	u := []float64{0, 0}
	for vid, _ := range dat.Ids {
		x := Dom.Msh.Verts[vid].C
		for j := 0; j < 2; j++ {
			u[j] = x[dat.Iu[j]]
		}
		err := dat.Ubins.Append(u, vid)
		if err != nil {
			chk.Panic("cannot append {u,v} coordinate to bins. u=%v", u)
		}
	}

	// allocate F
	dat.F = la.MatAlloc(dat.Nu[0], dat.Nu[1])
	return &dat
}

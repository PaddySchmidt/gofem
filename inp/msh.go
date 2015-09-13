// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"math"
	"path/filepath"

	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// constants
const Ztol = 1e-7

// Vert holds vertex data
type Vert struct {
	Id  int       // id
	Tag int       // tag
	C   []float64 // coordinates (size==2 or 3)
}

// Cell holds cell data
type Cell struct {

	// input data
	Id     int    // id
	Tag    int    // tag
	Geo    int    // geometry type (gemlab code)
	Type   string // geometry type (string)
	Part   int    // partition id
	Verts  []int  // vertices
	FTags  []int  // edge (2D) or face (3D) tags
	STags  []int  // seam tags (for 3D only; it is actually a 3D edge tag)
	JlinId int    // joint line id
	JsldId int    // joint solid id

	// neighbours
	Neighs []int // neighbours; e.g. [3, 7, -1, 11] => side:cid => 0:3, 1:7, 2:-1(no cell), 3:11

	// derived
	Shp         *shp.Shape // shape structure
	FaceBcs     FaceConds  // face boundary condition
	GoroutineId int        // go routine id

	// specific problems data
	IsJoint   bool         // cell represents joint element
	SeepVerts map[int]bool // local vertices ids of vertices on seepage faces
	LbbCell   *Cell        // copy of this Cell with less vertices if LBB (externally allocated)

	// NURBS
	Nrb  int   // index of NURBS patch to which this cell belongs to
	Span []int // knot indices indicating which span this cell is located
}

// CellFaceId structure
type CellFaceId struct {
	C   *Cell // cell
	Fid int   // face id
}

// CellSeamId structure
type CellSeamId struct {
	C   *Cell // cell
	Sid int   // seam id
}

// Mesh holds a mesh for FE analyses
type Mesh struct {

	// from JSON
	Verts []*Vert // vertices
	Cells []*Cell // cells

	// derived
	FnamePath  string  // complete filename path
	Ndim       int     // space dimension
	Xmin, Xmax float64 // min and max x-coordinate
	Ymin, Ymax float64 // min and max x-coordinate
	Zmin, Zmax float64 // min and max x-coordinate

	// derived: maps
	VertTag2verts map[int][]*Vert      // vertex tag => set of vertices
	CellTag2cells map[int][]*Cell      // cell tag => set of cells
	FaceTag2cells map[int][]CellFaceId // face tag => set of cells
	FaceTag2verts map[int][]int        // face tag => vertices on tagged face
	SeamTag2cells map[int][]CellSeamId // seam tag => set of cells
	Ctype2cells   map[string][]*Cell   // cell type => set of cells
	Part2cells    map[int][]*Cell      // partition number => set of cells

	// NURBS
	Nurbss   []gm.NurbsD   // all NURBS data (read from file)
	PtNurbs  []*gm.Nurbs   // all NURBS' structures (allocated here)
	NrbFaces [][]*gm.Nurbs // all NURBS's faces
}

// ReadMsh reads a mesh for FE analyses
//  Note: returns nil on errors
func ReadMsh(dir, fn string, goroutineId int) (o *Mesh, err error) {

	// new mesh
	o = new(Mesh)

	// read file
	o.FnamePath = filepath.Join(dir, fn)
	b, err := io.ReadFile(o.FnamePath)
	if err != nil {
		return
	}

	// decode
	err = json.Unmarshal(b, &o)
	if err != nil {
		return
	}

	// check
	if len(o.Verts) < 2 {
		err = chk.Err("at least 2 vertices are required in mesh\n")
		return
	}
	if len(o.Cells) < 1 {
		err = chk.Err("at least 1 cell is required in mesh\n")
		return
	}

	// variables for NURBS
	var controlpts [][]float64
	has_nurbs := false
	if len(o.Nurbss) > 0 {
		controlpts = make([][]float64, len(o.Verts))
		has_nurbs = true
	}

	// vertex related derived data
	o.Ndim = 2
	o.Xmin = o.Verts[0].C[0]
	o.Ymin = o.Verts[0].C[1]
	if len(o.Verts[0].C) > 2 {
		o.Zmin = o.Verts[0].C[2]
	}
	o.Xmax = o.Xmin
	o.Ymax = o.Ymin
	o.Zmax = o.Zmin
	o.VertTag2verts = make(map[int][]*Vert)
	for i, v := range o.Verts {

		// check vertex id
		if v.Id != i {
			err = chk.Err("vertices ids must coincide with order in \"verts\" list. %d != %d\n", v.Id, i)
			return
		}

		// ndim
		nd := len(v.C)
		if nd < 2 || nd > 4 {
			err = chk.Err("number of space dimensions must be 2, 3 or 4 (NURBS). %d is invalid\n", nd)
			return
		}
		if nd == 3 {
			if math.Abs(v.C[2]) > Ztol {
				o.Ndim = 3
			}
		}

		// tags
		if v.Tag < 0 {
			verts := o.VertTag2verts[v.Tag]
			o.VertTag2verts[v.Tag] = append(verts, v)
		}

		// limits
		o.Xmin = utl.Min(o.Xmin, v.C[0])
		o.Xmax = utl.Max(o.Xmax, v.C[0])
		o.Ymin = utl.Min(o.Ymin, v.C[1])
		o.Ymax = utl.Max(o.Ymax, v.C[1])
		if nd > 2 {
			o.Zmin = utl.Min(o.Zmin, v.C[2])
			o.Zmax = utl.Max(o.Zmax, v.C[2])
		}

		// control points to initialise NURBS
		if has_nurbs {
			controlpts[i] = make([]float64, 4)
			for j := 0; j < 4; j++ {
				controlpts[i][j] = v.C[j]
			}
		}
	}

	// allocate NURBSs
	o.PtNurbs = make([]*gm.Nurbs, len(o.Nurbss))
	o.NrbFaces = make([][]*gm.Nurbs, len(o.Nurbss))
	for i, d := range o.Nurbss {
		o.PtNurbs[i] = new(gm.Nurbs)
		o.PtNurbs[i].Init(d.Gnd, d.Ords, d.Knots)
		o.PtNurbs[i].SetControl(controlpts, d.Ctrls)
		o.NrbFaces[i] = o.PtNurbs[i].ExtractSurfaces()
	}

	// derived data
	o.CellTag2cells = make(map[int][]*Cell)
	o.FaceTag2cells = make(map[int][]CellFaceId)
	o.FaceTag2verts = make(map[int][]int)
	o.SeamTag2cells = make(map[int][]CellSeamId)
	o.Ctype2cells = make(map[string][]*Cell)
	o.Part2cells = make(map[int][]*Cell)
	for i, c := range o.Cells {

		// check id and tag
		if c.Id != i {
			err = chk.Err("cells ids must coincide with order in \"verts\" list. %d != %d\n", c.Id, i)
			return
		}
		if c.Tag >= 0 {
			err = chk.Err("cells tags must be negative. %d is incorrect\n", c.Tag)
			return
		}

		// get shape structure
		switch c.Type {
		case "joint":
			c.IsJoint = true
		case "nurbs":
			c.Shp = shp.GetShapeNurbs(o.PtNurbs[c.Nrb], o.NrbFaces[c.Nrb], c.Span)
			if c.Shp == nil {
				err = chk.Err("cannot allocate \"shape\" structure for cell type = %q\n", c.Type)
				return
			}
		default:
			c.Shp = shp.Get(c.Type, goroutineId)
			if c.Shp == nil {
				err = chk.Err("cannot allocate \"shape\" structure for cell type = %q\n", c.Type)
				return
			}
		}
		c.GoroutineId = goroutineId

		// face tags
		cells := o.CellTag2cells[c.Tag]
		o.CellTag2cells[c.Tag] = append(cells, c)
		for i, ftag := range c.FTags {
			if ftag < 0 {
				pairs := o.FaceTag2cells[ftag]
				o.FaceTag2cells[ftag] = append(pairs, CellFaceId{c, i})
				for _, l := range c.Shp.FaceLocalVerts[i] {
					utl.IntIntsMapAppend(&o.FaceTag2verts, ftag, o.Verts[c.Verts[l]].Id)
				}
			}
		}

		// seam tags
		if o.Ndim == 3 {
			for i, stag := range c.STags {
				if stag < 0 {
					pairs := o.SeamTag2cells[stag]
					o.SeamTag2cells[stag] = append(pairs, CellSeamId{c, i})
				}
			}
		}

		// cell type => cells
		cells = o.Ctype2cells[c.Type]
		o.Ctype2cells[c.Type] = append(cells, c)

		// partition => cells
		cells = o.Part2cells[c.Part]
		o.Part2cells[c.Part] = append(cells, c)
	}

	// remove duplicates
	for ftag, verts := range o.FaceTag2verts {
		o.FaceTag2verts[ftag] = utl.IntUnique(verts)
	}

	// results
	return
}

// String returns a JSON representation of *Vert
func (o *Vert) String() string {
	l := io.Sf("{\"id\":%4d, \"tag\":%6d, \"c\":[", o.Id, o.Tag)
	for i, x := range o.C {
		if i > 0 {
			l += ", "
		}
		l += io.Sf("%23.15e", x)
	}
	l += "] }"
	return l
}

// String returns a JSON representation of *Cell
func (o *Cell) String() string {
	l := io.Sf("{\"id\":%d, \"tag\":%d, \"type\":%q, \"part\":%d, \"verts\":[", o.Id, o.Tag, o.Type, o.Part)
	for i, x := range o.Verts {
		if i > 0 {
			l += ", "
		}
		l += io.Sf("%d", x)
	}
	l += "], \"ftags\":["
	for i, x := range o.FTags {
		if i > 0 {
			l += ", "
		}
		l += io.Sf("%d", x)
	}
	l += "] }"
	return l
}

// AllocLbb allocates Lbb cell
func (o *Cell) AllocLbb() {
	o.LbbCell = new(Cell)

	// input data
	o.LbbCell.Id = o.Id
	o.LbbCell.Tag = o.Tag
	o.LbbCell.Geo = o.Geo
	o.LbbCell.Type = o.Type
	o.LbbCell.Part = o.Part
	o.LbbCell.Verts = o.Verts
	o.LbbCell.FTags = o.FTags
	o.LbbCell.STags = o.STags
	o.LbbCell.JlinId = o.JlinId
	o.LbbCell.JsldId = o.JsldId

	// neighbours
	o.LbbCell.Neighs = o.Neighs

	// derived
	if o.Shp.Nurbs == nil {
		o.LbbCell.Shp = shp.Get(o.Shp.BasicType, o.GoroutineId)
		if o.LbbCell.Shp == nil {
			chk.Panic("cannot allocate \"shape\" structure for cell type = %q\n", o.Shp.BasicType)
		}
	} else {
		chk.Panic("cannot handle LBB cells with NURBS yet")
	}
	o.LbbCell.FaceBcs = o.FaceBcs
	o.LbbCell.GoroutineId = o.GoroutineId

	// specific problems data
	o.LbbCell.IsJoint = o.IsJoint
	o.LbbCell.SeepVerts = o.SeepVerts

	// NURBS
	o.LbbCell.Nrb = o.Nrb
	o.LbbCell.Span = o.Span
	return
}

// String returns a JSON representation of *Mesh
func (o Mesh) String() string {
	l := "{\n  \"verts\" : [\n"
	for i, x := range o.Verts {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("    %v", x)
	}
	l += "\n  ],\n  \"cells\" : [\n"
	for i, x := range o.Cells {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("    %v", x)
	}
	l += "\n  ]\n}"
	return l
}

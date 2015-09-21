// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strings"

	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
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

// GetNverts returns the number of vertices, whether LBB condition is on or not
func (o *Cell) GetNverts(lbb bool) int {
	if o.Type == "joint" {
		return len(o.Verts)
	}
	if lbb {
		return o.Shp.BasicNverts
	}
	return o.Shp.Nverts
}

// GetVtkInfo returns information about this cell for generating VTK files
func (o *Cell) GetVtkInfo(lbb bool) (nvtkverts, vtkcode int) {
	if o.Type == "joint" {
		nvtkverts = len(o.Verts)
		vtkcode = shp.VTK_POLY_VERTEX
		return
	}
	if lbb {
		nvtkverts = o.Shp.BasicNverts
		vtkcode = o.Shp.BasicVtkCode
		return
	}
	nvtkverts = o.Shp.VtkNverts
	vtkcode = o.Shp.VtkCode
	return
}

// GetSimilar allocates a copy of this cell
//  Note: the resulting cell with share the same slices as the original one => not a full copy
func (o *Cell) GetSimilar(lbb bool) (newcell *Cell) {

	// new cell
	newcell = new(Cell)

	// input data
	newcell.Id = o.Id
	newcell.Tag = o.Tag
	newcell.Geo = o.Geo
	newcell.Type = o.Type
	newcell.Part = o.Part
	newcell.Verts = o.Verts
	newcell.FTags = o.FTags
	newcell.STags = o.STags
	newcell.JlinId = o.JlinId
	newcell.JsldId = o.JsldId

	// neighbours
	newcell.Neighs = o.Neighs

	// new cell type
	ctype := o.Shp.Type
	if lbb {
		ctype = o.Shp.BasicType
	}

	// derived
	if o.Shp.Nurbs == nil {
		newcell.Shp = shp.Get(ctype, o.GoroutineId)
		if newcell.Shp == nil {
			chk.Panic("cannot allocate \"shape\" structure for cell type = %q\n", ctype)
		}
	} else {
		chk.Panic("cannot handle similar cells with NURBS yet")
	}
	newcell.FaceBcs = o.FaceBcs
	newcell.GoroutineId = o.GoroutineId

	// specific problems data
	newcell.IsJoint = o.IsJoint
	newcell.SeepVerts = o.SeepVerts

	// NURBS
	newcell.Nrb = o.Nrb
	newcell.Span = o.Span
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

// Draw2d draws 2D mesh
func (o *Mesh) Draw2d() {

	// auxiliary
	type triple struct{ a, b, c int }   // points on edge
	edgesdrawn := make(map[triple]bool) // edges drawn already
	var tri triple

	// loop over cells
	for _, cell := range o.Cells {

		// loop edges of cells
		for _, lvids := range cell.Shp.FaceLocalVerts {

			// set triple of nodes
			tri.a = cell.Verts[lvids[0]]
			tri.b = cell.Verts[lvids[1]]
			nv := len(lvids)
			if nv > 2 {
				tri.c = cell.Verts[lvids[2]]
			} else {
				tri.c = len(o.Verts) + 1 // indicator of not-available
			}
			utl.IntSort3(&tri.a, &tri.b, &tri.c)

			// draw edge if not drawn yet
			if _, drawn := edgesdrawn[tri]; !drawn {
				x := make([]float64, nv)
				y := make([]float64, nv)
				x[0] = o.Verts[tri.a].C[0]
				y[0] = o.Verts[tri.a].C[1]
				if nv == 3 {
					x[1] = o.Verts[tri.c].C[0]
					y[1] = o.Verts[tri.c].C[1]
					x[2] = o.Verts[tri.b].C[0]
					y[2] = o.Verts[tri.b].C[1]
				} else {
					x[1] = o.Verts[tri.b].C[0]
					y[1] = o.Verts[tri.b].C[1]
				}
				plt.Plot(x, y, "'k-o', ms=3, clip_on=0")
				edgesdrawn[tri] = true
			}
		}

		// add middle node
		if cell.Type == "qua9" {
			vid := cell.Verts[8]
			x := o.Verts[vid].C[0]
			y := o.Verts[vid].C[1]
			plt.PlotOne(x, y, "'ko', ms=3, clip_on=0")
		}

		// linear cells
		if strings.HasPrefix(cell.Type, "lin") {
			nv := len(cell.Verts)
			x := make([]float64, nv)
			y := make([]float64, nv)
			for i, vid := range cell.Verts {
				x[i] = o.Verts[vid].C[0]
				y[i] = o.Verts[vid].C[1]
			}
			plt.Plot(x, y, "'-o', ms=3, clip_on=0, color='#41045a', lw=2")
		}
	}

	// set up
	plt.Equal()
	plt.AxisRange(o.Xmin, o.Xmax, o.Ymin, o.Ymax)
	plt.AxisOff()
}

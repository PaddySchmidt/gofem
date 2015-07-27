// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"sort"

	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// FaceCond holds information of one single face boundary condition. Example:
//
//                      -12 => "qn", "seepH"
//   36    -12     35   -11 => "ux", "seepH"
//    (3)--------(2)
//     |    2     |     face id => conditions
//     |          |           0 => <nil>
//     |3        1| -11       1 => {"ux", "seepH"} => localVerts={1,2} => globalVerts={34,35}
//     |          |           2 => {"qn", "seepH"} => localVerts={2,3} => globalVerts={35,36}
//     |    0     |           3 => <nil>
//    (0)--------(1)
//   33            34    "seepH" => localVerts={1,2,3}
//                       "seepH" => globalVerts={34,35,36}
//
//   face id => condition
//         1 => "ux"    => localVerts={1,2} => globalVerts={34,35}
//         1 => "seepH" => localVerts={1,2} => globalVerts={34,35}
//         2 => "qn"    => localVerts={2,3} => globalVerts={35,36}
//         2 => "seepH" => localVerts={2,3} => globalVerts={35,36}
type FaceCond struct {
	FaceId      int      // msh: cell's face local id
	LocalVerts  []int    // msh: cell's face local vertices ids (sorted)
	GlobalVerts []int    // msh: global vertices ids (sorted)
	Cond        string   // sim: condition; e.g. "qn" or "seepH"
	Func        fun.Func // sim: function to compute boundary condition
	Extra       string   // sim: extra information
}

// FaceConds hold many face boundary conditions
type FaceConds []*FaceCond

// GetVerts gets all vertices with any of the given conditions
//  Example: "seepH" => localVerts={1,2,3}
func (o FaceConds) GetVerts(conds ...string) (verts []int) {

	// for each face condition
	for _, fc := range o {

		// check if fc has any of conds
		if utl.StrIndexSmall(conds, fc.Cond) < 0 {
			continue
		}

		// add local verts
		for _, lv := range fc.LocalVerts {

			// check if vert was added already
			if utl.IntIndexSmall(verts, lv) < 0 {
				verts = append(verts, lv) // add a vert
			}
		}
	}

	// results
	sort.Ints(verts)
	return
}

// SetFaceConds sets face boundary conditions map in cell
func (o *Cell) SetFaceConds(stg *Stage, functions FuncsData) (err error) {

	// for each face tag
	o.FaceBcs = make([]*FaceCond, 0)
	for faceId, faceTag := range o.FTags {

		// skip zero or positive tags
		if faceTag >= 0 {
			continue
		}

		// find face boundary condition and skip nil data
		faceBc := stg.GetFaceBc(faceTag)
		if faceBc == nil {
			continue
		}

		// local and global ids of vertices on face
		lverts := shp.GetFaceLocalVerts(o.Type, faceId)
		gverts := make([]int, len(lverts))
		for i, l := range lverts {
			gverts[i] = o.Verts[l]
		}

		// for each boundary key such as "ux", "qn", etc.
		for j, key := range faceBc.Keys {
			fcn := functions.Get(faceBc.Funcs[j])
			if fcn == nil {
				return chk.Err("cannot find function named %q corresponding to face tag %d (@ cell %d)", faceBc.Funcs[j], faceTag, o.Id)
			}
			o.FaceBcs = append(o.FaceBcs, &FaceCond{faceId, lverts, gverts, key, fcn, faceBc.Extra})
		}
	}
	return
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
)

func compute_extrapolated_values() {

	// auxiliary
	verts := Dom.Msh.Verts
	cells := Dom.Msh.Cells

	// allocate structures for extrapolation
	nverts := len(verts)
	ExVals = make([]map[string]float64, nverts)
	counts := make([]map[string]float64, nverts)
	for i := 0; i < nverts; i++ {
		ExVals[i] = make(map[string]float64)
		counts[i] = make(map[string]float64)
	}

	// loop over elements
	for _, ele := range Dom.Elems {

		// get shape and integration points from known elements
		var sha *shp.Shape
		var ips []*shp.Ipoint
		switch e := ele.(type) {
		case *fem.ElemP:
			sha = e.Shp
			ips = e.IpsElem
		case *fem.ElemU:
			sha = e.Shp
			ips = e.IpsElem
		case *fem.ElemUP:
			sha = e.U.Shp
			ips = e.U.IpsElem
		}
		if sha == nil {
			chk.Panic("cannot get shape structrue from element")
		}

		// compute Extrapolator matrix
		Emat := la.MatAlloc(sha.Nverts, len(ips))
		err := sha.Extrapolator(Emat, ips)
		if err != nil {
			chk.Panic("cannot compute extrapolator matrix: %v", err)
		}

		// get ips data
		dat := ele.OutIpsData()

		// perform extrapolation
		cell := cells[ele.Id()]
		for j := 0; j < len(ips); j++ {
			vals := dat[j].Calc(Dom.Sol)
			for _, key := range Extrap {
				if val, ok := vals[key]; ok {
					for i := 0; i < sha.Nverts; i++ {
						v := cell.Verts[i]
						ExVals[v][key] += Emat[i][j] * val
					}
				} else {
					chk.Panic("ip does not have key = %s", key)
				}
			}
		}

		// increment counter
		for i := 0; i < sha.Nverts; i++ {
			v := cell.Verts[i]
			for _, key := range Extrap {
				counts[v][key] += 1
			}
		}
	}

	// compute average
	for i := 0; i < nverts; i++ {
		for key, cnt := range counts[i] {
			ExVals[i][key] /= cnt
		}
	}
}

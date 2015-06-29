// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"sort"

	"github.com/cpmech/gosl/chk"
)

// PointLocator defines interface for locating space positions
type Locator interface {
	Locate() Points
}

// At implements locator at point => PointLocator
type At []float64

// AtIp implements locator at integration point => PointLocator
//  Note: this is useful when there are vertices overlapping ips
type AtIp []float64

// N implements node locator
// Ids or tags of vertices can be stored in Verts
type N []int

// P implements [element][integrationPoint] locator
// Pairs of ids or tags of cells and integration points indices can be stored in Cells
type P [][]int

// Along implements locator along line => LineLocator
//  Example: with 2 points in 3D: {{0,0,0}, {1,1,1}}
type Along [][]float64

// AlongX implements LineLocator with []float64{y_cte} or []float64{y_cte, z_cte}
type AlongX []float64

// AlongY implements LineLocator with []float64{x_cte} or []float64{x_cte, z_cte}
type AlongY []float64

// AlongZ implements LineLocator with []float64{x_cte, y_cte}
type AlongZ []float64

// Locate finds points
func (o At) Locate() Points {

	// node
	vid := NodBins.Find(o)
	if vid >= 0 {
		q := get_nod_point(vid, nil)
		if q != nil {
			return Points{q}
		}
	}

	// integration point
	ipid := IpsBins.Find(o)
	if ipid >= 0 {
		q := get_ip_point(ipid, nil)
		if q != nil {
			return Points{q}
		}
	}
	return nil
}

// Locate finds integration points
func (o AtIp) Locate() Points {
	ipid := IpsBins.Find(o)
	if ipid >= 0 {
		q := get_ip_point(ipid, nil)
		if q != nil {
			return Points{q}
		}
	}
	return nil
}

// Locate finds nodes
func (o N) Locate() (res Points) {
	var A []float64 // reference point
	for _, idortag := range o {
		if idortag < 0 {
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				q := get_nod_point(v.Id, A)
				if q != nil {
					res = append(res, q)
					if A == nil {
						A = q.X
					}
				}
			}
		} else {
			vid := idortag
			q := get_nod_point(vid, A)
			if q != nil {
				res = append(res, q)
				if A == nil {
					A = q.X
				}
			}
		}
	}
	if len(res) < len(o) {
		chk.Panic("cannot locate all nodes in %v", o)
	}
	return
}

// Locate finds points
func (o P) Locate() (res Points) {
	var A []float64 // reference point
	append_to_res := func(cid, idx int) {
		ipid := Cid2ips[cid][idx]
		q := get_ip_point(ipid, A)
		if q != nil {
			res = append(res, q)
			if A == nil {
				A = q.X
			}
		}
	}
	ncells := len(o)
	for i := 0; i < ncells; i++ {
		if len(o[i]) != 2 {
			continue
		}
		idortag := o[i][0]
		if idortag < 0 {
			tag := idortag
			cells := Dom.Msh.CellTag2cells[tag]
			for _, c := range cells {
				cid := c.Id
				idx := o[i][1]
				if idx >= len(Cid2ips[cid]) {
					continue
				}
				append_to_res(cid, idx)
			}
		} else {
			cid := idortag
			idx := o[i][1]
			if idx >= len(Cid2ips[cid]) {
				continue
			}
			append_to_res(cid, idx)
		}
	}
	if len(res) < len(o) {
		chk.Panic("cannot locate all points in %v", o)
	}
	return
}

// Locate finds points
func (o Along) Locate() (res Points) {

	// check if there are two points
	if len(o) != 2 {
		return
	}
	A := o[0]
	B := o[1]

	// node quantities
	ids := NodBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		q := get_nod_point(id, A)
		if q != nil {
			res = append(res, q)
		}
	}

	// integration point quantitites
	ids = IpsBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		q := get_ip_point(id, A)
		if q != nil {
			res = append(res, q)
		}
	}
	sort.Sort(res)
	return
}

// Locate finds points
func (o AlongX) Locate() (res Points) {
	y_cte, z_cte := o[0], 0.0
	if len(o) > 1 {
		z_cte = o[1]
	}
	return Along{{0, y_cte, z_cte}, {1, y_cte, z_cte}}.Locate()
}

// Locate finds points
func (o AlongY) Locate() (res Points) {
	x_cte, z_cte := o[0], 0.0
	if len(o) > 1 {
		z_cte = o[1]
	}
	return Along{{x_cte, 0, z_cte}, {x_cte, 1, z_cte}}.Locate()
}

// Locate finds points
func (o AlongZ) Locate() (res Points) {
	x_cte, y_cte := o[0], o[1]
	return Along{{x_cte, y_cte, 0}, {x_cte, y_cte, 1}}.Locate()
}

// AllIps returns all cell/ip indices
func AllIps() P {
	var p [][]int
	for i, ips := range Cid2ips {
		for j, _ := range ips {
			p = append(p, []int{i, j})
		}
	}
	return p
}

func AllNodes() N {
	var res []int
	for _, nod := range Dom.Nodes {
		res = append(res, nod.Vert.Id)
	}
	return res
}

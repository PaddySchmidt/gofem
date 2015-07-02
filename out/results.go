// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"strings"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/utl"
)

// Define defines aliases
//  alias -- an alias to a group of points, an individual point, or to a set of points.
//           Example: "A", "left-column" or "a b c". If the number of points found is different
//           than the number of aliases, a group is created.
//  Note:
//    To use spaces in aliases, prefix the alias with an exclamation mark; e.g "!right column"
func Define(alias string, loc Locator) {

	// check
	if len(alias) < 1 {
		chk.Panic("alias must have at least one character. %q is invalid", alias)
	}

	// set planes map
	if plane, ok := loc.(*PlaneData); ok {
		Planes[alias] = plane
		defer func() {
			plane.ConnectResults(alias)
		}()
	}

	// locate points
	pts := loc.Locate()
	if len(pts) < 1 {
		chk.Panic("cannot define entities with alias=%q and locator=%v", alias, loc)
	}

	// set results map
	if alias[0] == '!' {
		Results[alias[1:]] = pts
		return
	}
	lbls := strings.Fields(alias)
	if len(lbls) == len(pts) {
		for i, l := range lbls {
			Results[l] = []*Point{pts[i]}
		}
		return
	}
	Results[alias] = pts
}

// LoadResults loads all results after points are defined
//  times -- specified selected output times
//           use nil to indicate that all times are required
func LoadResults(times []float64) {

	// selected output times and indices
	if times == nil {
		times = Sum.OutTimes
	}
	TimeInds, Times = utl.GetITout(Sum.OutTimes, times, TolT)

	// for each selected output time
	for _, tidx := range TimeInds {

		// input results into domain
		if !Dom.In(Sum, tidx, true) {
			chk.Panic("cannot load results into domain; please check log file")
		}
		if Dom.Ny != len(Dom.Sol.Y) {
			chk.Panic("inconsistency of results detected: summary and simulation file might be different")
		}

		// extrapolation
		if Extrap != nil {
			ComputeExtrapolatedValues(Extrap)
		}

		// for each point
		for _, pts := range Results {
			for _, p := range pts {

				// node or integration point id
				vid := p.Vid
				pid := p.IpId

				// handle node
				if vid >= 0 {

					// add dofs to results map
					nod := Dom.Vid2node[vid]
					for _, dof := range nod.Dofs {
						if dof != nil {
							utl.StrDblsMapAppend(&p.Vals, dof.Key, Dom.Sol.Y[dof.Eq])
						}
					}

					// add extrapolated values to results map
					if ExVals != nil {
						for key, val := range ExVals[vid] {
							utl.StrDblsMapAppend(&p.Vals, key, val)
						}
					}
				}

				// handle integration point
				if pid >= 0 {
					dat := Ipoints[pid]
					vals := dat.Calc(Dom.Sol)
					for key, val := range vals {
						utl.StrDblsMapAppend(&p.Vals, key, val)
					}
				}
			}
		}
	}
}

// GetRes gets results as a time or space series corresponding to a given alias
// for a single point or set of points.
//  idxI -- index in TimeInds slice corresponding to selected output time; use -1 for the last item.
//          If alias defines a single point, the whole time series is returned and idxI is ignored.
func GetRes(key, alias string, idxI int) []float64 {
	if idxI < 0 {
		idxI = len(TimeInds) - 1
	}
	if pts, ok := Results[alias]; ok {
		if len(pts) == 1 {
			for k, v := range pts[0].Vals {
				if k == key {
					return v
				}
			}
		} else {
			var res []float64
			for _, p := range pts {
				for k, v := range p.Vals {
					if k == key {
						res = append(res, v[idxI])
					}
				}
			}
			return res
		}
	}
	chk.Panic("cannot get %q at %q", key, alias)
	return nil
}

// GetIds return the ids corresponding to alias
func GetIds(alias string) (vids, ipids []int) {
	if pts, ok := Results[alias]; ok {
		for _, p := range pts {
			if p.Vid >= 0 {
				vids = append(vids, p.Vid)
			}
			if p.IpId >= 0 {
				ipids = append(ipids, p.IpId)
			}
		}
	}
	return
}

// GetCoords returns the coordinates of a single point
func GetCoords(alias string) []float64 {
	if pts, ok := Results[alias]; ok {
		if len(pts) == 1 {
			return pts[0].X
		}
	}
	chk.Panic("cannot get coordinates of point with alias %q (make sure this alias corresponds to a single point)", alias)
	return nil
}

// GetDist returns the distance from a reference point on the given line with selected points
// if they contain a given key
//  key -- use any to get coordinates of points with any key such as "ux", "pl", etc.
func GetDist(key, alias string) (dist []float64) {
	any := key == "any"
	if pts, ok := Results[alias]; ok {
		for _, p := range pts {
			for k, _ := range p.Vals {
				if k == key || any {
					dist = append(dist, p.Dist)
					break
				}
			}
		}
		return
	}
	chk.Panic("cannot get distance with key %q and alias %q", key, alias)
	return
}

// GetXYZ returns the x-y-z coordinates of selected points that have a specified key
//  key -- use any to get coordinates of points with any key such as "ux", "pl", etc.
func GetXYZ(key, alias string) (x, y, z []float64) {
	any := key == "any"
	if pts, ok := Results[alias]; ok {
		for _, p := range pts {
			for k, _ := range p.Vals {
				if k == key || any {
					x = append(x, p.X[0])
					y = append(y, p.X[1])
					if len(p.X) == 3 {
						z = append(z, p.X[2])
					}
					break
				}
			}
		}
		return
	}
	chk.Panic("cannot get x-y-z coordinates with key %q and alias %q", key, alias)
	return
}

// Integrate integrates key along direction "x", "y", or "z"
//  idxI -- index in TimeInds slice corresponding to selected output time; use -1 for the last item.
//          If alias defines a single point, the whole time series is returned and idxI is ignored.
func Integrate(key, alias, along string, idxI int) float64 {
	if idxI < 0 {
		idxI = len(TimeInds) - 1
	}
	y := GetRes(key, alias, idxI)
	var x []float64
	switch along {
	case "x":
		x, _, _ = GetXYZ(key, alias)
	case "y":
		_, x, _ = GetXYZ(key, alias)
	case "z":
		_, _, x = GetXYZ(key, alias)
	}
	var err error
	_, x, y, _, err = utl.SortQuadruples(nil, x, y, nil, "x")
	if err != nil {
		chk.Panic("%q: cannot integrate %q along %q:\n%v\n", alias, key, along, err)
	}
	return num.Trapz(x, y)
}

// IntegOnPlane integrates the results of nodes located on a plane perpedicular to either {x,y,z}
// The mesh on the plane must be structured and with rectangle elements in order to build
// an {u,v} coordinates system
// Input:
//   key   -- node variable; e.g. "uz", "ex_sz"
//   alias -- alias of plane; e.g. "surf"
// Output:
//   res -- the result of the integration res = ∫_u ∫_v f(u,v) du dv for all output times
func IntegOnPlane(key, alias string) (res []float64) {

	// get plane
	plane, ok := Planes[alias]
	if !ok {
		chk.Panic("cannot get plane with alias=%q", alias)
	}

	// compute integral
	res = make([]float64, len(TimeInds))
	u := make([]float64, 2)
	for idxI, k := range TimeInds {
		for j := 0; j < plane.Nu[1]; j++ {
			u[1] = plane.Umin[1] + float64(j)*plane.Du[1]
			for i := 0; i < plane.Nu[0]; i++ {
				u[0] = plane.Umin[0] + float64(i)*plane.Du[0]
				id := plane.Ubins.Find(u)
				if id < 0 {
					chk.Panic("IntegOnPlane with key=%q and alias=%q\nError: cannot find {u,v} coordinate in any bin. u=%v vid/ipid=%d", key, alias, u, id)
				}
				vals, ok := plane.id2pt[id].Vals[key]
				if !ok {
					chk.Panic("cannot find results with key=%q for plane with alias=%q", key, alias)
				}
				plane.F[i][j] = vals[idxI]
			}
		}
		res[k] = num.Simps2D(plane.Du[0], plane.Du[1], plane.F)
	}
	return
}

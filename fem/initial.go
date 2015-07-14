// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import "github.com/cpmech/gofem/inp"

// SetInitial sets the initial state
func (o *Domain) SetInitial(stg *inp.Stage) (ok bool) {

	// get function
	fcn := Global.Sim.Functions.Get(stg.Initial.Fcn)
	if LogErrCond(fcn == nil, "cannot get function named %q", stg.Initial.Fcn) {
		return
	}

	// set nodes
	for _, nod := range o.Nodes {
		eq := nod.GetEq("h")
		if LogErrCond(eq < 0, "SetInitial: equation cannot be found for dof='phi' in node=%d", nod.Vert.Id) {
			return
		}
		o.Sol.Y[eq] = fcn.F(0, nod.Vert.C)
	}

	// success
	return true
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import "github.com/cpmech/gofem/inp"

// SetInitial sets the initial state
func (o *Domain) SetInitial(stg *inp.Stage) (ok bool) {

	// check
	if LogErrCond(len(stg.Initial.Fcns) != len(stg.Initial.Dofs), "number of functions (fcns) must be equal to number of dofs for setting initial values. %d != %d", len(stg.Initial.Fcns), len(stg.Initial.Dofs)) {
		return
	}

	// loop over functions
	for i, fname := range stg.Initial.Fcns {

		// get function
		fcn := Global.Sim.Functions.Get(fname)
		if LogErrCond(fcn == nil, "cannot get function named %q", fname) {
			return
		}

		// set nodes
		key := stg.Initial.Dofs[i]
		for _, nod := range o.Nodes {
			eq := nod.GetEq(key)
			if LogErrCond(eq < 0, "dof=%q cannot be found in node=%d for setting initial values", key, nod.Vert.Id) {
				return
			}
			o.Sol.Y[eq] = fcn.F(0, nod.Vert.C)
		}
	}

	// success
	return true
}

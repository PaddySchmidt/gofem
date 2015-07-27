// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
)

// SetInitial sets the initial state
func (o *Domain) SetInitial(stg *inp.Stage) (err error) {

	// check
	if len(stg.Initial.Fcns) != len(stg.Initial.Dofs) {
		return chk.Err("number of functions (fcns) must be equal to number of dofs for setting initial values. %d != %d", len(stg.Initial.Fcns), len(stg.Initial.Dofs))
	}

	// loop over functions
	for i, fname := range stg.Initial.Fcns {

		// get function
		fcn := o.Sim.Functions.Get(fname)
		if fcn == nil {
			return chk.Err("cannot get function named %q", fname)
		}

		// set nodes
		key := stg.Initial.Dofs[i]
		for _, nod := range o.Nodes {
			eq := nod.GetEq(key)
			if eq < 0 {
				return chk.Err("dof=%q cannot be found in node=%d for setting initial values", key, nod.Vert.Id)
			}
			o.Sol.Y[eq] = fcn.F(0, nod.Vert.C)
		}
	}
	return
}

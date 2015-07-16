// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func init() {
	io.Verbose = false
}

func verbose() {
	io.Verbose = true
	chk.Verbose = true
}

func get_nids_eqs(dom *Domain) (nids, eqs []int) {
	for _, nod := range dom.Nodes {
		nids = append(nids, nod.Vert.Id)
		for _, dof := range nod.Dofs {
			eqs = append(eqs, dof.Eq)
		}
	}
	return
}

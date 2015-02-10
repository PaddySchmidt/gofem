// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/utl"
)

// Dof holds information about a degree-of-freedom == solution variable
type Dof struct {
	Ukey string // primary variable key. e.g. "ux"
	Eq   int    // equation number
}

// String returns the string representation of this Dof
func (o *Dof) String() string {
	l := utl.Sf("{ \"Ukey\" : %s  \"Eq\" : %d } ", o.Ukey, o.Eq)
	return l
}

// Node holds node dofs information
type Node struct {
	dofs []*Dof    // degrees-of-freedom == solution variables
	vert *inp.Vert // pointer to Vertex
}

// NewNode allocates a new Node
func NewNode(v *inp.Vert) *Node {
	return &Node{[]*Dof{}, v}
}

// String returns the string representation of this node
func (o *Node) String() string {
	l := "{ "
	l += utl.Sf(" \"Id\" :  %d ", o.vert.Id)
	for _, dof := range o.dofs {
		l += dof.String()
	}
	l += " } "
	return l
}

// AddDof adds a new dof to thisnode; ignores it if it exists already
//  nexteq -- is the next equation number == eqnum + 1;
//            returns eqnum if dof exists already
func (o *Node) AddDofAndEq(ukey string, eqnum int) (nexteq int) {

	// check if ukey exists already
	for _, dof := range o.dofs {
		if ukey == dof.Ukey {
			return eqnum
		}
	}

	// add new Dof
	o.dofs = append(o.dofs, &Dof{ukey, eqnum})
	return eqnum + 1
}

// SetEq numbers a specific Dof with the equation number in the current (stage) global system
func (o *Node) SetEq(ukey string, eqNumber int) {
	// set Dof
}

// GetDof returns the Dof structure for given Dof name (ukey)
//  Note: returns nil if not found
func (o *Node) GetDof(ukey string) *Dof {
	for _, dof := range o.dofs {
		if dof.Ukey == ukey {
			return dof
		}
	}
	return nil
}

// GetDofOrPanic returns the Dof structure for given Dof name (ukey)
func (o *Node) GetDofOrPanic(ukey string) *Dof {
	for _, dof := range o.dofs {
		if dof.Ukey == ukey {
			return dof
		}
	}
	PanicOrNot(true, "cannot find dof named %q", ukey)
	return nil
}

// GetEq returns the equation number for given Dof name (ukey)
//  Note: returns -1 if not found
func (o *Node) GetEq(ukey string) (eqNumber int) {
	for _, dof := range o.dofs {
		if dof.Ukey == ukey {
			return dof.Eq
		}
	}
	return -1
}
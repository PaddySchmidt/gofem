// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// xfem_set_info sets extra information for XFEM elements
func xfem_set_info(info *Info, cell *inp.Cell, edat *inp.ElemData) (ykeys []string) {

	// flags
	xmat, xcrk, xfem := false, false, false
	if s_xmat, found := io.Keycode(edat.Extra, "xmat"); found {
		xmat = io.Atob(s_xmat)
		xfem = true
	}
	if s_xcrk, found := io.Keycode(edat.Extra, "xcrk"); found {
		xcrk = io.Atob(s_xcrk)
		xfem = true
	}
	_ = xcrk // TODO

	// skip if not XFEM
	if !xfem {
		return
	}

	// extra information
	nverts := cell.Shp.Nverts
	if xmat {
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = append(info.Dofs[m], "am")
		}
		info.Y2F["am"] = "nil"
		info.T2vars = append(info.T2vars, "am")
	}
	return
}

// xfem_init initialises variables need by xfem model
func (o *ElemU) xfem_init(edat *inp.ElemData) {

	// flags
	o.Xmat, o.Xcrk, o.Xfem = false, false, false
	if s_xmat, found := io.Keycode(edat.Extra, "xmat"); found {
		o.Xmat = io.Atob(s_xmat)
		o.Xfem = true
	}
	if s_xcrk, found := io.Keycode(edat.Extra, "xcrk"); found {
		o.Xcrk = io.Atob(s_xcrk)
		o.Xfem = true
	}

	// skip if not XFEM
	if !o.Xfem {
		return
	}

	// auxiliary variables
	switch {
	case o.Xmat:
		o.Na = 1
	case o.Xcrk:
		o.Na = 4 // TODO: this is 2D only
		if o.Ndim != 2 {
			chk.Panic("xfem with crack works in 2D only for now")
		}
	}

	// allocate extrem XFEM coupling matrices
	nverts := o.Cell.Shp.Nverts
	o.Amap = make([]int, o.Na*nverts)
	o.Kua = la.MatAlloc(o.Nu, o.Na)
	o.Kau = la.MatAlloc(o.Na, o.Nu)
	o.Kaa = la.MatAlloc(o.Na, o.Na)
}

// xfem_add_to_rhs adds contribution to rhs due to xfem model
func (o *ElemU) xfem_add_to_rhs(fb []float64, sol *Solution) (err error) {
	return
}

// contact_add_to_jac adds coupled equations due to xfem to Jacobian
func (o *ElemU) xfem_add_to_jac(Kb *la.Triplet, sol *Solution) (err error) {
	return
}

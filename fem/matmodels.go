// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/msolid"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// GetAndInitPorousModel get porous model from material name
// It returns nil on errors, after logging
func GetAndInitPorousModel(mdb *inp.MatDb, matname, simfnk string) (mdl *mporous.Model, err error) {

	// materials
	cndmat, lrmmat, pormat, err := mdb.GroupGet3(matname, "c", "l", "p")
	if err != nil {
		err = chk.Err("materials database failed on getting %q (porous) group:\n%v", matname, err)
		return
	}

	// conductivity models
	getnew := false
	cnd := mconduct.GetModel(simfnk, cndmat.Name, cndmat.Model, getnew)
	if cnd == nil {
		err = chk.Err("cannot allocate conductivity models with name=%q", cndmat.Model)
		return
	}

	// retention model
	lrm := mreten.GetModel(simfnk, lrmmat.Name, lrmmat.Model, getnew)
	if lrm == nil {
		err = chk.Err("cannot allocate liquid retention model with name=%q", lrmmat.Model)
		return
	}

	// porous model
	mdl = mporous.GetModel(simfnk, pormat.Name, getnew)
	if mdl == nil {
		err = chk.Err("cannot allocate model for porous medium with name=%q", pormat.Name)
		return
	}

	// initialise all models
	// TODO: initialise just once
	err = cnd.Init(cndmat.Prms)
	if err != nil {
		err = chk.Err("cannot initialise conductivity model:\n%v", err)
		return
	}
	err = lrm.Init(lrmmat.Prms)
	if err != nil {
		err = chk.Err("cannot initialise liquid retention model:\n%v", err)
		return
	}
	err = mdl.Init(pormat.Prms, cnd, lrm)
	if err != nil {
		err = chk.Err("cannot initialise porous model:\n%v", err)
		return
	}
	return
}

func GetAndInitSolidModel(mdb *inp.MatDb, matname, simfnk string, ndim int, pstress bool) (mdl msolid.Model, prms fun.Prms, err error) {

	// material name
	matdata := mdb.Get(matname)
	if matdata == nil {
		err = chk.Err("materials database failed on getting %q (solid) material\n", matname)
		return
	}
	mdlname := matdata.Model

	// handle groups
	if mdlname == "group" {
		if s_matname, found := io.Keycode(matdata.Extra, "s"); found {
			matname = s_matname
			matdata = mdb.Get(matname)
			if matdata == nil {
				err = chk.Err("materials database failed on getting %q (solid/sub) material\n", matname)
				return
			}
			mdlname = matdata.Model
		} else {
			err = chk.Err("cannot find solid model in grouped material data. 's' subkey needed in Extra field")
			return
		}
	}

	// initialise model
	mdl, existent := msolid.GetModel(simfnk, matname, mdlname, false)
	if mdl == nil {
		err = chk.Err("cannot find solid model named %q", mdlname)
		return
	}
	if !existent {
		err = mdl.Init(ndim, pstress, matdata.Prms)
		if err != nil {
			err = chk.Err("solid model initialisation failed:\n%v", err)
			return
		}
	}

	// results
	prms = matdata.Prms
	return
}

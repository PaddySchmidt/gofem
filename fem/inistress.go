// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

// SetIniStress sets the initial state with initial stresses
func (o *Domain) SetIniStress(stg *inp.Stage) (err error) {

	// set elements with homogeneous stress state
	dat := stg.IniStress
	if dat.Hom {

		// isotropic state
		if dat.Iso {
			for _, e := range o.ElemIntvars {

				// build map with isotropic and homogeneus state
				coords := e.Ipoints()
				nip := len(coords)
				v := utl.DblVals(nip, dat.S0)
				ivs := map[string][]float64{"sx": v, "sy": v, "sz": v}

				// set element's states
				err = e.SetIniIvs(o.Sol, ivs)
				if err != nil {
					return chk.Err("homogeneous/isotropic: element's internal values setting failed:\n%v", err)
				}
			}
			return
		}

		// plane-strain state
		if dat.Psa {
			sz := dat.Nu * (dat.Sh + dat.Sv)
			for _, e := range o.ElemIntvars {

				// build map with plane-strain and homogeneus state
				coords := e.Ipoints()
				nip := len(coords)
				vx := utl.DblVals(nip, dat.Sh)
				vy := utl.DblVals(nip, dat.Sv)
				vz := utl.DblVals(nip, sz)
				ivs := map[string][]float64{"sx": vx, "sy": vy, "sz": vz}

				// set element's states
				err = e.SetIniIvs(o.Sol, ivs)
				if err != nil {
					return chk.Err("homogeneous/plane-strain: element's internal values setting failed:\n%v", err)
				}
			}
			return
		}
	}
	return
}

// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// contact_set_info sets extra information for elements with contact
func contact_set_info(info *Info, cell *inp.Cell, edat *inp.ElemData) {
	if len(cell.FaceBcs) > 0 { // vertices on faces with contact
		lverts := cell.FaceBcs.GetVerts("contact")
		for _, m := range lverts {
			info.Dofs[m] = append(info.Dofs[m], "qb")
		}
		if len(lverts) > 0 {
			info.Y2F["qb"] = "nil"
			info.T2vars = append(info.T2vars, "qb")
		}
	}
	return
}

// contact_init initialises variables need by contact model
func (o *ElemU) contact_init(edat *inp.ElemData) {

	// vertices on faces with contact
	var contactverts []int
	if len(o.Cell.FaceBcs) > 0 {
		lverts := o.Cell.FaceBcs.GetVerts("contact")
		for _, m := range lverts {
			contactverts = append(contactverts, m)
		}
	}
	o.Nq = len(contactverts)

	// contact flag
	o.HasContact = o.Nq > 0
	if !o.HasContact {
		return
	}

	// vertices on contact face; numbering
	o.ContactId2vid = contactverts
	o.Vid2contactId = utl.IntVals(o.Nu, -1)
	o.Qmap = make([]int, o.Nq)
	for μ, m := range o.ContactId2vid {
		o.Vid2contactId[m] = μ
	}

	// flags
	o.Macaulay, o.βrmp, o.κ = GetContactFaceFlags(edat.Extra)

	// allocate coupling matrices
	o.Kuq = la.MatAlloc(o.Nu, o.Nq)
	o.Kqu = la.MatAlloc(o.Nq, o.Nu)
	o.Kqq = la.MatAlloc(o.Nq, o.Nq)
}

// contact_add_to_rhs adds contribution to rhs due to contact modelling
func (o *ElemU) contact_add_to_rhs(fb []float64, sol *Solution) (err error) {

	// compute surface integral
	var qb, db, rmp, rx, rq float64
	for _, nbc := range o.NatBcs {

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Cell.Shp.Sf
			nvec := o.Cell.Shp.Fnvec

			// select natural boundary condition type
			switch nbc.Key {

			// contact
			case "contact":

				// variables extrapolated to face
				qb = o.fipvars(iface, sol)
				xf := o.Cell.Shp.FaceIpRealCoords(o.X, ipf, iface)
				la.VecAdd(xf, 1, o.us) // add displacement: x = X + u
				db = o.contact_g(xf)

				// compute residuals
				coef := ipf[3] * o.Thickness
				Jf := la.VecNorm(nvec)
				rmp = o.contact_ramp(qb + o.κ*db)
				rx = rmp
				rq = qb - rmp
				for j, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2contactId[m]
					fb[o.Qmap[μ]] -= coef * Sf[j] * rq * Jf // -residual
					for i := 0; i < o.Ndim; i++ {
						r := o.Umap[i+m*o.Ndim]
						fb[r] -= coef * Sf[j] * rx * nvec[i] // -extra term
					}
				}
			}
		}
	}
	return
}

// contact_add_to_jac adds coupled equations due to contact modelling to Jacobian
func (o *ElemU) contact_add_to_jac(Kb *la.Triplet, sol *Solution) (err error) {

	// clear matrices
	for i := 0; i < o.Nq; i++ {
		for j := 0; j < o.Nq; j++ {
			o.Kqq[i][j] = 0
		}
		for j := 0; j < o.Nu; j++ {
			o.Kqu[i][j] = 0
			o.Kuq[j][i] = 0
		}
	}

	// compute surface integral
	var qb, db, Hb float64
	dddu := make([]float64, o.Ndim)
	for _, nbc := range o.NatBcs {

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// contact
			switch nbc.Key {
			case "contact":

				// interpolation functions and gradients @ face
				iface := nbc.IdxFace
				err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
				if err != nil {
					return
				}
				Sf := o.Cell.Shp.Sf
				nvec := o.Cell.Shp.Fnvec
				coef := ipf[3] * o.Thickness
				Jf := la.VecNorm(nvec)

				// variables extrapolated to face
				qb = o.fipvars(iface, sol)
				xf := o.Cell.Shp.FaceIpRealCoords(o.X, ipf, iface)
				la.VecAdd(xf, 1, o.us) // add displacement: x = X + u
				db = o.contact_g(xf)
				o.contact_dgdx(dddu, xf)

				// compute derivatives
				Hb = o.contact_rampD1(qb + o.κ*db)
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2contactId[m]
					for j, n := range o.Cell.Shp.FaceLocalVerts[iface] {
						ν := o.Vid2contactId[n]
						o.Kqq[μ][ν] += coef * Jf * Sf[i] * Sf[j] * (1.0 - Hb)
						for k := 0; k < o.Ndim; k++ {
							r := k + m*o.Ndim
							c := k + n*o.Ndim
							val := coef * Sf[i] * Sf[j] * Hb * o.κ * dddu[k] * Jf
							o.Kuq[r][ν] += coef * Sf[i] * Sf[j] * Hb * nvec[k]
							o.Kqu[μ][c] -= val
							o.K[r][c] += val
						}
					}
				}
			}
		}
	}

	// add Ks to sparse matrix Kb
	for i, I := range o.Qmap {
		for j, J := range o.Qmap {
			Kb.Put(I, J, o.Kqq[i][j])
		}
		for j, J := range o.Umap {
			Kb.Put(I, J, o.Kqu[i][j])
			Kb.Put(J, I, o.Kuq[j][i])
		}
	}
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// contact_ramp implements the ramp function
func (o *ElemU) contact_ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.βrmp)
}

// contact_rampderiv returns the ramp function first derivative
func (o *ElemU) contact_rampD1(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.βrmp)
}

// TODO: improve these
func (o *ElemU) contact_f(x []float64) float64 {
	r := make([]float64, 3)
	o.Cell.Shp.InvMap(r, x, o.X)
	return o.Cell.Shp.CellBryDist(r)
}

func (o *ElemU) contact_g(x []float64) float64 {
	if false {
		return x[0] - 1.025
	}

	r := make([]float64, 3)
	Y := o.contact_get_Y()
	qua4 := shp.Get("qua4", 0) //o.Sim.GoroutineId)
	qua4.InvMap(r, x, Y)
	δ := o.Cell.Shp.CellBryDist(r)
	return δ
}

func (o *ElemU) contact_dgdx(dgdx, x []float64) {
	if false {
		dgdx[0], dgdx[1] = 1.0, 0.0
	}

	r := make([]float64, 3)
	Y := o.contact_get_Y()
	qua4 := shp.Get("qua4", 0) //o.Sim.GoroutineId)
	qua4.InvMap(r, x, Y)
	dfdR := make([]float64, 2)
	o.Cell.Shp.CellBryDistDeriv(dfdR, r)
	qua4.CalcAtR(Y, r, true)
	dgdx[0], dgdx[1] = 0.0, 0.0
	for i := 0; i < 2; i++ {
		for k := 0; k < 2; k++ {
			dgdx[i] += dfdR[k] * qua4.DRdx[k][i]
		}
	}
}

func (o *ElemU) contact_get_Y() (Y [][]float64) {
	test := 3
	switch test {
	case 1:
		m, l := 1.0, 1.0
		Y = [][]float64{
			{2.0 / m, 2.0/m + l, l, 0},
			{0, l / m, 2 + l/m, 2},
		}
	case 2:
		Y = [][]float64{
			{1.025, 2, 2, 1.025},
			{0, 0, 2, 2},
		}
	case 3:
		Y = [][]float64{
			{1.15, 2, 2, 1.025},
			{0, 0, 2, 2},
		}
	}
	return
}

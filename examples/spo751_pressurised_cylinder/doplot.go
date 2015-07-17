// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"math"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// filename
	filename, fnkey := io.Args0toFilename("spo751", ".sim", true)

	// constants
	nidx := 20 // selected node at outer surface
	didx := 0  // selected  dof index for plot
	nels := 4  // number of elements
	nips := 4  // number of ips

	// selected P values for stress plot
	Psel := []float64{100, 140, 180, 190}
	tolPsel := 2.0    // tolerance to compare P
	GPa2MPa := 1000.0 // conversion factor

	// input data
	Pcen := 200.0         // [Mpa]
	a, b := 100.0, 200.0  // [mm], [mm]
	E, ν := 210000.0, 0.3 // [MPa], [-]
	σy := 240.0           // [MPa]

	// analytical solution
	var sol ana.PressCylin
	sol.Init([]*fun.Prm{
		&fun.Prm{N: "a", V: a}, &fun.Prm{N: "b", V: b},
		&fun.Prm{N: "E", V: E}, &fun.Prm{N: "ν", V: ν},
		&fun.Prm{N: "σy", V: σy},
	})
	np := 41
	P_ana, Ub_ana := sol.CalcPressDisp(np)
	R_ana, Sr_ana, St_ana := sol.CalcStresses(Psel, np)

	// fem
	if !fem.Start(filename, false, false, false) {
		io.PfRed("Start failed\n")
		return
	}
	dom, sum, ok := fem.AllocSetAndInit(0, true, true)
	if !ok {
		io.PfRed("AllocSetAndInit failed\n")
		return
	}

	// gofem results
	nto := len(sum.OutTimes)
	P := make([]float64, nto)
	Ub := make([]float64, nto)
	R := utl.Deep3alloc(len(Psel), nels, nips)
	Sr := utl.Deep3alloc(len(Psel), nels, nips)
	St := utl.Deep3alloc(len(Psel), nels, nips)
	i := 0
	for tidx, t := range sum.OutTimes {

		// read results from file
		if !dom.In(sum, tidx, true) {
			io.PfRed("plot_spo751: cannot read solution\n")
			return
		}

		// collect results for load versus displacement plot
		nod := dom.Nodes[nidx]
		eq := nod.Dofs[didx].Eq
		P[tidx] = t * Pcen
		Ub[tidx] = dom.Sol.Y[eq]

		// stresses
		if isPsel(Psel, P[tidx], tolPsel) {
			for j, ele := range dom.ElemIntvars {
				e := ele.(*fem.ElemU)
				ipsdat := e.OutIpsData()
				for k, dat := range ipsdat {
					res := dat.Calc(dom.Sol)
					x, y := dat.X[0], dat.X[1]
					sx := res["sx"] * GPa2MPa
					sy := res["sy"] * GPa2MPa
					sxy := res["sxy"] * GPa2MPa / math.Sqrt2
					R[i][j][k], Sr[i][j][k], St[i][j][k], _ = ana.PolarStresses(x, y, sx, sy, sxy)
				}
			}
			i++
		}
	}

	// auxiliary data for plotting stresses
	colors := []string{"r", "m", "g", "k", "y", "c", "r", "m"}
	markers1 := []string{"o", "s", "x", ".", "^", "*", "o", "s"}
	markers2 := []string{"+", "+", "+", "+", "+", "+", "+", "+"}

	// plot load displacements
	plt.SetForPng(0.8, 400, 200)
	if true {
		//if false {
		plt.Plot(Ub_ana, P_ana, "'b-', ms=2, label='solution', clip_on=0")
		plt.Plot(Ub, P, "'r.--', label='fem: outer', clip_on=0")
		plt.Gll("$u_x\\;\\mathrm{[mm]}$", "$P\\;\\mathrm{[MPa]}$", "")
		plt.SaveD("/tmp", io.Sf("gofem_%s_disp.png", fnkey))
	}

	// plot radial stresses
	if true {
		//if false {
		plt.Reset()
		for i, Pval := range Psel {
			plt.Plot(R_ana, Sr_ana[i], "'b-'")
			for k := 0; k < nips; k++ {
				for j := 0; j < nels; j++ {
					args := io.Sf("'%s%s'", colors[i], markers1[i])
					if k > 1 {
						args = io.Sf("'k%s', ms=10", markers2[i])
					}
					if k == 0 && j == 0 {
						args += io.Sf(", label='P=%g'", Pval)
					}
					plt.PlotOne(R[i][j][k], Sr[i][j][k], args)
				}
			}
		}
		plt.Gll("$r\\;\\mathrm{[mm]}$", "$\\sigma_r\\;\\mathrm{[MPa]}$", "leg_loc='lower right'")
		plt.AxisXrange(a, b)
		plt.SaveD("/tmp", io.Sf("gofem_%s_sr.png", fnkey))
	}

	// plot tangential stresses
	if true {
		//if false {
		plt.Reset()
		for i, Pval := range Psel {
			plt.Plot(R_ana, St_ana[i], "'b-'")
			for k := 0; k < nips; k++ {
				for j := 0; j < nels; j++ {
					args := io.Sf("'%s%s'", colors[i], markers1[i])
					if k > 1 {
						args = io.Sf("'k%s', ms=10", markers2[i])
					}
					if k == 0 && j == 0 {
						args += io.Sf(", label='P=%g'", Pval)
					}
					plt.PlotOne(R[i][j][k], St[i][j][k], args)
				}
			}
		}
		plt.Gll("$r\\;\\mathrm{[mm]}$", "$\\sigma_t\\;\\mathrm{[MPa]}$", "leg_loc='upper left'")
		plt.SaveD("/tmp", io.Sf("gofem_%s_st.png", fnkey))
	}
}

func isPsel(Psel []float64, p, tol float64) bool {
	for _, pp := range Psel {
		if math.Abs(p-pp) < tol {
			return true
		}
	}
	return false
}

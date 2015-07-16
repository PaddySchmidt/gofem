// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ana"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_spo751a(tst *testing.T) {

	/* de Souza Neto, Perić and Owen, ex 7.5.1 p244
	 *
	 *                       22
	 *                        .
	 *                  19  ,' `.
	 *                    ,'     '.
	 *              17  ,'         \
	 *                .'            \
	 *           14 ,' `.            \ 21
	 *         12 ,'     \            '
	 *       9  .'        \            '
	 *     7  ,' `.        \ 16         '
	 *   4  .'     \        .           `
	 *  2  ' `.     \ 11     .          |
	 *    `.   \ 6   .       |          |
	 *     1.   .    |       |          |
	 *      |   |    |       |          |
	 *      -----------------------------
	 *      0 3 5 8 10  13  15    18   20
	 */

	//verbose()
	chk.PrintTitle("spo751a")

	// start simulation
	if !Start("data/spo751.sim", true, chk.Verbose, false) {
		tst.Errorf("test failed\n")
		return
	}

	// allocate domain and others
	if !Alloc(true) {
		tst.Errorf("Alloc failed\n")
		return
	}

	// set stage
	if !SetStage(0, true) {
		tst.Errorf("SetStage failed\n")
		return
	}

	// domain
	dom := Global.Domains[0]

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 23)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 5, 7, 2, 3, 6, 4, 1, 10, 12, 8, 11, 9, 15, 17, 13, 16, 14, 20, 22, 18, 21, 19})
	chk.Ints(tst, "eqs", eqs, utl.IntRange(23*2))

	// check solution arrays
	ny := 23 * 2
	nλ := 9 + 9
	nyb := ny + nλ
	chk.IntAssert(len(dom.Sol.Y), ny)
	chk.IntAssert(len(dom.Sol.Dydt), 0)
	chk.IntAssert(len(dom.Sol.D2ydt2), 0)
	chk.IntAssert(len(dom.Sol.Psi), 0)
	chk.IntAssert(len(dom.Sol.Zet), 0)
	chk.IntAssert(len(dom.Sol.Chi), 0)
	chk.IntAssert(len(dom.Sol.L), nλ)
	chk.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	chk.IntAssert(len(dom.Fb), nyb)
	chk.IntAssert(len(dom.Wb), nyb)

	// check umap
	umaps := [][]int{
		{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
		{2, 3, 16, 17, 18, 19, 4, 5, 20, 21, 22, 23, 24, 25, 10, 11},
		{16, 17, 26, 27, 28, 29, 18, 19, 30, 31, 32, 33, 34, 35, 22, 23},
		{26, 27, 36, 37, 38, 39, 28, 29, 40, 41, 42, 43, 44, 45, 32, 33},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemU)
		io.Pforan("e%d.umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "umap", e.Umap, umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_uy_eqs []int // constrained uy equations [sorted]
	var ct_incsup_xeqs []int
	var ct_incsup_yeqs []int
	αrad := 120.0 * math.Pi / 180.0
	cα, sα := math.Cos(αrad), math.Sin(αrad)
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssertLessThanOrEqualTo(1, len(c.Eqs)) // 1 ≤ neqs
		io.Pforan("c.Key=%s c.Eqs=%v\n", c.Key, c.Eqs)
		if len(c.Eqs) == 1 {
			if c.Key == "uy" {
				ct_uy_eqs = append(ct_uy_eqs, c.Eqs[0])
				continue
			}
		} else {
			if c.Key == "incsup" {
				ct_incsup_xeqs = append(ct_incsup_xeqs, c.Eqs[0])
				ct_incsup_yeqs = append(ct_incsup_yeqs, c.Eqs[1])
				chk.AnaNum(tst, "cos(α)", 1e-15, c.ValsA[0], cα, false)
				chk.AnaNum(tst, "sin(α)", 1e-15, c.ValsA[1], sα, false)
				continue
			}
		}
		tst.Errorf("key %s is incorrect", c.Key)
	}
	sort.Ints(ct_uy_eqs)
	sort.Ints(ct_incsup_xeqs)
	sort.Ints(ct_incsup_yeqs)
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 3, 9, 17, 21, 27, 31, 37, 41})
	chk.Ints(tst, "incsup x equations", ct_incsup_xeqs, []int{4, 6, 12, 18, 24, 28, 34, 38, 44})
	chk.Ints(tst, "incsup y equations", ct_incsup_yeqs, []int{5, 7, 13, 19, 25, 29, 35, 39, 45})
}

func Test_spo751b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("spo751b")

	// run simulation
	if !Start("data/spo751.sim", true, chk.Verbose, false) {
		tst.Errorf("test failed\n")
		return
	}

	// for debugging Kb
	//if true {
	if false {
		defer u_DebugKb(&testKb{
			tst: tst, eid: 3, tol: 1e-5, verb: chk.Verbose,
			ni: 1, nj: 1, itmin: 1, itmax: -1, tmin: 0.89, tmax: 0.96,
		})()
	}

	// run simulation
	if !RunAll() {
		tst.Errorf("test failed\n")
		return
	}

	// check
	if true {
		verb := false
		skipK := true
		tolK := 1e-17
		tolu := 1e-12
		tols := 1e-14
		TestingCompareResultsU(tst, "data/spo751.sim", "cmp/spo751.cmp", tolK, tolu, tols, skipK, verb)
	}

	// plot
	//if true {
	if false {
		plot_spo751("spo751")
	}
}

func Test_spo751re(tst *testing.T) {

	verbose()
	chk.PrintTitle("spo751re. Richardson extrapolation")

	// run simulation
	if !Start("data/spo751re.sim", true, chk.Verbose, false) {
		io.Pfred("start failed\n")
		return
	}

	// run simulation
	if !RunAll() {
		io.Pfred("run failed\n")
	}

	// plot
	//if true {
	if false {
		plot_spo751("spo751re")
	}
}

func plot_spo751(fnkey string) {

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

	// read summary
	sum := ReadSum(Global.Dirout, Global.Fnkey)

	// allocate domain
	distr := false
	d := NewDomain(Global.Sim.Regions[0], distr)
	if !d.SetStage(0, Global.Sim.Stages[0], distr) {
		io.PfRed("plot_spo751: SetStage failed\n")
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
		if !d.In(sum, tidx, true) {
			io.PfRed("plot_spo751: cannot read solution\n")
			return
		}

		// collect results for load versus displacement plot
		nod := d.Nodes[nidx]
		eq := nod.Dofs[didx].Eq
		P[tidx] = t * Pcen
		Ub[tidx] = d.Sol.Y[eq]

		// stresses
		if isPsel(Psel, P[tidx], tolPsel) {
			for j, ele := range d.ElemIntvars {
				e := ele.(*ElemU)
				ipsdat := e.OutIpsData()
				for k, dat := range ipsdat {
					res := dat.Calc(d.Sol)
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
	plt.SetForEps(0.8, 300)
	if true {
		//if false {
		plt.Plot(Ub_ana, P_ana, "'b-', ms=2, label='solution', clip_on=0")
		plt.Plot(Ub, P, "'r.--', label='fem: outer', clip_on=0")
		plt.Gll("$u_x\\;\\mathrm{[mm]}$", "$P\\;\\mathrm{[MPa]}$", "")
		plt.SaveD("/tmp", io.Sf("gofem_%s_disp.eps", fnkey))
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
		plt.SaveD("/tmp", io.Sf("gofem_%s_sr.eps", fnkey))
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
		plt.SaveD("/tmp", io.Sf("gofem_%s_st.eps", fnkey))
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

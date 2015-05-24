// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"path"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// input data
	simfn := "elast.sim"
	matname := "lrm1"
	pcmax := 30.0
	npts := 101

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		matname = flag.Arg(1)
	}
	if len(flag.Args()) > 2 {
		pcmax = io.Atof(flag.Arg(2))
	}
	if len(flag.Args()) > 3 {
		npts = io.Atoi(flag.Arg(3))
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfn   = %30s // simulation filename\n", simfn)
	io.Pf("  matname = %30s // material name\n", matname)
	io.Pf("  pcmax   = %30v // max pc\n", pcmax)
	io.Pf("  npts    = %30v // number of points\n", npts)
	io.Pf("\n")

	// load simulation
	sim := inp.ReadSim("", simfn, "lrm_", false)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.Mdb.Get(matname)
	if mat == nil {
		io.PfRed("cannot get material\n")
		return
	}
	io.Pforan("mat = %v\n", mat)

	// get and initialise model
	mdl := mreten.GetModel(simfn, matname, mat.Model, false)
	if mdl == nil {
		io.PfRed("cannot allocate model\n")
		return
	}
	mdl.Init(mat.Prms)

	// plot drying path
	d_Pc := utl.LinSpace(0, pcmax, npts)
	d_Sl := make([]float64, npts)
	d_Sl[0] = 1
	var err error
	for i := 1; i < npts; i++ {
		d_Sl[i], err = mreten.Update(mdl, d_Pc[i-1], d_Sl[i-1], d_Pc[i]-d_Pc[i-1])
		if err != nil {
			io.PfRed("drying: cannot updated model\n%v\n", err)
			return
		}
	}
	plt.Plot(d_Pc, d_Sl, io.Sf("'b-', label='%s (dry)', clip_on=0", matname))

	// plot wetting path
	w_Pc := utl.LinSpace(pcmax, 0, npts)
	w_Sl := make([]float64, npts)
	w_Sl[0] = d_Sl[npts-1]
	for i := 1; i < npts; i++ {
		w_Sl[i], err = mreten.Update(mdl, w_Pc[i-1], w_Sl[i-1], w_Pc[i]-w_Pc[i-1])
		if err != nil {
			io.PfRed("wetting: cannot updated model\n%v\n", err)
			return
		}
	}
	plt.Plot(w_Pc, w_Sl, io.Sf("'c-', label='%s (wet)', clip_on=0", matname))

	// save results
	type Results struct{ Pc, Sl []float64 }
	res := Results{append(d_Pc, w_Pc...), append(d_Sl, w_Sl...)}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	err = enc.Encode(&res)
	if err != nil {
		io.PfRed("cannot encode results\n")
		return
	}
	fn := path.Join(sim.Data.DirOut, matname+".dat")
	io.WriteFile(fn, &buf)
	io.Pf("file <[1;34m%s[0m> written\n", fn)

	// show figure
	plt.AxisYrange(0, 1)
	plt.Cross()
	plt.Gll("$p_c$", "$s_{\\ell}$", "")
	plt.Show()
}

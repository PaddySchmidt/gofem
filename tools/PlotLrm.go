// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"
	"encoding/json"
	"path"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data
	simfn, _ := io.ArgToFilename(0, "elast", ".sim", true)
	matname := io.ArgToString(1, "lrm1")
	pcmax := io.ArgToFloat(2, 30.0)
	npts := io.ArgToInt(3, 101)

	// print input table
	io.Pf("\n%s\n", io.ArgsTable(
		"simulation filename", "simfn", simfn,
		"material name", "matname", matname,
		"max pc", "pcmax", pcmax,
		"number of points", "npts", npts,
	))

	// load simulation
	sim := inp.ReadSim(simfn, "lrm", false, 0)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.MatParams.Get(matname)
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
	plt.Cross("")
	plt.Gll("$p_c$", "$s_{\\ell}$", "")
	plt.Show()
}

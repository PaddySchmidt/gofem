// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// Plot plots retention model
//  args1 -- arguments for model computed by solving differential equation; e.g. "'b*-'"
//           if args1 == "", plot is skiped
//  args2 -- arguments for model computed by directly calling sl(pc), if available
//           if args2 == "", plot is skiped
func Plot(mdl Model, pc0, sl0, pcf float64, npts int, args1, args2, label string) (err error) {

	// plot using Update
	Pc := utl.LinSpace(pc0, pcf, npts)
	Sl := make([]float64, npts)
	if args1 != "" {
		Sl[0] = sl0
		for i := 1; i < npts; i++ {
			Sl[i], err = Update(mdl, Pc[i-1], Sl[i-1], Pc[i]-Pc[i-1])
			if err != nil {
				return
			}
		}
		plt.Plot(Pc, Sl, io.Sf("%s, label='%s', clip_on=0", args1, label))
	}

	// plot using Sl function
	if args2 != "" {
		if m, ok := mdl.(Nonrate); ok {
			Pc = utl.LinSpace(pc0, pcf, 101)
			Sl = make([]float64, 101)
			for i, pc := range Pc {
				Sl[i] = m.Sl(pc)
			}
			plt.Plot(Pc, Sl, io.Sf("%s, label='%s_direct', clip_on=0", args2, label))
		}
	}
	return
}

// PlotEnd ends plot and show figure, if show==true
func PlotEnd(show bool) {
	plt.AxisYrange(0, 1)
	plt.Cross()
	plt.Gll("$p_c$", "$s_{\\ell}$", "")
	if show {
		plt.Show()
	}
}

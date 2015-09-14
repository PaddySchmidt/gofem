// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ana implements analytical solutions
package ana

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// PlateHole implements Kirsch's solution to 2D Plate with hole
//
//      y ^
//        |    qnV
//        ↓↓↓↓↓↓↓↓↓↓↓↓↓
//        ------------- ←
//        |           | ←
//       ▷|      E    | ←
//        |      ν    | ← qnH
//        `'-.        | ←
//            \       | ←
//         r   -------- ← ---> x
//                △
type PlateHole struct {

	// input
	r       float64 // radius
	E       float64 // Young's modulus
	ν       float64 // Poisson's coefficient
	qnV     float64 // vertical distributed load
	qnH     float64 // horizontal distributed load
	pstress bool    // plane stress
}

// Init initialises this structure
func (o *PlateHole) Init(prms fun.Prms) {

	// default values
	o.r = 1.0
	o.E = 1e5
	o.ν = 0.3
	o.qnV = 0.0
	o.qnH = 10.0

	// parameters
	for _, p := range prms {
		switch p.N {
		case "r":
			o.r = p.V
		case "E":
			o.E = p.V
		case "nu":
			o.ν = p.V
		case "qnV":
			o.qnV = p.V
		case "qnH":
			o.qnH = p.V
		case "pstress":
			o.pstress = p.V > 0
		}
	}
}

// Stress computes stresses @ (x,y) and time t
func (o *PlateHole) Stress(t float64, x []float64) (sx, sy, sz, sxy float64) {

	// auxiliary
	ph := o.qnH * t
	pv := o.qnV * t

	// polar coordinates
	d := math.Sqrt(x[0]*x[0] + x[1]*x[1])
	c, s := x[0]/d, x[1]/d
	cc, ss := c*c, s*s
	cs := c * s
	c2t := cc - ss
	s2t := 2.0 * c * s

	// solution in polar coordinates
	pm := (ph + pv) / 2.0
	pd := (ph - pv) / 2.0
	b := o.r * o.r / (d * d)
	sr := pm*(1.0-b) + pd*(1.0-4.0*b+3.0*b*b)*c2t
	st := pm*(1.0+b) - pd*(1.0+3.0*b*b)*c2t
	srt := -pd * (1.0 + 2.0*b - 3.0*b*b) * s2t

	// rotation to x-y coordinates
	sx = cc*sr + ss*st - 2.0*cs*srt
	sy = ss*sr + cc*st + 2.0*cs*srt
	sxy = cs*sr - cs*st + (cc-ss)*srt

	// out-of-plane stress
	if o.pstress {
		sz = o.ν * (sx + sy)
	}
	return
}

// CompareStress compares stresses
//  Output:
//   e -- L² error for each component
func (o *PlateHole) CompareStress(t float64, x, σ []float64, tol float64, verbose bool) (e []float64) {

	// analytical solution
	sx, sy, sz, sxy := o.Stress(t, x)

	// message
	if verbose {
		chk.PrintAnaNum("σx ", tol, sx, σ[0], verbose)
		chk.PrintAnaNum("σy ", tol, sy, σ[1], verbose)
		chk.PrintAnaNum("σz ", tol, sz, σ[2], verbose)
		chk.PrintAnaNum("σxy", tol, sxy, σ[3], verbose)
	}

	// check stresses
	e = []float64{
		math.Abs(sx - σ[0]),
		math.Abs(sy - σ[1]),
		math.Abs(sz - σ[2]),
		math.Abs(sxy - σ[3]),
	}
	return
}

// PlotStress plots stresses along y=0 (horizontal line) and x=0 (vertical line)
func (o *PlateHole) PlotStress(t, L float64, npts int) {

	d := utl.LinSpace(o.r, L, npts)
	Sx := make([]float64, npts)
	Sy := make([]float64, npts)
	Sxy := make([]float64, npts)

	plt.Subplot(2, 1, 1)
	for i := 0; i < npts; i++ {
		Sx[i], Sy[i], _, Sxy[i] = o.Stress(t, []float64{d[i], 0}) // y=0
	}
	plt.Plot(d, Sx, "color='r',label='$\\sigma_x$ @ $y=0$'")
	plt.Plot(d, Sy, "color='g',label='$\\sigma_y$ @ $y=0$'")
	plt.Plot(d, Sxy, "color='b',label='$\\sigma_{xy}$ @ $y=0$'")
	plt.Gll("$x$", "stresses", "")

	plt.Subplot(2, 1, 2)
	for i := 0; i < npts; i++ {
		Sx[i], Sy[i], _, Sxy[i] = o.Stress(t, []float64{0, d[i]}) // x=0
	}
	plt.Plot(Sx, d, "color='r',label='$\\sigma_x$ @ $x=0$'")
	plt.Plot(Sy, d, "color='g',label='$\\sigma_y$ @ $x=0$'")
	plt.Plot(Sxy, d, "color='b',label='$\\sigma_{xy}$ @ $x=0$'")
	plt.Gll("stresses", "$y$", "")
}

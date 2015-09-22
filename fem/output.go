// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// PlotAllBendingMoments plots all bending moments
//  Input:
//   dom       -- Domain
//   nstations -- number of stations
//   withtext  -- show bending moment values
//   numfmt    -- number format for values
//   tolM      -- tolerance to clip absolute values of M
//   coef      -- coefficient to scale max(dimension) divided by max(Y); e.g. 0.1
//  Output:
//   beams -- all beam elements
//   allM  -- bending moments corresponding to all beams
func PlotAllBendingMoments(dom *Domain, nstations int, withtext bool, numfmt string, tolM, coef float64) (beams []*Beam, allM [][]float64) {

	// collect beams
	for _, elem := range dom.Elems {
		if beam, ok := elem.(*Beam); ok {
			beams = append(beams, beam)
		}
	}

	// compute bending moments
	allM = make([][]float64, len(beams))
	for i, beam := range beams {
		_, allM[i] = beam.CalcVandM(dom.Sol, 0, nstations)
	}

	// scaling factor
	maxAbsM := la.MatLargest(allM, 1)
	dist := utl.Max(dom.Msh.Xmax-dom.Msh.Xmin, dom.Msh.Ymax-dom.Msh.Ymin)
	sf := 1.0
	if maxAbsM > 1e-7 {
		sf = coef * dist / maxAbsM
	}

	// draw
	dom.Msh.Draw2d()
	for i, beam := range beams {
		beam.PlotDiagMoment(allM[i], withtext, numfmt, tolM, sf)
	}
	return
}

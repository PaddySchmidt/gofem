// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ana implements analytical solutions
package ana

import "math"

// PolarStresses computes stress components w.r.t polar system from given
// Cartesian components
func PolarStresses(x, y, sx, sy, sxy float64) (r, sr, st, srt float64) {
	r = math.Sqrt(x*x + y*y)
	β := math.Atan2(y, x)
	si, co := math.Sin(β), math.Cos(β)
	ss, cc, cs := si*si, co*co, co*si
	sr = cc*sx + ss*sy + 2.0*cs*sxy
	st = ss*sx + cc*sy - 2.0*cs*sxy
	srt = -cs*sx + cs*sy + (cc-ss)*sxy
	return
}

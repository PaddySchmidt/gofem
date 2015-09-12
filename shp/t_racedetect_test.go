// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_race01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("Test race01")

	nchan := 2
	done := make(chan int, nchan)

	shapes := make([]*Shape, nchan)
	for i := 0; i < nchan; i++ {
		shapes[i] = Get("tri3", i)
	}
	io.Pforan("shapes = %v\n", shapes)

	for i := 0; i < nchan; i++ {
		go func(shape *Shape) {
			shape.CalcAtR([][]float64{
				{0, 1, 0},
				{0, 0, 1},
			}, []float64{0.5, 0.5, 0}, true)
			done <- 1
		}(shapes[i])
	}

	for i := 0; i < nchan; i++ {
		<-done
	}
}

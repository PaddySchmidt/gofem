// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// PltEntity stores all data for a plot entity (X vs Y)
type PltEntity struct {
	Alias string    // alias
	X     []float64 // x-values
	Y     []float64 // y-values
	Xlbl  string    // horizontal axis label (raw; e.g. "t")
	Ylbl  string    // vertical axis label (raw; e.g. "pl")
	Style plt.Fmt   // style
}

// SplotDat stores all data for one subplot
type SplotDat struct {
	Title  string       // title of subplot
	Topts  string       // title options
	Xscale float64      // x-axis scale
	Yscale float64      // y-axis scale
	Xlbl   string       // x-axis label (formatted; e.g. "$t$")
	Ylbl   string       // y-axis label (formatted; e.g. "$p_{\ell}$")
	Data   []*PltEntity // data and styles to be plotted
}

// Splot activates a new subplot window
func Splot(splotTitle string) {
	s := &SplotDat{Title: splotTitle}
	Splots = append(Splots, s)
	Csplot = s
}

// SplotConfig configures units and scales of axes
func SplotConfig(xunit, yunit string, xscale, yscale float64) {
	if Csplot != nil {
		var xlabel, ylabel string
		if len(Csplot.Data) > 0 {
			xlabel = Csplot.Data[0].Xlbl
			ylabel = Csplot.Data[0].Ylbl
		}
		Csplot.Xlbl = GetTexLabel(xlabel, xunit)
		Csplot.Ylbl = GetTexLabel(ylabel, yunit)
		Csplot.Xscale = xscale
		Csplot.Yscale = yscale
	}
}

// Plot plots data
//  xHandle -- can be a string, e.g. "t" or a slice, e.g. pc = []float64{0, 1, 2}
//  yHandle -- can be a string, e.g. "pl" or a slice, e.g. sl = []float64{0, 1, 2}
//  alias   -- alias such as "centre"
//  fm      -- formatting codes; e.g. plt.Fmt{C:"blue", L:"label"}
//  idxI    -- index of time; use -1 for all times
func Plot(xHandle, yHandle interface{}, alias string, fm plt.Fmt, idxI int) {
	var e PltEntity
	e.Alias = alias
	e.Style = fm
	e.X, e.Xlbl = get_vals_and_labels(xHandle, yHandle, alias, idxI)
	e.Y, e.Ylbl = get_vals_and_labels(yHandle, xHandle, alias, idxI)
	if len(e.X) != len(e.Y) {
		chk.Panic("lengths of x- and y-series are different. len(x)=%d, len(y)=%d, x=%v, y=%v", len(e.X), len(e.Y), xHandle, yHandle)
	}
	if Csplot == nil {
		Splot("")
	}
	Csplot.Data = append(Csplot.Data, &e)
	SplotConfig("", "", 1, 1)
}

// ExtraPlt defines a callback function for extra plt commands
//  Note: i and j are indices as in Subplot
type ExtraPlt func(i, j, nplots int)

// Draw draws or save figure with plot
//  dirout -- directory to save figure
//  fname  -- file name; e.g. myplot.eps or myplot.png. Use "" to skip saving
//  show   -- shows figure
//  extra  -- is called just after Subplot command and before any plotting
func Draw(dirout, fname string, show bool, extra ExtraPlt) {
	nplots := len(Splots)
	nr, nc := utl.BestSquare(nplots)
	var k int
	for i := 0; i < nr; i++ {
		for j := 0; j < nc; j++ {
			plt.Subplot(nr, nc, k+1)
			if extra != nil {
				extra(i+1, j+1, nplots)
			}
			if Splots[k].Title != "" {
				plt.Title(Splots[k].Title, Splots[k].Topts)
			}
			data := Splots[k].Data
			for _, d := range data {
				if d.Style.L == "" {
					d.Style.L = d.Alias
				}
				plt.Plot(d.X, d.Y, d.Style.GetArgs("clip_on=0"))
			}
			plt.Gll(Splots[k].Xlbl, Splots[k].Ylbl, "")
			k += 1
		}
	}
	if fname != "" {
		plt.SaveD(dirout, fname)
	}
	if show {
		plt.Show()
	}
}

// auxiliary /////////////////////////////////////////////////////////////////////////////////////////

func get_vals_and_labels(handle, otherHandle interface{}, alias string, idxI int) ([]float64, string) {
	otherKey := "any"
	if key, ok := otherHandle.(string); ok {
		otherKey = key
	}
	switch hnd := handle.(type) {
	case []float64:
		return hnd, io.Sf("%s-type", alias)
	case string:
		switch hnd {
		case "t":
			return T, "t"
		case "x":
			xcoords, _, _ := GetXYZ(otherKey, alias)
			return xcoords, "x"
		case "y":
			_, ycoords, _ := GetXYZ(otherKey, alias)
			return ycoords, "y"
		case "z":
			_, _, zcoords := GetXYZ(otherKey, alias)
			return zcoords, "z"
		case "dist":
			return GetDist(otherKey, alias), "dist"
		}
		return GetRes(hnd, alias, idxI), hnd
	}
	chk.Panic("cannot get values slice with handle = %v", handle)
	return nil, ""
}

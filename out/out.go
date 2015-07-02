// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package out implements FE simulation output handling for analyses and plotting
package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// constants
var (
	TolC = 1e-8 // tolerance to compare x-y-z coordinates
	TolT = 1e-3 // tolerance to compare times
	Ndiv = 20   // bins n-division
)

// ResultsMap maps aliases to points
type ResultsMap map[string]Points

// Global variables
var (

	// data set by Start
	Sum       *fem.Summary     // summary
	Dom       *fem.Domain      // FE domain
	Ipoints   []*fem.OutIpData // all integration points. ipid == index in Ipoints
	Cid2ips   [][]int          // [ncells][nip] maps cell id to index in Ipoints
	Ipkey2ips map[string][]int // maps ip keys to indices in Ipoints
	Ipkeys    map[string]bool  // all ip keys
	NodBins   gm.Bins          // bins for nodes
	IpsBins   gm.Bins          // bins for integration points
	IpsMin    []float64        // [ndim] {x,y,z}_min among all ips
	IpsMax    []float64        // [ndim] {x,y,z}_max among all ips

	// defined entities and results loaded by LoadResults
	Planes   map[string]*PlaneData // for points defined on planes. maps aliases to data
	Results  ResultsMap            // maps labels => points
	TimeInds []int                 // selected output indices
	Times    []float64             // selected output times

	// extrapolated values
	Extrap []string             // keys to be extrapolated; e.g. []string{"nwlx", "nwly"}
	ExVals []map[string]float64 // [nverts][nkeys] extrapolated values

	// subplots
	Splots []*SplotDat // all subplots
	Csplot *SplotDat   // current subplot
)

// End must be called at the end to flush log file
func End() {
	if err := recover(); err != nil {
		io.PfRed("ERROR: %v\n", err)
	} else {
		fem.End()
	}
}

// Start starts handling of results given a simulation input file
func Start(simfnpath string, stageIdx, regionIdx int) {

	// start FE global structure
	fem.Global.LogPrefix = "out_"
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		chk.Panic("cannot start analysis process with simfnpath=%q\n", simfnpath)
	}

	// read summary
	Sum = fem.ReadSum(fem.Global.Dirout, fem.Global.Fnkey)
	if Sum == nil {
		chk.Panic("cannot read summary file for simulation=%q\n", simfnpath)
	}

	// allocate domain
	distr := false
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx], distr)
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx], distr) {
		chk.Panic("cannot allocate domain\n")
	}

	// clear previous data
	Ipoints = make([]*fem.OutIpData, 0)
	Cid2ips = make([][]int, len(Dom.Msh.Cells))
	Ipkey2ips = make(map[string][]int)
	Ipkeys = make(map[string]bool)
	Planes = make(map[string]*PlaneData)
	Results = make(map[string]Points)
	TimeInds = make([]int, 0)
	Times = make([]float64, 0)
	Splots = make([]*SplotDat, 0)

	// bins
	m := Dom.Msh
	δ := TolC * 2
	xi := []float64{m.Xmin - δ, m.Ymin - δ}
	xf := []float64{m.Xmax + δ, m.Ymax + δ}
	if m.Ndim == 3 {
		xi = append(xi, m.Zmin-δ)
		xf = append(xf, m.Zmax+δ)
	}
	err := NodBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for nodes: %v", err)
	}
	err = IpsBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for integration points: %v", err)
	}

	// add nodes to bins
	for _, nod := range Dom.Nodes {
		err := NodBins.Append(nod.Vert.C, nod.Vert.Id)
		if err != nil {
			return
		}
	}

	// to find limits
	ndim := Dom.Msh.Ndim
	IpsMin = make([]float64, ndim)
	IpsMax = make([]float64, ndim)
	first := true

	// add integration points to slice of ips and to bins
	for cid, ele := range Dom.Cid2elem {
		if ele == nil {
			continue
		}
		dat := ele.OutIpsData()
		nip := len(dat)
		ids := make([]int, nip)
		for i, d := range dat {

			// add to bins
			ipid := len(Ipoints)
			ids[i] = ipid
			Ipoints = append(Ipoints, d)
			err = IpsBins.Append(d.X, ipid)
			if err != nil {
				chk.Panic("cannot append to bins of integration points: %v", err)
			}

			// set auxiliary map
			vals := d.Calc(Dom.Sol)
			for key, _ := range vals {
				utl.StrIntsMapAppend(&Ipkey2ips, key, ipid)
				Ipkeys[key] = true
			}

			// limits
			if first {
				for j := 0; j < ndim; j++ {
					IpsMin[j] = d.X[j]
					IpsMax[j] = d.X[j]
				}
				first = false
			} else {
				for j := 0; j < ndim; j++ {
					IpsMin[j] = min(IpsMin[j], d.X[j])
					IpsMax[j] = max(IpsMax[j], d.X[j])
				}
			}
		}
		Cid2ips[cid] = ids
	}
}

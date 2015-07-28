// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			if mpi.Rank() == 0 {
				chk.Verbose = true
				for i := 8; i > 3; i-- {
					chk.CallerInfo(i)
				}
				io.PfRed("ERROR: %v\n", err)
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// default input parameters

	// read input parameters
	fnamepath, _ := io.ArgToFilename(0, "", ".sim", true)
	verbose := io.ArgToBool(1, true)
	erasePrev := io.ArgToBool(2, true)
	saveSummary := io.ArgToBool(3, true)
	allowParallel := io.ArgToBool(4, true)
	alias := io.ArgToString(5, "")

	// message
	if mpi.Rank() == 0 && verbose {
		io.PfWhite("\nGofem v3 -- Go Finite Element Method\n\n")
		io.Pf("Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.\n")
		io.Pf("Use of this source code is governed by a BSD-style\n")
		io.Pf("license that can be found in the LICENSE file.\n\n")

		io.Pf("\n%v\n", io.ArgsTable(
			"filename path", "fnamepath", fnamepath,
			"show messages", "verbose", verbose,
			"erase previous results", "erasePrev", erasePrev,
			"save summary", "saveSummary", saveSummary,
			"allow parallel run", "allowParallel", allowParallel,
			"word to add to results", "alias", alias,
		))
	}

	// profiling?
	defer utl.DoProf(false)()

	// analysis data
	readSummary := false
	analysis := fem.NewFEM(fnamepath, alias, erasePrev, saveSummary, readSummary, allowParallel, verbose)

	// run simulation
	err := analysis.Run()
	if err != nil {
		chk.Panic("Run failed:\n%v", err)
	}
}

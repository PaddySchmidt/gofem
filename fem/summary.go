// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"os"
	"path"

	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// Summary records summary of outputs
type Summary struct {

	// main data
	Nproc    int          // number of processors used in last last run; equal to 1 if not distributed
	OutTimes []float64    // [nOutTimes] output times
	Resids   utl.DblSlist // residuals (if Stat is on; includes all stages)
	Dirout   string       // directory where results are stored
	Fnkey    string       // filename key of simulation

	// auxiliary
	tidx int // time output index
}

// SaveResults save the results from all domains (nodes and elements)
func (o *Summary) SaveResults() (ok bool) {

	// output results from all domains
	for _, d := range Global.Domains {
		if LogErrCond(!d.Out(o.tidx), "SaveResults failed") {
			break
		}
	}
	if Stop() {
		return
	}

	// update internal structures
	o.OutTimes = append(o.OutTimes, Global.Time)
	o.tidx += 1

	// success
	return true
}

// SaveSums saves summary to disc
func (o Summary) Save() (ok bool) {

	// set flags before saving
	o.Nproc = Global.Nproc
	o.Dirout = Global.Dirout
	o.Fnkey = Global.Fnkey

	// skip if not root
	if !Global.Root {
		return true
	}

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode summary
	if LogErr(enc.Encode(o), "SaveSum") {
		return
	}

	// save file
	fn := out_sum_path(Global.Dirout, Global.Fnkey, Global.Rank)
	return save_file("SaveSum", "summary", fn, &buf)
}

// ReadSum reads summary back
//  Note: returns nil on errors
func ReadSum(dir, fnkey string) (o *Summary) {

	// open file
	fn := out_sum_path(dir, fnkey, 0) // reading always from proc # 0
	fil, err := os.Open(fn)
	if LogErr(err, "ReadSum") {
		return nil
	}
	defer func() {
		LogErr(fil.Close(), "ReadSum: cannot close file")
	}()

	// decode summary
	var sum Summary
	dec := GetDecoder(fil)
	err = dec.Decode(&sum)
	if LogErr(err, "ReadSum") {
		return nil
	}
	return &sum
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_sum_path(dir, fnkey string, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_sum.%s", fnkey, proc, Global.Enc))
}

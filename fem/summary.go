// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"os"
	"path"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// Summary records summary of outputs
type Summary struct {

	// main data
	Dirout   string       // directory where results are stored
	Fnkey    string       // filename key of simulation
	Nproc    int          // number of processors
	OutTimes []float64    // [nOutTimes] output times
	Resids   utl.DblSlist // residuals (if Stat is on; includes all stages)

	// auxiliary
	tidx int // time output index
}

// SaveDomains save the results from all domains (nodes and elements)
func (o *Summary) SaveDomains(time float64, doms []*Domain, verbose bool) (err error) {

	// output results from all domains
	for _, d := range doms {
		err = d.Save(o.tidx, verbose)
		if err != nil {
			return chk.Err("SaveResults failed:\n%v", err)
		}
	}

	// update internal structures
	o.OutTimes = append(o.OutTimes, time)
	o.tidx += 1
	return
}

// SaveSums saves summary to disc
func (o Summary) Save(dirout, fnkey, enctype string, nproc, proc int, verbose bool) (err error) {

	// skip if not root
	if proc != 0 {
		return
	}

	// set data
	o.Dirout = dirout
	o.Fnkey = fnkey
	o.Nproc = nproc

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf, enctype)

	// encode summary
	err = enc.Encode(o)
	if err != nil {
		return chk.Err("encoding of summary failed:\n%v", err)
	}

	// save file
	fn := out_sum_path(dirout, fnkey, enctype, proc)
	return save_file(fn, &buf, verbose)
}

// Read reads summary back
func (o *Summary) Read(dir, fnkey, enctype string) (err error) {

	// open file
	fn := out_sum_path(dir, fnkey, enctype, 0) // reading always from proc # 0
	fil, err := os.Open(fn)
	if err != nil {
		return
	}
	defer func() { err = fil.Close() }()

	// decode summary
	dec := GetDecoder(fil, enctype)
	err = dec.Decode(o)
	if err != nil {
		return chk.Err("cannot decode summary:\n%v", err)
	}
	return
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_sum_path(dir, fnkey, enctype string, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_sum.%s", fnkey, proc, enctype))
}

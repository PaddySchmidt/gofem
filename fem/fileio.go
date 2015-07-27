// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"encoding/gob"
	"encoding/json"
	goio "io"
	"os"
	"path"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

// Encoder defines encoders; e.g. gob or json
type Encoder interface {
	Encode(e interface{}) error
}

// Decoder defines decoders; e.g. gob or json
type Decoder interface {
	Decode(e interface{}) error
}

// GetEncoder returns a new encoder
func GetEncoder(w goio.Writer, enctype string) Encoder {
	if enctype == "json" {
		return json.NewEncoder(w)
	}
	return gob.NewEncoder(w)
}

// GetDecoder returns a new decoder
func GetDecoder(r goio.Reader, enctype string) Decoder {
	if enctype == "json" {
		return json.NewDecoder(r)
	}
	return gob.NewDecoder(r)
}

// SaveSol saves solution (o.Sol) to a file which name is set with tidx (time output index)
func (o Domain) SaveSol(tidx int, verbose bool) (err error) {

	// skip if root
	if o.Proc != 0 {
		return
	}

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf, o.Sim.EncType)

	// encode Sol
	err = enc.Encode(o.Sol.T)
	if err != nil {
		return chk.Err("cannot encode Domain.Sol.T\n%v", err)
	}
	err = enc.Encode(o.Sol.Y)
	if err != nil {
		return chk.Err("cannot encode Domain.Sol.Y\n%v", err)
	}
	err = enc.Encode(o.Sol.Dydt)
	if err != nil {
		return chk.Err("cannot encode Domain.Sol.Dydt\n%v", err)
	}
	err = enc.Encode(o.Sol.D2ydt2)
	if err != nil {
		return chk.Err("cannot encode Domain.Sol.D2ydt2\n%v", err)
	}

	// save file
	fn := out_nod_path(o.Sim.DirOut, o.Sim.Key, o.Sim.EncType, tidx, o.Proc)
	return save_file(fn, &buf, verbose)
}

// ReadSol reads Solution from a file which name is set with tidx (time output index)
func (o *Domain) ReadSol(dir, fnkey, enctype string, tidx int) (err error) {

	// open file
	fn := out_nod_path(dir, fnkey, enctype, tidx, 0) // 0 => reading always from proc # 0
	fil, err := os.Open(fn)
	if err != nil {
		return
	}
	defer func() { err = fil.Close() }()

	// get decoder
	dec := GetDecoder(fil, enctype)

	// decode Sol
	err = dec.Decode(&o.Sol.T)
	if err != nil {
		return chk.Err("cannot decode Domain.Sol.T\n%v", err)
	}
	err = dec.Decode(&o.Sol.Y)
	if err != nil {
		return chk.Err("cannot decode Domain.Sol.Y\n%v", err)
	}
	err = dec.Decode(&o.Sol.Dydt)
	if err != nil {
		return chk.Err("cannot decode Domain.Sol.Dydt\n%v", err)
	}
	err = dec.Decode(&o.Sol.D2ydt2)
	if err != nil {
		return chk.Err("cannot decode Domain.Sol.D2ydt2\n%v", err)
	}
	return
}

// SaveIvs saves elements's internal values to a file which name is set with tidx (time output index)
func (o Domain) SaveIvs(tidx int, verbose bool) (err error) {

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf, o.Sim.EncType)

	// elements that go to file
	enc.Encode(o.MyCids)

	// encode internal variables
	for _, e := range o.Elems {
		err = e.Encode(enc)
		if err != nil {
			return
		}
	}

	// save file
	fn := out_ele_path(o.Sim.DirOut, o.Sim.Key, o.Sim.EncType, tidx, o.Proc)
	return save_file(fn, &buf, verbose)
}

// ReadIvs reads elements's internal values from a file which name is set with tidx (time output index)
func (o *Domain) ReadIvs(dir, fnkey, enctype string, tidx, proc int) (err error) {

	// open file
	fn := out_ele_path(dir, fnkey, enctype, tidx, proc)
	fil, err := os.Open(fn)
	if err != nil {
		return
	}
	defer func() { err = fil.Close() }()

	// decoder
	dec := GetDecoder(fil, enctype)

	// elements that are in file
	err = dec.Decode(&o.MyCids)
	if err != nil {
		return chk.Err("cannot decode elements ids:\n%v", err)
	}

	// decode internal variables
	for _, cid := range o.MyCids {
		elem := o.Cid2elem[cid]
		if elem == nil {
			return chk.Err("cannot find element with cid=%d", cid)
		}
		err = elem.Decode(dec)
		if err != nil {
			return chk.Err("cannot decode element:\n%v", err)
		}
	}
	return
}

// Out performs output of Solution and Internal values to files
func (o *Domain) Save(tidx int, verbose bool) (err error) {
	err = o.SaveSol(tidx, verbose)
	if err != nil {
		return
	}
	return o.SaveIvs(tidx, verbose)
}

// In performs the inverse operation of Out()
//
//  allInOne -- indicates that all results must be read into the root processor only
//              For example when plotting or generating VTU files (or testing)
//
//  If allInOne is false, each processor will read its part as described by Summary.
//  Thus, recoreving the state as in the previous simulation.
//
func (o *Domain) Read(sum *Summary, tidx, proc int, allInOne bool) (err error) {

	// serial run
	if allInOne {
		for i := 0; i < sum.Nproc; i++ {
			err = o.ReadIvs(sum.Dirout, sum.Fnkey, o.Sim.EncType, tidx, i)
			if err != nil {
				return
			}
		}
		return o.ReadSol(sum.Dirout, sum.Fnkey, o.Sim.EncType, tidx)
	}

	// parallel run
	err = o.ReadIvs(sum.Dirout, sum.Fnkey, o.Sim.EncType, tidx, proc)
	if err != nil {
		return
	}
	return o.ReadSol(sum.Dirout, sum.Fnkey, o.Sim.EncType, tidx)
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_nod_path(dir, fnkey, enctype string, tidx, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_nod_%010d.%s", fnkey, proc, tidx, enctype))
}

func out_ele_path(dir, fnkey, enctype string, tidx, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_ele_%010d.%s", fnkey, proc, tidx, enctype))
}

func save_file(filename string, buf *bytes.Buffer, verbose bool) (err error) {
	fil, err := os.Create(filename)
	if err != nil {
		return
	}
	defer func() { err = fil.Close() }()
	_, err = fil.Write(buf.Bytes())
	if verbose {
		io.Pfblue2("file <%s> written\n", filename)
	}
	return
}

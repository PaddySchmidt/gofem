// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

// global variables
var (
	ndim   int         // space dimension
	verts  []*inp.Vert // all vertices
	cells  []*inp.Cell // all cells
	dirout string      // directory for output
	fnkey  string      // filename key
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data
	var mshfn string
	mshfn, fnkey = io.ArgToFilename(0, "data/d2-coarse", ".msh", true)
	io.Pf("\n%s\n", io.ArgsTable(
		"mesh filename", "mshfn", mshfn,
	))

	// read mesh
	msh, err := inp.ReadMsh("", mshfn, 0)
	if err != nil {
		io.PfRed("cannot read mesh:\n%v", err)
		return
	}
	ndim = msh.Ndim
	verts = msh.Verts
	cells = msh.Cells
	dirout = "/tmp/gofem"

	// buffers
	geo := new(bytes.Buffer)
	vtu := new(bytes.Buffer)

	// generate topology
	topology(geo)

	// points data
	pdata_write(vtu)

	// cells data
	cdata_write(vtu)

	// write vtu file
	vtu_write(geo, vtu)
}

// headers and footers ///////////////////////////////////////////////////////////////////////////////

func vtu_write(geo, dat *bytes.Buffer) {
	if geo == nil || dat == nil {
		return
	}
	nv := len(verts)
	nc := len(cells)
	var hdr, foo bytes.Buffer
	io.Ff(&hdr, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n")
	io.Ff(&hdr, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nv, nc)
	io.Ff(&foo, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
	io.WriteFileVD(dirout, fnkey+".vtu", &hdr, geo, dat, &foo)
}

// topology ////////////////////////////////////////////////////////////////////////////////////////

func topology(buf *bytes.Buffer) {
	if buf == nil {
		return
	}

	// coordinates
	io.Ff(buf, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	var z float64
	for _, v := range verts {
		if ndim == 3 {
			z = v.C[2]
		}
		io.Ff(buf, "%23.15e %23.15e %23.15e ", v.C[0], v.C[1], z)
	}
	io.Ff(buf, "\n</DataArray>\n</Points>\n")

	// connectivities
	io.Ff(buf, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	for _, c := range cells {
		nverts, _ := c.GetVtkInfo(false)
		for j := 0; j < nverts; j++ {
			io.Ff(buf, "%d ", c.Verts[j])
		}
	}

	// offsets of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	for _, c := range cells {
		nverts, _ := c.GetVtkInfo(false)
		offset += nverts
		io.Ff(buf, "%d ", offset)
	}

	// types of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
	for _, c := range cells {
		_, vtkcode := c.GetVtkInfo(false)
		if vtkcode < 0 {
			chk.Panic("cannot handle cell type %q", c.Shp.Type)
		}
		io.Ff(buf, "%d ", vtkcode)
	}
	io.Ff(buf, "\n</DataArray>\n</Cells>\n")
	return
}

// points data /////////////////////////////////////////////////////////////////////////////////////

func pdata_write(buf *bytes.Buffer) {

	// open
	io.Ff(buf, "<PointData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Int32\" Name=\"nid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, v := range verts {
		io.Ff(buf, "%d ", v.Id)
	}

	// positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, v := range verts {
		io.Ff(buf, "%d ", iabs(v.Tag))
	}

	// close
	io.Ff(buf, "\n</DataArray>\n</PointData>\n")
}

func cdata_write(buf *bytes.Buffer) {

	// open
	io.Ff(buf, "<CellData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Int32\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, c := range cells {
		io.Ff(buf, "%d ", c.Id)
	}

	// cells positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, c := range cells {
		ptag := iabs(c.Tag)
		io.Ff(buf, "%d ", ptag)
	}

	// close
	io.Ff(buf, "\n</DataArray>\n</CellData>\n")
}

func iabs(val int) int {
	if val < 0 {
		return -val
	}
	return val
}

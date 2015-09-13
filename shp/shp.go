// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package shp implements shape structures/routines
package shp

import (
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// constants
const MINDET = 1.0e-14 // minimum determinant allowed for dxdR

// ShpFunc is the shape functions callback function
type ShpFunc func(S []float64, dSdR [][]float64, r []float64, derivs bool, idxface int)

// Shape holds geometry data
type Shape struct {

	// geometry
	Type           string      // name; e.g. "lin2"
	Func           ShpFunc     // shape/derivs function callback function
	FaceFunc       ShpFunc     // face shape/derivs function callback function
	BasicType      string      // geometry of basic element; e.g. "qua8" => "qua4"
	FaceType       string      // geometry of face; e.g. "qua8" => "lin3"
	Gndim          int         // geometry of shape; e.g. "lin3" => gnd == 1 (even in 3D simulations)
	Nverts         int         // number of vertices in cell; e.g. "qua8" => 8
	VtkCode        int         // VTK code
	FaceNvertsMax  int         // max number of vertices on face
	FaceLocalVerts [][]int     // face local vertices [nfaces][...]
	NatCoords      [][]float64 // natural coordinates [gndim][nverts]

	// geometry: for seams (3D-edges)
	SeamType       int     // geometry of seam (3D-edge); e.g. "hex8" => "lin2"
	SeamLocalVerts [][]int // seam (3d-edge) local vertices [nseams][nVertsOnSeam]

	// scratchpad: volume
	S    []float64   // [nverts] shape functions
	G    [][]float64 // [nverts][gndim] G == dSdx. derivative of shape function
	J    float64     // Jacobian: determinant of dxdr
	DSdR [][]float64 // [nverts][gndim] derivatives of S w.r.t natural coordinates
	DxdR [][]float64 // [gndim][gndim] derivatives of real coordinates w.r.t natural coordinates
	DRdx [][]float64 // [gndim][gndim] dRdx == inverse(dxdR)

	// scratchpad: line
	Jvec3d []float64 // Jacobian: norm of dxdr for line elements (size==3)
	Gvec   []float64 // [nverts] G == dSdx. derivative of shape function

	// scratchpad: face
	Sf     []float64   // [FaceNvertsMax] shape functions values
	Fnvec  []float64   // [gndim] face normal vector multiplied by Jf (and Ju if NURBS)
	DSfdRf [][]float64 // [FaceNvertsMax][gndim-1] derivatives of Sf w.r.t natural coordinates
	DxfdRf [][]float64 // [gndim][gndim-1] derivatives of real coordinates w.r.t natural coordinates

	// NURBS
	Nurbs      *gm.Nurbs   // pointer to NURBS structure => indicates that this shape strucutre is based on NURBS
	NurbsFaces []*gm.Nurbs // boundaries (surfaces) of NURBS [normalTo0, normalTo0, normalTo1, normalTo1, normalTo2, normalTo2]
	Span       []int       // NURBS knots' indices defining cell/element; e.g. [2, 3, 1, 2] for x-quad/y-lin cell
	Ibasis     []int       // indices of basis functions corresponding to Span == local indices of control points
	SpanFace   [][]int     // NURBS knots' indices defining cell/element; e.g. [2, 3, 1, 2] for x-quad/y-lin cell
	IbasisFace [][]int     // indices of basis functions corresponding to Span == local indices of control points
	U          []float64   // [gndim] NURBS' parametric space coordinates
	Ju         float64     // parametric-natural mapping Jacobian: determinant of dudr
}

// GetCopy returns a new copy of this shape structure
func (o Shape) GetCopy() *Shape {

	// new structure
	var p Shape

	// geometry
	p.Type = o.Type
	p.Func = o.Func
	p.FaceFunc = o.FaceFunc
	p.BasicType = o.BasicType
	p.FaceType = o.FaceType
	p.Gndim = o.Gndim
	p.Nverts = o.Nverts
	p.VtkCode = o.VtkCode
	p.FaceNvertsMax = o.FaceNvertsMax
	p.FaceLocalVerts = utl.IntsClone(o.FaceLocalVerts)
	p.NatCoords = la.MatClone(o.NatCoords)

	// geometry: for seams (3D-edges)
	p.SeamType = o.SeamType
	p.SeamLocalVerts = utl.IntsClone(o.SeamLocalVerts)

	// scratchpad: volume
	p.S = la.VecClone(o.S)
	p.G = la.MatClone(o.G)
	p.J = o.J
	p.DSdR = la.MatClone(o.DSdR)
	p.DxdR = la.MatClone(o.DxdR)
	p.DRdx = la.MatClone(o.DRdx)

	// scratchpad: line
	p.Jvec3d = la.VecClone(o.Jvec3d)
	p.Gvec = la.VecClone(o.Gvec)

	// scratchpad: face
	p.Sf = la.VecClone(o.Sf)
	p.Fnvec = la.VecClone(o.Fnvec)
	p.DSfdRf = la.MatClone(o.DSfdRf)
	p.DxfdRf = la.MatClone(o.DxfdRf)
	return &p
}

// factory holds all Shapes available
var factory = make(map[string]*Shape)

// Get returns an existent Shape structure
//  Note: 1) returns nil on errors
//        2) use goroutineId > 0 to get a copy
func Get(geoType string, goroutineId int) *Shape {
	s, ok := factory[geoType]
	if !ok {
		return nil
	}
	if goroutineId > 0 {
		return s.GetCopy()
	}
	return s
}

// IpRealCoords returns the real coordinates (y) of an integration point
func (o *Shape) IpRealCoords(x [][]float64, ip Ipoint) (y []float64) {
	ndim := len(x)
	y = make([]float64, ndim)
	o.Func(o.S, o.DSdR, ip, false, -1)
	for i := 0; i < ndim; i++ {
		for m := 0; m < o.Nverts; m++ {
			y[i] += o.S[m] * x[i][m]
		}
	}
	return
}

// FaceIpRealCoords returns the real coordinates (y) of an integration point @ face
// TODO: check this function
func (o *Shape) FaceIpRealCoords(x [][]float64, ipf Ipoint, idxface int) (y []float64) {
	ndim := len(x)
	y = make([]float64, ndim)
	o.FaceFunc(o.Sf, o.DSfdRf, ipf, false, idxface)
	for i := 0; i < ndim; i++ {
		for k, n := range o.FaceLocalVerts[idxface] {
			y[i] += o.Sf[k] * x[i][n]
		}
	}
	return
}

// CalcAtIp calculates volume data such as S and G at natural coordinate r
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   ip                -- integration point
//  Output:
//   S, DSdR, DxdR, DRdx, G, and J
func (o *Shape) CalcAtIp(x [][]float64, ip Ipoint, derivs bool) (err error) {

	// S and dSdR
	o.Func(o.S, o.DSdR, ip, derivs, -1)
	if !derivs {
		return
	}

	if o.Gndim == 1 {
		// calculate Jvec3d == dxdR
		for i := 0; i < len(x); i++ {
			o.Jvec3d[i] = 0.0
			for m := 0; m < o.Nverts; m++ {
				o.Jvec3d[i] += x[i][m] * o.DSdR[m][0] // dxdR := x * dSdR
			}
		}

		// calculate J = norm of Jvec3d
		o.J = la.VecNorm(o.Jvec3d)

		// calculate G
		for m := 0; m < o.Nverts; m++ {
			o.Gvec[m] = o.DSdR[m][0] / o.J
		}

		return
	}

	// dxdR := sum_n x * dSdR   =>  dx_i/dR_j := sum_n x^n_i * dS^n/dR_j
	for i := 0; i < len(x); i++ {
		for j := 0; j < o.Gndim; j++ {
			o.DxdR[i][j] = 0.0
			for n := 0; n < o.Nverts; n++ {
				o.DxdR[i][j] += x[i][n] * o.DSdR[n][j]
			}
		}
	}

	// dRdx := inv(dxdR)
	o.J, err = la.MatInv(o.DRdx, o.DxdR, MINDET)
	if err != nil {
		return
	}

	// fix J if NURBS
	if o.Nurbs != nil {
		o.J *= o.Ju
	}

	// G == dSdx := dSdR * dRdx  =>  dS^m/dR_i := sum_i dS^m/dR_i * dR_i/dx_j
	la.MatMul(o.G, 1, o.DSdR, o.DRdx)
	return
}

// CalcAtR calculates volume data such as S and G at natural coordinate r
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   R[3]              -- local/natural coordinates
//  Output:
//   S, DSdR, DxdR, DRdx, G, and J
func (o *Shape) CalcAtR(x [][]float64, R []float64, derivs bool) (err error) {
	return o.CalcAtIp(x, R, derivs)
}

// CalcAtFaceIp calculates face data such as Sf and Fnvec
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   ipf               -- local/natural coordinates of face
//   idxface           -- local index of face
//  Output:
//   Sf and Fnvec
func (o *Shape) CalcAtFaceIp(x [][]float64, ipf Ipoint, idxface int) (err error) {

	// skip 1D elements
	if o.Gndim == 1 {
		return
	}

	// Sf and dSfdR
	o.Ju = 1.0
	o.FaceFunc(o.Sf, o.DSfdRf, ipf, true, idxface)

	// dxfdRf := sum_n x * dSfdRf   =>  dxf_i/dRf_j := sum_n xf^n_i * dSf^n/dRf_j
	for i := 0; i < len(x); i++ {
		for j := 0; j < o.Gndim-1; j++ {
			o.DxfdRf[i][j] = 0.0
			for k, n := range o.FaceLocalVerts[idxface] {
				o.DxfdRf[i][j] += x[i][n] * o.DSfdRf[k][j]
			}
		}
	}

	// face normal vector
	if o.Gndim == 2 {
		o.Fnvec[0] = o.DxfdRf[1][0] * o.Ju
		o.Fnvec[1] = -o.DxfdRf[0][0] * o.Ju
		return
	}
	o.Fnvec[0] = (o.DxfdRf[1][0]*o.DxfdRf[2][1] - o.DxfdRf[2][0]*o.DxfdRf[1][1]) * o.Ju
	o.Fnvec[1] = (o.DxfdRf[2][0]*o.DxfdRf[0][1] - o.DxfdRf[0][0]*o.DxfdRf[2][1]) * o.Ju
	o.Fnvec[2] = (o.DxfdRf[0][0]*o.DxfdRf[1][1] - o.DxfdRf[1][0]*o.DxfdRf[0][1]) * o.Ju
	return
}

// AxisymGetRadius returns the x0 == radius for axisymmetric computations
//  Note: must be called after CalcAtIp
func (o *Shape) AxisymGetRadius(x [][]float64) (radius float64) {
	for m := 0; m < o.Nverts; m++ {
		radius += o.S[m] * x[0][m]
	}
	return
}

// AxisymGetRadiusF (face) returns the x0 == radius for axisymmetric computations
//  Note: must be called after CalcAtFaceIp
func (o *Shape) AxisymGetRadiusF(x [][]float64, idxface int) (radius float64) {
	for m := 0; m < o.FaceNvertsMax; m++ {
		radius += o.Sf[m] * x[0][o.FaceLocalVerts[idxface][m]]
	}
	return
}

// init_scratchpad initialise volume data (scratchpad)
func (o *Shape) init_scratchpad() {

	// volume data
	o.S = make([]float64, o.Nverts)
	o.DSdR = la.MatAlloc(o.Nverts, o.Gndim)
	o.DxdR = la.MatAlloc(o.Gndim, o.Gndim)
	o.DRdx = la.MatAlloc(o.Gndim, o.Gndim)
	o.G = la.MatAlloc(o.Nverts, o.Gndim)

	// face data
	if o.Gndim > 1 {
		o.Sf = make([]float64, o.FaceNvertsMax)
		o.DSfdRf = la.MatAlloc(o.FaceNvertsMax, o.Gndim-1)
		o.DxfdRf = la.MatAlloc(o.Gndim, o.Gndim-1)
		o.Fnvec = make([]float64, o.Gndim)
	}

	// lin data
	if o.Gndim == 1 {
		o.Jvec3d = make([]float64, 3)
		o.Gvec = make([]float64, o.Nverts)
	}
}

// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package inp implements the input data read from a (.sim) JSON file
package inp

import (
	"encoding/json"
	"io"
	"log"
	"math"
	"os"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

// Data holds global data for simulations
type Data struct {

	// global information
	Desc    string `json:"desc"`    // description of simulation
	Matfile string `json:"matfile"` // materials file path
	DirOut  string `json:"dirout"`  // directory for output; e.g. /tmp/gofem
	Encoder string `json:"encoder"` // encoder name; e.g. "gob" "json" "xml"

	// problem definition and options
	Steady  bool `json:"steady"`  // steady simulation
	Pstress bool `json:"pstress"` // plane-stress
	Axisym  bool `json:"axisym"`  // axisymmetric

	// options
	React bool `json:"react"` // indicates whether or not reaction forces must be computed
	ShowR bool `json:"showr"` // show residual
	NoDiv bool `json:"nodiv"` // disregard divergence control in both fb or Lδu
	CteTg bool `json:"ctetg"` // use constant tangent (modified Newton) during iterations

	// derived
	FnameKey string // simulation filename key; e.g. mysim01.sim => mysim01
}

// SetDefault sets defaults values
func (o *Data) SetDefault() {
	o.DirOut = "/tmp/gofem"
	o.Encoder = "gob"
}

// PostProcess performs a post-processing of the just read json file
func (o *Data) PostProcess(simfilepath string, erasefiles bool) {
	if o.DirOut == "" {
		o.DirOut = "/tmp/gofem"
	}
	if o.Encoder == "" {
		o.Encoder = "gob"
	}
	o.FnameKey = utl.FnKey(simfilepath)
	os.MkdirAll(o.DirOut, 0777)
	if erasefiles {
		utl.RemoveAll(utl.Sf("%s/%s_*.gob", o.DirOut, o.FnameKey))
		utl.RemoveAll(utl.Sf("%s/%s_*.json", o.DirOut, o.FnameKey))
	}
}

// LinSolData holds data for linear solvers
type LinSolData struct {
	Name      string `json:"name"`      // "mumps" or "umfpack"
	Symmetric bool   `json:"symmetric"` // use symmetric solver
	Verbose   bool   `json:"verbose"`   // verbose?
	Timing    bool   `json:"timing"`    // show timing statistics
	Ordering  string `json:"ordering"`  // ordering scheme
	Scaling   string `json:"scaling"`   // scaling scheme
}

// SetDefault sets defaults values
func (o *LinSolData) SetDefault() {
	o.Name = "umfpack"
	o.Ordering = "amf"
	o.Scaling = "rcit"
}

// PostProcess performs a post-processing of the just read json file
func (o *LinSolData) PostProcess() {
	if mpi.IsOn() {
		if mpi.Size() > 1 {
			o.Name = "mumps"
		}
	} else {
		o.Name = "umfpack"
	}
}

// SolverData holds FEM solver data
type SolverData struct {

	// constants
	Eps float64 // smallest number satisfying 1.0 + ϵ > 1.0

	// nonlinear solver
	NmaxIt int     `json:"nmaxit"` // number of max iterations
	Atol   float64 `json:"atol"`   // absolute tolerance
	Rtol   float64 `json:"rtol"`   // relative tolerance
	FbTol  float64 `json:"fbtol"`  // tolerance for convergence on fb
	FbMin  float64 `json:"fbmin"`  // minimum value of fb

	// transient analyses
	DtMin float64 `json:"dtmin"` // minium value of Dt for transient (θ and Newmark / Dyn coefficients)
	Theta float64 `json:"theta"` // θ-method

	// dynamics
	Theta1 float64 `json:"theta1"` // Newmark's method parameter
	Theta2 float64 `json:"theta2"` // Newmark's method parameter
	HHT    bool    `json:"hht"`    // use Hilber-Hughes-Taylor method
	HHTalp float64 `json:"hhtalp"` // HHT α parameter
	RayM   float64 `json:"raym"`   // Rayleigh damping coefficient
	RayK   float64 `json:"rayk"`   // Rayleigh damping coefficient

	// derived
	Itol float64 // iterations tolerance
}

// SetDefault set defaults values
func (o *SolverData) SetDefault() {

	// constants
	o.Eps = 1e-16

	// nonlinear solver
	o.NmaxIt = 20
	o.Atol = 1e-6
	o.Rtol = 1e-6
	o.FbTol = 1e-8
	o.FbMin = 1e-14

	// transient analyses
	o.DtMin = 1e-8
	o.Theta = 0.5

	// dynamics
	o.Theta1 = 0.5
	o.Theta2 = 0.5
	o.HHTalp = 0.5
}

// PostProcess performs a post-processing of the just read json file
func (o *SolverData) PostProcess() {
	o.Itol = max(10.0*o.Eps/o.Rtol, min(0.01, math.Sqrt(o.Rtol)))
}

// ElemData holds element data
type ElemData struct {
	Tag   int    `json:"tag"`   // tag of element
	Mat   string `json:"mat"`   // material name
	Type  string `json:"type"`  // type of element. ex: u, p, up, rod, beam, rjoint
	Extra string `json:"extra"` // extra flags (in keycode format). ex: "!thick:0.2 !nip:4"
	Inact bool   `json:"inact"` // whether element starts inactive or not
}

// Region holds region data
type Region struct {

	// input data
	Desc      string      `json:"desc"`      // description of region. ex: ground, indenter, etc.
	Mshfile   string      `json:"mshfile"`   // file path of file with mesh data
	ElemsData []*ElemData `json:"elemsdata"` // list of elements data

	// derived data
	etag2idx map[int]int // maps element tag to element index in ElemsData slice
}

// FaceBc holds face boundary condition
type FaceBc struct {
	Tag   int      `json:"tag"`   // tag of face
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: qn, pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// SeamBc holds seam (3D edge) boundary condition
type SeamBc struct {
	Tag   int      `json:"tag"`   // tag of seam
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: qn, pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// NodeBc holds node boundary condition
type NodeBc struct {
	Tag   int      `json:"tag"`   // tag of node
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// EleCond holds element condition
type EleCond struct {
	Tag   int      `json:"tag"`   // tag of cell/element
	Keys  []string `json:"keys"`  // key indicating type of condition. ex: "g" (gravity), "qn" for beams, etc.
	Funcs []string `json:"funcs"` // name of function. ex: grav, none
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// TimeControl holds data for defining the simulation time stepping
type TimeControl struct {
	Tf     float64 `json:"tf"`     // final time
	Dt     float64 `json:"dt"`     // time step size (if constant)
	DtOut  float64 `json:"dtout"`  // time step size for output
	DtFcn  string  `json:"dtfcn"`  // time step size (function name)
	DtoFcn string  `json:"tdofcn"` // time step size for output (function name)

	// derived
	DtFunc  fun.Func // time step function
	DtoFunc fun.Func // output time step function
}

// Stage holds stage data
type Stage struct {

	// main
	Desc       string `json:"desc"`       // description of simulation stage. ex: activation of top layer
	Activate   []int  `json:"activate"`   // array of tags of elements to be activated
	Deactivate []int  `json:"deactivate"` // array of tags of elements to be deactivated
	Save       bool   `json:"save"`       // save stage data to binary file
	Load       string `json:"load"`       // load stage data (filename) from binary file

	// conditions
	EleConds []*EleCond `json:"eleconds"` // element conditions. ex: gravity or beam distributed loads
	FaceBcs  []*FaceBc  `json:"facebcs"`  // face boundary conditions
	SeamBcs  []*SeamBc  `json:"seambcs"`  // seam (3D) boundary conditions
	NodeBcs  []*NodeBc  `json:"nodebcs"`  // node boundary conditions

	// timecontrol
	Control TimeControl `json:"control"` // time control
}

// Simulation holds all simulation data
type Simulation struct {
	Data      Data       `json:"data"`      // stores global simulation data
	Functions FuncsData  `json:"functions"` // stores all boundary condition functions
	Regions   []*Region  `json:"regions"`   // stores all regions
	LinSol    LinSolData `json:"linsol"`    // linear solver data
	Solver    SolverData `json:"solver"`    // FEM solver data
	Stages    []*Stage   `json:"stages"`    // stores all stages
}

// ReadSim reads all simulation data from a .sim JSON file
func ReadSim(simfilepath string, erasefiles, dolog bool) (o *Simulation) {

	// new sim
	o = new(Simulation)

	// read file
	b, err := utl.ReadFile(simfilepath)
	if err != nil {
		utl.Panic("%v", err.Error())
	}

	// set default values
	o.Data.SetDefault()
	o.Solver.SetDefault()
	o.LinSol.SetDefault()

	// decode
	err = json.Unmarshal(b, o)
	if err != nil {
		utl.Panic("%v", err.Error())
	}

	// derived data
	o.Data.PostProcess(simfilepath, erasefiles)
	o.LinSol.PostProcess()
	o.Solver.PostProcess()

	// for all regions
	for _, reg := range o.Regions {

		// dependent variables
		reg.etag2idx = make(map[int]int)
		for j, ed := range reg.ElemsData {
			reg.etag2idx[ed.Tag] = j
		}
	}

	// for all stages
	var t float64
	for _, stg := range o.Stages {

		// fix Tf
		if stg.Control.Tf < 1e-14 {
			stg.Control.Tf = 1
		}

		// fix Dt
		if stg.Control.DtFcn == "" {
			if stg.Control.Dt < 1e-14 {
				stg.Control.Dt = 1
			}
			stg.Control.DtFunc = &fun.Cte{C: stg.Control.Dt}
		} else {
			stg.Control.DtFunc = o.Functions.GetOrPanic(stg.Control.DtFcn)
			stg.Control.Dt = stg.Control.DtFunc.F(t, nil)
		}

		// fix DtOut
		if stg.Control.DtoFcn == "" {
			if stg.Control.DtOut < 1e-14 {
				stg.Control.DtOut = stg.Control.Dt
				stg.Control.DtoFunc = stg.Control.DtFunc
			} else {
				if stg.Control.DtOut < stg.Control.Dt {
					stg.Control.DtOut = stg.Control.Dt
				}
				stg.Control.DtoFunc = &fun.Cte{C: stg.Control.DtOut}
			}
		} else {
			stg.Control.DtoFunc = o.Functions.GetOrPanic(stg.Control.DtoFcn)
			stg.Control.DtOut = stg.Control.DtoFunc.F(t, nil)
		}

		// update time
		t += stg.Control.Tf
	}

	// log
	InitLogFile(o.Data.DirOut, o.Data.FnameKey)
	if dolog {
		log.Printf("sim: fn=%s desc=%s nfunctions=%d nregions=%d nstages=%d linsol=%s itol=%g\n", simfilepath, o.Data.Desc, len(o.Functions), len(o.Regions), len(o.Stages), o.LinSol.Name, o.Solver.Itol)
	}
	return
}

// Etag2data returns the ElemData corresponding to element tag
func (d *Region) Etag2data(etag int) *ElemData {
	idx, ok := d.etag2idx[etag]
	if !ok {
		utl.Panic("cannot find element data with etag = %d", etag)
	}
	return d.ElemsData[idx]
}

// GetInfo returns formatted information
func (o *Simulation) GetInfo(w io.Writer) {
	b, err := json.MarshalIndent(o, "", "  ")
	if err != nil {
		utl.Panic("cannot marshal Simulation data")
	}
	w.Write(b)
}
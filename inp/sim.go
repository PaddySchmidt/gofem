// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package inp implements the input data read from a (.sim) JSON file
package inp

import (
	"encoding/json"
	goio "io"
	"log"
	"math"
	"os"
	"path/filepath"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
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
	Debug   bool `json:"debug"`   // activate debugging
	Stat    bool `json:"stat"`    // activate statistics

	// options
	React bool `json:"react"` // indicates whether or not reaction forces must be computed
	ShowR bool `json:"showr"` // show residual
	NoDiv bool `json:"nodiv"` // disregard divergence control in both fb or Lδu
	CteTg bool `json:"ctetg"` // use constant tangent (modified Newton) during iterations

	// derived
	FnameDir string // directory where .sim filename is locatd
	FnameKey string // simulation filename key; e.g. mysim01.sim => mysim01
}

// SetDefault sets defaults values
func (o *Data) SetDefault() {
	o.Encoder = "gob"
}

// PostProcess performs a post-processing of the just read json file
func (o *Data) PostProcess(dir, fn string, erasefiles bool) {
	o.FnameDir = os.ExpandEnv(dir)
	o.FnameKey = io.FnKey(fn)
	if o.DirOut == "" {
		o.DirOut = "/tmp/gofem/" + o.FnameKey
	}
	if o.Encoder != "gob" || o.Encoder != "json" {
		o.Encoder = "gob"
	}
	err := os.MkdirAll(o.DirOut, 0777)
	if err != nil {
		chk.Panic("cannot create directory for output results (%s): %v", o.DirOut, err)
	}
	if erasefiles {
		io.RemoveAll(io.Sf("%s/%s_*.gob", o.DirOut, o.FnameKey))
		io.RemoveAll(io.Sf("%s/%s_*.json", o.DirOut, o.FnameKey))
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

// HydroStData holds data to set hydrostatic initial state
type HydroStData struct {
	GamW   float64 `json:"gamw"`   // unit weight of water to use in HydroSt
	Zwater float64 `json:"zwater"` // water elevation to set ponding or unsaturated condition
	Unsat  bool    `json:"unsat"`  // consider unsaturated state with zwater < zmax
}

// IniStressData holds data for setting initial stresses
type IniStressData struct {
	Hom bool    `json:"hom"` // homogeneous stress distribution
	Iso bool    `json:"iso"` // isotropic state
	Psa bool    `json:"psa"` // plane-strain state
	S0  float64 `json:"s0"`  // Iso => stress value to use in homogeneous and isotropic distribution
	Sh  float64 `json:"sh"`  // Psa => horizontal stress
	Sv  float64 `json""sv"`  // Psa => vertical stress
	Nu  float64 `json:"nu"`  // Psa => Poisson's coefficient for plane-strain state
}

// GeoStData holds data for setting initial geostatic state (hydrostatic as well)
type GeoStData struct {
	HydroStData
	Hom   bool    `json:"hom"`   // homogeneous stress distribution
	Nu    float64 `json:"nu"`    // Poisson's coefficient to compute effective horizontal state
	K0    float64 `json:"K0"`    // Earth pressure coefficient at rest to compute effective horizontal stresses
	UseK0 bool    `json:"useK0"` // use K0 to compute effective horizontal stresses instead of "nu"
}

// Stage holds stage data
type Stage struct {

	// main
	Desc       string `json:"desc"`       // description of simulation stage. ex: activation of top layer
	Activate   []int  `json:"activate"`   // array of tags of elements to be activated
	Deactivate []int  `json:"deactivate"` // array of tags of elements to be deactivated
	Save       bool   `json:"save"`       // save stage data to binary file
	Load       string `json:"load"`       // load stage data (filename) from binary file

	// specific problems data
	HydroSt   *HydroStData   `json:"hydrost"`   // hydrostatic data
	SeepFaces []int          `json:"seepfaces"` // face tags corresponding to seepage faces
	IniStress *IniStressData `json:"inistress"` // initial stress data
	GeoSt     *GeoStData     `json:"geost"`     // initial geostatic state data (hydrostatic as well)

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
//  Notes:  1) this function initialises log file
//          2) returns nil on errors
func ReadSim(dir, fn string, erasefiles bool) *Simulation {

	// new sim
	var o Simulation

	// read file
	b, err := io.ReadFile(filepath.Join(dir, fn))
	if err != nil {
		io.PfRed("sim: cannot read simulation file %s/%s\n%v\n", dir, fn, err)
		return nil
	}

	// set default values
	o.Data.SetDefault()
	o.Solver.SetDefault()
	o.LinSol.SetDefault()

	// decode
	err = json.Unmarshal(b, &o)
	if err != nil {
		io.PfRed("sim: cannot unmarshal simulation file %s/%s\n%v\n", dir, fn, err)
		return nil
	}

	// derived data
	o.Data.PostProcess(dir, fn, erasefiles)
	o.LinSol.PostProcess()
	o.Solver.PostProcess()

	// init log file
	err = InitLogFile(o.Data.DirOut, o.Data.FnameKey)
	if err != nil {
		io.PfRed("sim: cannot create log file\n%v", err)
		return nil
	}

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
			stg.Control.DtFunc = o.Functions.Get(stg.Control.DtFcn)
			if LogErrCond(stg.Control.DtFunc == nil, "sim: cannot get Dt function named %s\n", stg.Control.DtFcn) {
				return nil
			}
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
			stg.Control.DtoFunc = o.Functions.Get(stg.Control.DtoFcn)
			if LogErrCond(stg.Control.DtoFunc == nil, "sim: cannot get DtOut function named %s\n", stg.Control.DtoFcn) {
				return nil
			}
			stg.Control.DtOut = stg.Control.DtoFunc.F(t, nil)
		}

		// update time
		t += stg.Control.Tf
	}

	// log
	log.Printf("sim: file=%s/%s desc=%q nfunctions=%d nregions=%d nstages=%d linsol=%s itol=%g\n", dir, fn, o.Data.Desc, len(o.Functions), len(o.Regions), len(o.Stages), o.LinSol.Name, o.Solver.Itol)
	return &o
}

// Etag2data returns the ElemData corresponding to element tag
//  Note: returns nil if not found
func (d *Region) Etag2data(etag int) *ElemData {
	idx, ok := d.etag2idx[etag]
	if !ok {
		return nil
	}
	return d.ElemsData[idx]
}

// GetInfo returns formatted information
func (o *Simulation) GetInfo(w goio.Writer) (err error) {
	b, err := json.MarshalIndent(o, "", "  ")
	if err != nil {
		return err
	}
	_, err = w.Write(b)
	return
}

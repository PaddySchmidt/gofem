// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package inp implements the input data read from a (.sim) JSON file
package inp

import (
	"encoding/json"
	goio "io"
	"math"
	"os"
	"path/filepath"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
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
	Steady  bool    `json:"steady"`  // steady simulation
	Axisym  bool    `json:"axisym"`  // axisymmetric
	Pstress bool    `json:"pstress"` // plane-stress
	NoLBB   bool    `json:"nolbb"`   // do not satisfy Ladyženskaja-Babuška-Brezzi condition; i.e. do not use [qua8,qua4] for u-p formulation
	Debug   bool    `json:"debug"`   // activate debugging
	Stat    bool    `json:"stat"`    // activate statistics
	Wlevel  float64 `json:"wlevel"`  // water level; 0 means use max elevation
	Surch   float64 `json:"surch"`   // surcharge load at surface == qn0
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

// SolverData holds FEM solver data
type SolverData struct {

	// nonlinear solver
	Type    string  `json:"type"`    // nonlinear solver type: {imp, exp, rex} => implicit, explicit, Richardson extrapolation
	NmaxIt  int     `json:"nmaxit"`  // number of max iterations
	Atol    float64 `json:"atol"`    // absolute tolerance
	Rtol    float64 `json:"rtol"`    // relative tolerance
	FbTol   float64 `json:"fbtol"`   // tolerance for convergence on fb
	FbMin   float64 `json:"fbmin"`   // minimum value of fb
	DvgCtrl bool    `json:"dvgctrl"` // use divergence control
	NdvgMax int     `json:"ndvgmax"` // max number of continued divergence
	CteTg   bool    `json:"ctetg"`   // use constant tangent (modified Newton) during iterations
	ShowR   bool    `json:"showr"`   // show residual

	// Richardson's extrapolation
	REnogus  bool    `json:"renogus"`  // Richardson extrapolation: no Gustafsson's step control
	REnssmax int     `json:"renssmax"` // Richardson extrapolation: max number of substeps
	REatol   float64 `json:"reatol"`   // Richardson extrapolation: absolute tolerance
	RErtol   float64 `json:"rertol"`   // Richardson extrapolation: relative tolerance
	REmfac   float64 `json:"remfac"`   // Richardson extrapolation: multiplier factor
	REmmin   float64 `json:"remmin"`   // Richardson extrapolation: min multiplier
	REmmax   float64 `json:"remmax"`   // Richardson extrapolation: max multiplier

	// transient analyses
	DtMin      float64 `json:"dtmin"`      // minium value of Dt for transient (θ and Newmark / Dyn coefficients)
	Theta      float64 `json:"theta"`      // θ-method
	ThGalerkin bool    `json:"thgalerkin"` // use θ = 2/3
	ThLiniger  bool    `json:"thliniger"`  // use θ = 0.878

	// dynamics
	Theta1 float64 `json:"theta1"` // Newmark's method parameter
	Theta2 float64 `json:"theta2"` // Newmark's method parameter
	HHT    bool    `json:"hht"`    // use Hilber-Hughes-Taylor method
	HHTalp float64 `json:"hhtalp"` // HHT α parameter

	// combination of coefficients
	ThCombo1 bool `json:"thcombo1"` // use θ=2/3, θ1=5/6 and θ2=8/9 to avoid oscillations

	// constants
	Eps float64 `json:"eps"` // smallest number satisfying 1.0 + ϵ > 1.0

	// derived
	Itol float64 // iterations tolerance
}

// ElemData holds element data
type ElemData struct {

	// input data
	Tag   int    `json:"tag"`   // tag of element
	Mat   string `json:"mat"`   // material name
	Type  string `json:"type"`  // type of element. ex: u, p, up, rod, beam, rjoint
	Nip   int    `json:"nip"`   // number of integration points; 0 => use default
	Nipf  int    `json:"nipf"`  // number of integration points on face; 0 => use default
	Extra string `json:"extra"` // extra flags (in keycode format). ex: "!thick:0.2 !nip:4"
	Inact bool   `json:"inact"` // whether element starts inactive or not

	// derived
	Lbb bool // LBB element; e.g. if "up", "upp", etc., unless NoLBB is true
}

// Region holds region data
type Region struct {

	// input data
	Desc      string      `json:"desc"`      // description of region. ex: ground, indenter, etc.
	Mshfile   string      `json:"mshfile"`   // file path of file with mesh data
	ElemsData []*ElemData `json:"elemsdata"` // list of elements data
	AbsPath   bool        `json:"abspath"`   // mesh filename is given in absolute path

	// derived
	Msh      *Mesh       // the mesh
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
	DtoFcn string  `json:"dtofcn"` // time step size for output (function name)

	// derived
	DtFunc  fun.Func // time step function
	DtoFunc fun.Func // output time step function
}

// GeoStData holds data for setting initial geostatic state (hydrostatic as well)
type GeoStData struct {
	Nu     []float64 `json:"nu"`     // [nlayers] Poisson's coefficient to compute effective horizontal state for each layer
	K0     []float64 `json:"K0"`     // [nlayers] Earth pressure coefficient at rest to compute effective horizontal stresses
	UseK0  []bool    `json:"useK0"`  // [nlayers] use K0 to compute effective horizontal stresses instead of "nu"
	Layers [][]int   `json:"layers"` // [nlayers][ntagsInLayer]; e.g. [[-1,-2], [-3,-4]] => 2 layers
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

// InitialData holds data for setting initial solution values such as Y, dYdt and d2Ydt2
type InitialData struct {
	File string   `json:"file"` // file with values at each node is given; filename with path is provided
	Fcns []string `json:"fcns"` // functions F(t, x) are given; from functions database
	Dofs []string `json:"dofs"` // degrees of freedom corresponding to "fcns"
}

// ImportRes holds definitions for importing results from a previous simulation
type ImportRes struct {
	Dir    string `json:"dir"`    // output directory with previous simulation files
	Fnk    string `json:"fnk"`    // previous simulation file name key (without .sim)
	ResetU bool   `json:"resetu"` // reset/zero u (displacements)
}

// Stage holds stage data
type Stage struct {

	// main
	Desc       string `json:"desc"`       // description of simulation stage. ex: activation of top layer
	Activate   []int  `json:"activate"`   // array of tags of elements to be activated
	Deactivate []int  `json:"deactivate"` // array of tags of elements to be deactivated
	Save       bool   `json:"save"`       // save stage data to binary file
	Load       string `json:"load"`       // load stage data (filename) from binary file
	Skip       bool   `json:"skip"`       // do not run stage

	// specific problems data
	HydroSt   bool           `json:"hydrost"`   // hydrostatic initial condition
	SeepFaces []int          `json:"seepfaces"` // face tags corresponding to seepage faces
	IniStress *IniStressData `json:"inistress"` // initial stress data
	GeoSt     *GeoStData     `json:"geost"`     // initial geostatic state data (hydrostatic as well)
	Import    *ImportRes     `json:"import"`    // import results from another previous simulation
	Initial   *InitialData   `json:"initial"`   // set initial solution values such as Y, dYdt and d2Ydt2

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

	// input
	Data      Data       `json:"data"`      // stores global simulation data
	Functions FuncsData  `json:"functions"` // stores all boundary condition functions
	PlotF     *PlotFdata `json:"plotf"`     // plot functions
	Regions   []*Region  `json:"regions"`   // stores all regions
	LinSol    LinSolData `json:"linsol"`    // linear solver data
	Solver    SolverData `json:"solver"`    // FEM solver data
	Stages    []*Stage   `json:"stages"`    // stores all stages

	// derived
	GoroutineId int      // id of goroutine to avoid race problems
	DirOut      string   // directory to save results
	Key         string   // simulation key; e.g. mysim01.sim => mysim01 or mysim01-alias
	EncType     string   // encoder type
	MatParams   *MatDb   // materials' parameters
	Ndim        int      // space dimension
	MaxElev     float64  // maximum elevation
	Gravity     fun.Func // first stage: gravity constant function
	WaterRho0   float64  // first stage: intrinsic density of water corresponding to pressure pl=0
	WaterBulk   float64  // first stage: bulk modulus of water
	WaterLevel  float64  // first stage: water level == max(Wlevel, MaxElev)
}

// Simulation //////////////////////////////////////////////////////////////////////////////////////

// ReadSim reads all simulation data from a .sim JSON file
func ReadSim(simfilepath, alias string, erasefiles bool, goroutineId int) *Simulation {

	// new sim
	var o Simulation
	o.GoroutineId = goroutineId

	// read file
	b, err := io.ReadFile(simfilepath)
	if err != nil {
		chk.Panic("ReadSim: cannot read simulation file %q", simfilepath)
	}

	// set default values
	o.Solver.SetDefault()
	o.LinSol.SetDefault()

	// decode
	err = json.Unmarshal(b, &o)
	if err != nil {
		chk.Panic("ReadSim: cannot unmarshal simulation file %q", simfilepath)
	}

	// input directory and filename key
	dir := filepath.Dir(simfilepath)
	fn := filepath.Base(simfilepath)
	dir = os.ExpandEnv(dir)
	fnkey := io.FnKey(fn)
	o.Key = fnkey
	if alias != "" {
		o.Key += "-" + alias
	}

	// output directory
	o.DirOut = o.Data.DirOut
	if o.DirOut == "" {
		o.DirOut = "/tmp/gofem/" + fnkey
	}

	// encoder type
	o.EncType = o.Data.Encoder
	if o.EncType != "gob" && o.EncType != "json" {
		o.EncType = "gob"
	}

	// create directory and erase previous simulation results
	if erasefiles {
		err = os.MkdirAll(o.DirOut, 0777)
		if err != nil {
			chk.Panic("cannot create directory for output results (%s): %v", o.DirOut, err)
		}
		io.RemoveAll(io.Sf("%s/%s*", o.DirOut, fnkey))
	}

	// set solver constants
	o.Solver.PostProcess()

	// read materials database
	o.MatParams = ReadMat(dir, o.Data.Matfile)
	if o.MatParams == nil {
		chk.Panic("ReadSim: cannot read materials database\n")
	}

	// for all regions
	for i, reg := range o.Regions {

		// read mesh
		ddir := dir
		if reg.AbsPath {
			ddir = ""
		}
		reg.Msh, err = ReadMsh(ddir, reg.Mshfile, goroutineId)
		if err != nil {
			chk.Panic("ReadSim: cannot read mesh file:\n%v", err)
		}

		// dependent variables
		reg.etag2idx = make(map[int]int)
		for j, ed := range reg.ElemsData {
			reg.etag2idx[ed.Tag] = j
		}

		// get ndim and max elevation
		if i == 0 {
			o.Ndim = reg.Msh.Ndim
			o.MaxElev = reg.Msh.Ymax
			if o.Ndim == 3 {
				o.MaxElev = reg.Msh.Zmax
			}
		} else {
			if reg.Msh.Ndim != o.Ndim {
				chk.Panic("ReadSim: Ndim value is inconsistent: %d != %d", reg.Msh.Ndim, o.Ndim)
			}
			if o.Ndim == 2 {
				o.MaxElev = utl.Max(o.MaxElev, reg.Msh.Ymax)
			} else {
				o.MaxElev = utl.Max(o.MaxElev, reg.Msh.Zmax)
			}
		}

		// get water data
		for _, mat := range o.MatParams.Materials {
			if mat.Model == "porous" {
				if prm := mat.Prms.Find("RhoL0"); prm != nil {
					o.WaterRho0 = prm.V
				}
				if prm := mat.Prms.Find("BulkL"); prm != nil {
					o.WaterBulk = prm.V
				}
			}
		}

		// set LBB flag
		if !o.Data.NoLBB {
			for _, ed := range reg.ElemsData {
				if ed.Type == "up" {
					ed.Lbb = true
				}
			}
		}
	}

	// water level
	o.WaterLevel = utl.Max(o.Data.Wlevel, o.MaxElev)

	// for all stages
	var t float64
	for i, stg := range o.Stages {

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
			if stg.Control.DtFunc == nil {
				chk.Panic("ReadSim: cannot find DtFunc named %q", stg.Control.DtFcn)
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
			if stg.Control.DtoFunc == nil {
				chk.Panic("ReadSim: cannot find DtoFunc named %q", stg.Control.DtoFcn)
			}
			stg.Control.DtOut = stg.Control.DtoFunc.F(t, nil)
		}

		// first stage
		if i == 0 {

			// gravity
			for _, econd := range stg.EleConds {
				for j, key := range econd.Keys {
					if key == "g" {
						if o.Gravity == nil {
							o.Gravity = o.Functions.Get(econd.Funcs[j])
							if o.Gravity == nil {
								chk.Panic("ReadSim: cannot find function named %q", econd.Funcs[j])
							}
							break
						}
					}
				}
			}
			if o.Gravity == nil {
				o.Gravity = &fun.Cte{C: 10}
			}
		}

		// update time
		t += stg.Control.Tf
	}

	// results
	return &o
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

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

// GetEleCond returns element condition structure by giving an elem tag
//  Note: returns nil if not found
func (o Stage) GetEleCond(elemtag int) *EleCond {
	for _, ec := range o.EleConds {
		if elemtag == ec.Tag {
			return ec
		}
	}
	return nil
}

// GetNodeBc returns node boundary condition structure by giving a node tag
//  Note: returns nil if not found
func (o Stage) GetNodeBc(nodetag int) *NodeBc {
	for _, nbc := range o.NodeBcs {
		if nodetag == nbc.Tag {
			return nbc
		}
	}
	return nil
}

// GetFaceBc returns face boundary condition structure by giving a face tag
//  Note: returns nil if not found
func (o Stage) GetFaceBc(facetag int) *FaceBc {
	for _, fbc := range o.FaceBcs {
		if facetag == fbc.Tag {
			return fbc
		}
	}
	return nil
}

// extra settings //////////////////////////////////////////////////////////////////////////////////

// SetDefault sets defaults values
func (o *LinSolData) SetDefault() {
	o.Name = "umfpack"
	o.Ordering = "amf"
	o.Scaling = "rcit"
}

// SetDefault set defaults values
func (o *SolverData) SetDefault() {

	// nonlinear solver
	o.Type = "imp"
	o.NmaxIt = 20
	o.Atol = 1e-6
	o.Rtol = 1e-6
	o.FbTol = 1e-8
	o.FbMin = 1e-14
	o.NdvgMax = 20

	// Richardson's extrapolation
	o.REnssmax = 10000
	o.REatol = 1e-6
	o.RErtol = 1e-6
	o.REmfac = 0.9
	o.REmmin = 0.1
	o.REmmax = 2.0

	// transient analyses
	o.DtMin = 1e-8
	o.Theta = 0.5

	// dynamics
	o.Theta1 = 0.5
	o.Theta2 = 0.5
	o.HHTalp = 0.5

	// constants
	o.Eps = 1e-16
}

// PostProcess performs a post-processing of the just read json file
func (o *SolverData) PostProcess() {

	// coefficients for transient analyses
	if o.ThGalerkin {
		o.Theta = 2.0 / 3.0
	}
	if o.ThLiniger {
		o.Theta = 0.878
	}
	if o.ThCombo1 {
		o.Theta = 2.0 / 3.0
		o.Theta1 = 5.0 / 6.0
		o.Theta2 = 8.0 / 9.0
	}

	// iterations tolerance
	o.Itol = utl.Max(10.0*o.Eps/o.Rtol, utl.Min(0.01, math.Sqrt(o.Rtol)))
}

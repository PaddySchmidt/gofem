// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/utl"
)

// constants
var (
	PlotSet1 = []string{"ed,q", "p,q", "ed,ev", "p,ev"}
	PlotSet2 = []string{"ed,q/p", "p,q", "ed,ev", "p,ev"}
	PlotSet3 = []string{"ed,q/p", "p,q", "ed,ev", "log(p),ev"}
	PlotSet4 = []string{"ed,q", "i,f", "p,q,ys", "ed,ev", "p,ev", "oct"}
	PlotSet5 = []string{"ed,q", "i,f", "p,q,ys", "ed,ev", "p,ev", "i,alp"}
	PlotSet6 = []string{"ed,q", "ed,ev", "p,q,ys", "p,ev", "s3,s1,ys", "oct,ys"}
	PlotSet7 = []string{"ed,q", "i,f", "p,q,ys", "ed,ev", "p,ev", "s3,s1,ys", "i,alp", "Dgam,f", "oct,ys"}
	PlotSet8 = []string{"ed,q", "i,f", "p,q,ys", "ed,ev", "log(p),ev", "s3,s1,ys", "i,alp", "Dgam,f", "oct,ys"}
)

type RampFcn_t func(x float64) float64
type PlotFcn_t func()

type Plotter struct {

	// optional variables
	PlotFcn  PlotFcn_t // callback function to call before saving plot
	PngRes   int       // resolution for .png files
	Split    bool      // split graphs instead of using subplot
	NptsPq   int       // number of points to draw yield surface (in 2D/p-q)
	NptsOct  int       // number of points to draw yield surface (in 2D/octahedral)
	NptsSig  int       // number of points to draw yield surface (in 2D/s3-s1)
	NptsRmp  int       // number of poitns to plot ramp function (contour)
	ArgsYs   string    // extra argumetns to plot yield surface
	Rmpf     RampFcn_t // ramp function
	SaveDir  string    // directory to put figure
	SaveFnk  string    // save figure after plot (filename)
	UseEps   bool      // save eps figure instead of png
	QdivP    bool      // q/p in ed, q plot?
	LogP     bool      // use log(p) in εv,log(p) plot
	WithYs   bool      // with yield surface (if model != nil)
	NoAlp    bool      // do not plot alphas in εv,log(p) plot
	Multq    bool      // multply q and εd according to Lode angle (for extension)
	YsRangeM float64   // multiplier to increase and decrease range of x,y values when drawing yield surfaces
	ArrWid   int       // width of arrows for predictor-corrector graph
	ClrPC    string    // color for predictor-corrector arrows
	Clr      string    // curve color
	Mrk      string    // curve marker
	Lbl      string    // curve label
	Ls       string    // curve linestyle
	LsAlt    string    // curve linestyle (alternative, e.g. graph with two y scales)
	SpMrk    string    // start-point marker
	EpMrk    string    // end-point marker
	SpClr    string    // start-point marker color
	EpClr    string    // end-point marker color
	SpMs     int       // start-point marker size
	EpMs     int       // end-point marker size
	YsClr0   string    // color for yield surface line (inner)
	YsClr1   string    // color for yield surface line (outer)
	YsLs0    string    // yield surface line style (inner)
	YsLs1    string    // yield surface line style (outer)
	YsLw0    float64   // yield surface line width (inner)
	YsLw1    float64   // yield surface line width (outer)
	Hspace   float64   // subplot horizontal spacing between rows
	Vspace   float64   // subplot vertical spacing between columns
	AxLblX   string    // axis y-label. "" => use default
	AxLblY   string    // axis x-label. "" => use default

	// SMP coefficients
	SMPa  float64 // SMP coefficient
	SMPb  float64 // SMP coefficient
	SMPβ  float64 // SMP coefficient
	SMPϵ  float64 // SMP coefficient
	SMPon bool    // SMP plotting is on

	// subplots
	Nrow int // subplot number of rows
	Ncol int // subplot number of cols
	Pidx int // subplot index

	// internal variables
	m     EPmodel // the model
	nalp  int     // number of α
	nsurf int     // number of yield surfaces

	// for computation of x(p)
	Pt float64 // tensile p
	Pr float64 // reference p

	// results
	P, Q, W []float64 // stress invariants
	Ev, Ed  []float64 // strain invariants

	// rosette
	maxR float64 // max radius for octrahedral rosette
	Phi  float64 // φ for plotting reference criteria in rosette

	// predictor-corrector
	PreCor [][]float64 // [npath][neps] predictor-corrector states

	// limits
	Pmin     float64              // min p value to use when drawing yield surfaces
	Pmax     float64              // max p value to use when drawing yield surfaces
	UsePmin  bool                 // use Pmin in yield surfaces drawing
	UsePmax  bool                 // use Pmax in yield surfaces drawing
	Lims     map[string][]float64 // limits to be used with a particular plotset; if not nil => to set plot area
	UseOct   bool                 // use octahedral invariants (poct,qoct) in p,q plot
	OctAxOff bool                 // turn off axes in ocahedral plane
	OctLims  []float64            // octahedral limits if not nil (to compute contour; not to set plot area)
	PqLims   []float64            // p-q limits if not nil (to compute contour; not to set plot area)
	S3s1Lims []float64            // s1-s3 limits if not nil (to compute contour; not to set plot area)
}

// SetFig sets figure space for plotting
// Note: this method is optional
func (o *Plotter) SetFig(split, epsfig bool, prop, width float64, savedir, savefnk string) {
	plt.Reset()
	if o.PngRes < 150 {
		o.PngRes = 150
	}
	o.Split = split
	o.UseEps = epsfig
	if o.UseEps {
		plt.SetForEps(prop, width)
	} else {
		plt.SetForPng(prop, width, o.PngRes)
	}
	o.SaveDir = savedir
	o.SaveFnk = io.FnKey(savefnk)
	o.maxR = -1
	// colors and markers
	o.set_default_clr_mrk()
}

// SetModel sets a solid model in Plotter
// Note: this method is optional
func (o *Plotter) SetModel(m EPmodel) {
	o.m = m
	o.nalp, o.nsurf = o.m.Info()
}

// Title addes title to plot
func (o *Plotter) Title(text string) {
	plt.SupTitle(text, "size=10")
}

// Plot runs the plot generation (basic set)
func (o *Plotter) Plot(keys []string, res []*State, sts [][]float64, first, last bool) {

	// auxiliary variables
	nr := imax(len(res), len(sts))
	if nr < 1 {
		return
	}
	x := make([]float64, nr)
	y := make([]float64, nr)
	o.P = make([]float64, nr)
	o.Q = make([]float64, nr)
	o.W = make([]float64, nr)
	o.Ev = make([]float64, nr)
	o.Ed = make([]float64, nr)

	// compute invariants
	for i := 0; i < len(res); i++ {
		if len(res[i].Sig) < 4 {
			chk.Panic("number of stress components is incorrect: %d", len(res[i].Sig))
		}
		o.P[i], o.Q[i], o.W[i] = tsr.M_pqw(res[i].Sig)
	}
	nsig := len(res[0].Sig)
	devε := make([]float64, nsig)
	for i := 0; i < len(sts); i++ {
		if len(sts[i]) < 4 {
			chk.Panic("number of strain components is incorrect: %d", len(sts[i]))
		}
		_, o.Ev[i], o.Ed[i] = tsr.M_devε(devε, sts[i])
	}

	// clear previous figure
	if first {
		plt.Clf()
		plt.SplotGap(0.35, 0.35)
		if o.Hspace > 0 {
			plt.SetHspace(o.Hspace)
		}
		if o.Vspace > 0 {
			plt.SetVspace(o.Vspace)
		}
	}

	// number of points for contour
	if o.NptsPq < 2 {
		o.NptsPq = 61
	}
	if o.NptsOct < 2 {
		o.NptsOct = 41
	}
	if o.NptsSig < 2 {
		o.NptsSig = 41
	}

	// subplot variables
	o.Pidx = 1
	o.Ncol, o.Nrow = utl.BestSquare(len(keys))
	if len(keys) == 2 {
		o.Ncol, o.Nrow = 1, 2
	}
	if len(keys) == 3 {
		o.Ncol, o.Nrow = 1, 3
	}

	// do plot
	for _, key := range keys {
		o.Subplot()
		switch key {
		case "ed,q":
			o.QdivP = false
			o.Plot_ed_q(x, y, res, sts, last)
		case "ed,q/p":
			o.QdivP = true
			o.Plot_ed_q(x, y, res, sts, last)
		case "p,q":
			o.WithYs = false
			o.Plot_p_q(x, y, res, sts, last)
		case "p,q,ys":
			o.WithYs = true
			o.Plot_p_q(x, y, res, sts, last)
		case "ed,ev":
			o.Plot_ed_ev(x, y, res, sts, last)
		case "p,ev":
			o.LogP = false
			o.Plot_p_ev(x, y, res, sts, last)
		case "log(p),ev":
			o.LogP = true
			o.Plot_p_ev(x, y, res, sts, last)
		case "i,f":
			o.Plot_i_f(x, y, res, sts, last)
		case "i,alp":
			o.Plot_i_alp(x, y, res, sts, last)
		case "Dgam,f":
			o.Plot_Dgam_f(x, y, res, sts, last)
		case "oct":
			o.WithYs = false
			o.Plot_oct(x, y, res, sts, last)
		case "oct,ys":
			o.WithYs = true
			o.Plot_oct(x, y, res, sts, last)
		case "s3,s1":
			o.WithYs = false
			o.Plot_s3_s1(x, y, res, sts, last)
		case "s3,s1,ys":
			o.WithYs = true
			o.Plot_s3_s1(x, y, res, sts, last)
		case "empty":
			continue
		default:
			chk.Panic("cannot handle key=%q", key)
		}
		if o.Split && last {
			o.Save("_", key)
		}
	}

	// save figure
	if !o.Split && last {
		o.Save("", "")
	}
}

func (o *Plotter) Plot_ed_q(x, y []float64, res []*State, sts [][]float64, last bool) {
	nr := len(res)
	if len(sts) != nr {
		return
	}
	k := nr - 1
	for i := 0; i < nr; i++ {
		x[i] = o.Ed[i] * 100.0
		if o.QdivP {
			y[i] = o.Q[i] / o.P[i]
		} else {
			y[i] = o.Q[i]
		}
		if o.Multq {
			y[i] *= fun.Sign(o.W[i])
		}
	}
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	if last {
		ylbl := "$q$"
		if o.QdivP {
			ylbl = "$q/p$"
		}
		plt.Gll("$\\varepsilon_d\\;[\\%]$", ylbl, "leg_out=1, leg_ncol=4, leg_hlen=1.5")
		if lims, ok := o.Lims["ed,q"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_ed_ev(x, y []float64, res []*State, sts [][]float64, last bool) {
	nr := len(sts)
	k := nr - 1
	for i := 0; i < nr; i++ {
		x[i], y[i] = o.Ed[i]*100.0, o.Ev[i]*100.0
	}
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	if last {
		plt.Gll("$\\varepsilon_d\\;[\\%]$", "$\\varepsilon_v\\;[\\%]$", "leg_out=1, leg_ncol=4, leg_hlen=1.5")
		if lims, ok := o.Lims["ed,ev"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_p_ev(x, y []float64, res []*State, sts [][]float64, last bool) {
	nr := len(res)
	if len(sts) != nr {
		return
	}
	k := nr - 1
	var x0, x1 []float64
	if !o.NoAlp {
		x0, x1 = make([]float64, nr), make([]float64, nr)
	}
	withα := false
	if o.LogP {
		xmin := o.calc_x(o.P[0])
		xmax := xmin
		for i := 0; i < nr; i++ {
			x[i], y[i] = o.calc_x(o.P[i]), o.Ev[i]*100.0
			if !o.NoAlp && len(res[i].Alp) > 0 {
				withα = true
				x0[i] = o.calc_x(res[i].Alp[0])
				if o.nsurf > 1 {
					x1[i] = o.calc_x(res[i].Alp[1])
				}
			}
			xmin = min(xmin, x[i])
			xmax = max(xmax, x[i])
		}
	} else {
		xmin := o.P[0]
		xmax := xmin
		for i := 0; i < nr; i++ {
			x[i], y[i] = o.P[i], o.Ev[i]*100.0
			if !o.NoAlp && len(res[i].Alp) > 0 {
				withα = true
				x0[i] = res[i].Alp[0]
				if o.nsurf > 1 {
					x1[i] = res[i].Alp[1]
				}
			}
			xmin = min(xmin, x[i])
			xmax = max(xmax, x[i])
		}
	}
	if withα {
		plt.Plot(x0, y, io.Sf("'r-', ls='--', lw=3, clip_on=0, color='grey', label=r'%s'", o.Lbl+" $\\alpha_0$"))
		if o.nsurf > 1 {
			plt.Plot(x1, y, io.Sf("'r-', ls=':', lw=3, clip_on=0, color='grey', label=r'%s'", o.Lbl+" $\\alpha_1$"))
		}
	}
	plt.Plot(x, y, io.Sf("'r.', ls='-', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	if last {
		xlbl := "$p$"
		if o.LogP {
			xlbl = "$\\log{[1+(p+p_t)/p_r]}$"
		}
		plt.Gll(xlbl, "$\\varepsilon_v\\;[\\%]$", "leg_out=1, leg_ncol=4, leg_hlen=2")
		if lims, ok := o.Lims["p,ev"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_i_f(x, y []float64, res []*State, sts [][]float64, last bool) {
	if o.m == nil {
		o.set_empty()
		return
	}
	nr := len(res)
	var y2 []float64
	if o.nsurf > 1 {
		y2 = make([]float64, nr)
	}
	for i := 0; i < nr; i++ {
		ys := o.m.YieldFuncs(res[i])
		y[i] = ys[0]
		if o.nsurf > 1 {
			y2[i] = ys[1]
		}
		x[i] = float64(i)
	}
	lbl := "f " + o.Lbl
	plt.Plot(x, y, io.Sf("'r.', ls='-', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Clr, o.Mrk, lbl))
	if o.nsurf > 1 {
		lbl = "F " + o.Lbl
		plt.Plot(x, y2, io.Sf("'b+', ls=':', lw=2, clip_on=0, color='%s', marker='%s', label=r'%s'", o.Clr, o.Mrk, lbl))
	}
	if last {
		plt.Gll("$i$", "$f,\\;F$", "leg_out=1, leg_ncol=4, leg_hlen=2")
		if lims, ok := o.Lims["i,f"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_i_alp(x, y []float64, res []*State, sts [][]float64, last bool) {
	nr := len(res)
	nα := len(res[0].Alp)
	if nα == 0 {
		o.set_empty()
		return
	}
	yy := la.MatAlloc(nα, nr)
	for i := 0; i < nr; i++ {
		x[i] = float64(i)
		for j := 0; j < nα; j++ {
			yy[j][i] = res[i].Alp[j]
		}
	}
	for j := 0; j < nα; j++ {
		lbl := io.Sf("$\\alpha_%d$ "+o.Lbl, j)
		plt.Plot(x, yy[j], io.Sf("'r-', ls='-', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Clr, o.Mrk, lbl))
	}
	if last {
		plt.Gll("$i$", "$\\alpha_k$", "leg_out=1, leg_ncol=4, leg_hlen=2")
		if lims, ok := o.Lims["i,alp"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) set_empty() {
	plt.AxisOff()
}

func (o *Plotter) Plot_Dgam_f(x, y []float64, res []*State, sts [][]float64, last bool) {
	if o.m == nil {
		o.set_empty()
		return
	}
	nr := len(res)
	k := nr - 1
	ys := o.m.YieldFuncs(res[0])
	fc0 := ys[0]
	xmi, xma, ymi, yma := res[0].Dgam, res[0].Dgam, fc0, fc0
	for i := 0; i < nr; i++ {
		x[i] = res[i].Dgam
		ys = o.m.YieldFuncs(res[i])
		y[i] = ys[0]
		xmi = min(xmi, x[i])
		xma = max(xma, x[i])
		ymi = min(ymi, y[i])
		yma = max(yma, y[i])
	}
	//o.DrawRamp(xmi, xma, ymi, yma)
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	if last {
		plt.Gll("$\\Delta\\gamma$", "$f$", "")
		if lims, ok := o.Lims["Dgam,f"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_p_q(x, y []float64, res []*State, sts [][]float64, last bool) {
	// stress path
	nr := len(res)
	k := nr - 1
	var xmi, xma, ymi, yma float64
	for i := 0; i < nr; i++ {
		x[i], y[i] = o.P[i], o.Q[i]
		if o.Multq {
			mult := fun.Sign(o.W[i])
			y[i] *= mult
		}
		if o.UseOct {
			x[i] *= tsr.SQ3
			y[i] *= tsr.SQ2by3
		}
		if i == 0 {
			xmi, xma = x[i], x[i]
			ymi, yma = y[i], y[i]
		} else {
			xmi = min(xmi, x[i])
			xma = max(xma, x[i])
			ymi = min(ymi, y[i])
			yma = max(yma, y[i])
		}
		if o.SMPon {
			x[i], y[i], _ = tsr.M_pq_smp(res[i].Sig, o.SMPa, o.SMPb, o.SMPβ, o.SMPϵ)
		}
	}
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	// yield surface
	if o.WithYs && o.m != nil {
		mx, my := 1.0, 1.0
		if o.UseOct {
			mx, my = tsr.SQ3, tsr.SQ2by3
		}
		if o.UsePmin {
			xmi = min(xmi, o.Pmin*mx)
		}
		if o.UsePmax {
			xma = max(xma, o.Pmax*mx)
			yma = max(yma, o.Pmax*my)
		}
		xmi, xma, ymi, yma = o.fix_range(xmi, xmi, xma, ymi, yma)
		if o.PqLims != nil {
			xmi, xma, ymi, yma = o.PqLims[0], o.PqLims[1], o.PqLims[2], o.PqLims[3]
		}
		//io.Pforan("xmi,xma ymi,yma = %v,%v %v,%v\n", xmi,xma, ymi,yma)
		dx := (xma - xmi) / float64(o.NptsPq-1)
		dy := (yma - ymi) / float64(o.NptsPq-1)
		xx := la.MatAlloc(o.NptsPq, o.NptsPq)
		yy := la.MatAlloc(o.NptsPq, o.NptsPq)
		za := la.MatAlloc(o.NptsPq, o.NptsPq)
		zb := la.MatAlloc(o.NptsPq, o.NptsPq)
		var p, q, σa, σb, σc, λ0, λ1, λ2 float64
		v := NewState(len(res[0].Sig), len(res[0].Alp), false)
		for k := 0; k < nr; k++ {
			copy(v.Alp, res[k].Alp)
			v.Dgam = res[k].Dgam
			for i := 0; i < o.NptsPq; i++ {
				for j := 0; j < o.NptsPq; j++ {
					xx[i][j] = xmi + float64(i)*dx
					yy[i][j] = ymi + float64(j)*dy
					p, q = xx[i][j], yy[i][j]
					if o.UseOct {
						p /= tsr.SQ3
						q /= tsr.SQ2by3
					}
					σa, σb, σc = tsr.PQW2O(p, q, o.W[k])
					λ0, λ1, λ2 = tsr.O2L(σa, σb, σc)
					v.Sig[0], v.Sig[1], v.Sig[2] = λ0, λ1, λ2
					ys := o.m.YieldFuncs(v)
					za[i][j] = ys[0]
					if o.nsurf > 1 {
						zb[i][j] = ys[1]
					}
					if o.SMPon {
						xx[i][j], yy[i][j], _ = tsr.M_pq_smp(v.Sig, o.SMPa, o.SMPb, o.SMPβ, o.SMPϵ)
					}
				}
			}
			plt.ContourSimple(xx, yy, za, io.Sf("colors=['%s'], levels=[0], linestyles=['%s'], linewidths=[%g], clip_on=0", o.YsClr0, o.YsLs0, o.YsLw0)+o.ArgsYs)
			if o.nsurf > 1 {
				plt.ContourSimple(xx, yy, zb, io.Sf("colors=['%s'], levels=[0], linestyles=['%s'], linewidths=[%g], clip_on=0", o.YsClr1, o.YsLs1, o.YsLw1)+o.ArgsYs)
			}
		}
	}
	// predictor-corrector
	if len(o.PreCor) > 1 {
		var p, q, pnew, qnew float64
		for i := 1; i < len(o.PreCor); i++ {
			p = tsr.M_p(o.PreCor[i-1])
			q = tsr.M_q(o.PreCor[i-1])
			pnew = tsr.M_p(o.PreCor[i])
			qnew = tsr.M_q(o.PreCor[i])
			if o.UseOct {
				p *= tsr.SQ3
				pnew *= tsr.SQ3
				q *= tsr.SQ2by3
				qnew *= tsr.SQ2by3
			}
			if o.SMPon {
				p, q, _ = tsr.M_pq_smp(o.PreCor[i-1], o.SMPa, o.SMPb, o.SMPβ, o.SMPϵ)
				pnew, qnew, _ = tsr.M_pq_smp(o.PreCor[i], o.SMPa, o.SMPb, o.SMPβ, o.SMPϵ)
			}
			if math.Abs(pnew-p) > 1e-10 || math.Abs(qnew-q) > 1e-10 {
				plt.Arrow(p, q, pnew, qnew, io.Sf("sc=%d, fc='%s', ec='%s'", o.ArrWid, o.ClrPC, o.ClrPC))
			}
		}
	}
	// settings
	if last {
		plt.Equal()
		xl, yl := "$p_{cam}$", "$q_{cam}$"
		if o.UseOct {
			xl, yl = "$p_{oct}$", "$q_{oct}$"
		}
		if o.SMPon {
			xl, yl = "$p_{smp}$", "$q_{smp}$"
		}
		if o.AxLblX != "" {
			xl = o.AxLblX
		}
		if o.AxLblY != "" {
			yl = o.AxLblY
		}
		plt.Gll(xl, yl, "leg_out=1, leg_ncol=4, leg_hlen=1.5")
		if lims, ok := o.Lims["p,q"]; ok {
			plt.AxisLims(lims)
		}
		if lims, ok := o.Lims["p,q,ys"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_oct(x, y []float64, res []*State, sts [][]float64, last bool) {
	// stress path
	nr := len(res)
	k := nr - 1
	var σa, σb, xmi, xma, ymi, yma float64
	for i := 0; i < nr; i++ {
		σa, σb, _ = tsr.PQW2O(o.P[i], o.Q[i], o.W[i])
		x[i], y[i] = σa, σb
		o.maxR = max(o.maxR, math.Sqrt(σa*σa+σb*σb))
		if i == 0 {
			xmi, xma = x[i], x[i]
			ymi, yma = y[i], y[i]
		} else {
			xmi = min(xmi, x[i])
			xma = max(xma, x[i])
			ymi = min(ymi, y[i])
			yma = max(yma, y[i])
		}
	}
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'%s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	// fix range and max radius
	xmi, xma, ymi, yma = o.fix_range(0, xmi, xma, ymi, yma)
	rr := math.Sqrt((xma-xmi)*(xma-xmi) + (yma-ymi)*(yma-ymi))
	if o.maxR < rr {
		o.maxR = rr
	}
	if o.maxR < 1e-10 {
		o.maxR = 1
	}
	if yma > -xmi {
		xmi = -yma
	}
	if o.OctLims != nil {
		xmi, xma, ymi, yma = o.OctLims[0], o.OctLims[1], o.OctLims[2], o.OctLims[3]
	}
	//xmi, xma, ymi, yma = -20000, 20000, -20000, 20000
	// yield surface
	var σcmax float64
	if o.WithYs && o.m != nil {
		//io.Pforan("xmi,xma ymi,yma = %v,%v %v,%v\n", xmi,xma, ymi,yma)
		dx := (xma - xmi) / float64(o.NptsOct-1)
		dy := (yma - ymi) / float64(o.NptsOct-1)
		xx := la.MatAlloc(o.NptsOct, o.NptsOct)
		yy := la.MatAlloc(o.NptsOct, o.NptsOct)
		zz := la.MatAlloc(o.NptsOct, o.NptsOct)
		var λ0, λ1, λ2, σc float64
		v := NewState(len(res[0].Sig), len(res[0].Alp), false)
		for k := 0; k < nr; k++ {
			copy(v.Alp, res[k].Alp)
			v.Dgam = res[k].Dgam
			σc = tsr.M_p(res[k].Sig) * tsr.SQ3
			//σc = 30000
			σcmax = max(σcmax, σc)
			for i := 0; i < o.NptsOct; i++ {
				for j := 0; j < o.NptsOct; j++ {
					xx[i][j] = xmi + float64(i)*dx
					yy[i][j] = ymi + float64(j)*dy
					λ0, λ1, λ2 = tsr.O2L(xx[i][j], yy[i][j], σc)
					v.Sig[0], v.Sig[1], v.Sig[2] = λ0, λ1, λ2
					ys := o.m.YieldFuncs(v)
					zz[i][j] = ys[0]
				}
			}
			plt.ContourSimple(xx, yy, zz, io.Sf("colors=['%s'], levels=[0], linestyles=['%s'], linewidths=[%g], clip_on=0", o.YsClr0, o.YsLs0, o.YsLw0)+o.ArgsYs)

		}
	}
	// predictor-corrector
	if len(o.PreCor) > 1 {
		var σa, σb, σanew, σbnew float64
		for i := 1; i < len(o.PreCor); i++ {
			σa, σb, _ = tsr.M_oct(o.PreCor[i-1])
			σanew, σbnew, _ = tsr.M_oct(o.PreCor[i])
			if math.Abs(σanew-σa) > 1e-7 || math.Abs(σbnew-σb) > 1e-7 {
				//plt.Plot([]float64{σa,σanew}, []float64{σb,σbnew}, "'k+', ms=3, color='k'")
				plt.Arrow(σa, σb, σanew, σbnew, io.Sf("sc=%d, fc='%s', ec='%s'", o.ArrWid, o.ClrPC, o.ClrPC))
			}
			o.maxR = max(o.maxR, math.Sqrt(σa*σa+σb*σb))
			o.maxR = max(o.maxR, math.Sqrt(σanew*σanew+σbnew*σbnew))
		}
	}
	// rosette and settings
	if last {
		tsr.PlotRefOct(o.Phi, σcmax, true)
		tsr.PlotRosette(o.maxR, false, true, true, 6)
		if o.OctAxOff {
			plt.AxisOff()
		}
		plt.Gll("$\\sigma_a$", "$\\sigma_b$", "")
		if lims, ok := o.Lims["oct"]; ok {
			plt.AxisLims(lims)
		}
		if lims, ok := o.Lims["oct,ys"]; ok {
			plt.AxisLims(lims)
		}
	}
}

func (o *Plotter) Plot_s3_s1(x, y []float64, res []*State, sts [][]float64, last bool) {
	// stress path
	nr := len(res)
	k := nr - 1
	x2 := make([]float64, nr)
	var xmi, xma, ymi, yma float64
	for i := 0; i < nr; i++ {
		σ1, σ2, σ3, err := tsr.M_PrincValsNum(res[i].Sig)
		if err != nil {
			chk.Panic("computation of eigenvalues failed", err)
		}
		x[i], y[i] = -tsr.SQ2*σ3, -σ1
		x2[i] = -tsr.SQ2 * σ2
		if i == 0 {
			xmi, xma = x[i], x[i]
			ymi, yma = y[i], y[i]
		} else {
			xmi = min(min(xmi, x[i]), x2[i])
			xma = max(max(xma, x[i]), x2[i])
			ymi = min(ymi, y[i])
			yma = max(yma, y[i])
		}
	}
	plt.Plot(x, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'$\\sigma_3$ %s'", o.Ls, o.Clr, o.Mrk, o.Lbl))
	plt.Plot(x2, y, io.Sf("'r.', ls='%s', clip_on=0, color='%s', marker='%s', label=r'$\\sigma_2$ %s'", o.LsAlt, o.Clr, o.Mrk, o.Lbl))
	plt.PlotOne(x[0], y[0], io.Sf("'bo', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.SpMrk, o.SpMs))
	plt.PlotOne(x[k], y[k], io.Sf("'bs', clip_on=0, color='%s', marker='%s', ms=%d", o.SpClr, o.EpMrk, o.EpMs))
	// yield surface
	if o.WithYs && o.m != nil {
		if o.UsePmin {
			xmi = min(xmi, o.Pmin*tsr.SQ2)
			ymi = min(ymi, o.Pmin)
		}
		if o.UsePmax {
			xma = max(xma, o.Pmax*tsr.SQ2)
			yma = max(yma, o.Pmax)
		}
		xmi, xma, ymi, yma = o.fix_range(0, xmi, xma, ymi, yma)
		if o.S3s1Lims != nil {
			xmi, xma, ymi, yma = o.S3s1Lims[0], o.S3s1Lims[1], o.S3s1Lims[2], o.S3s1Lims[3]
		}
		//io.Pforan("xmi,xma ymi,yma = %v,%v %v,%v\n", xmi,xma, ymi,yma)
		dx := (xma - xmi) / float64(o.NptsSig-1)
		dy := (yma - ymi) / float64(o.NptsSig-1)
		xx := la.MatAlloc(o.NptsSig, o.NptsSig)
		yy := la.MatAlloc(o.NptsSig, o.NptsSig)
		zz := la.MatAlloc(o.NptsSig, o.NptsSig)
		v := NewState(len(res[0].Sig), len(res[0].Alp), false)
		for k := 0; k < nr; k++ {
			copy(v.Alp, res[k].Alp)
			v.Dgam = res[k].Dgam
			for i := 0; i < o.NptsSig; i++ {
				for j := 0; j < o.NptsSig; j++ {
					xx[i][j] = xmi + float64(i)*dx
					yy[i][j] = ymi + float64(j)*dy
					v.Sig[0], v.Sig[1], v.Sig[2] = -yy[i][j], -xx[i][j]/tsr.SQ2, -xx[i][j]/tsr.SQ2
					ys := o.m.YieldFuncs(v)
					zz[i][j] = ys[0]
				}
			}
			plt.ContourSimple(xx, yy, zz, io.Sf("colors=['%s'], levels=[0], linestyles=['%s'], linewidths=[%g], clip_on=0", o.YsClr0, o.YsLs0, o.YsLw0)+o.ArgsYs)
		}
	}
	// predictor-corrector
	if len(o.PreCor) > 1 {
		var σ3, σ1, σ3new, σ1new float64
		for i := 1; i < len(o.PreCor); i++ {
			σ1, _, σ3, _ = tsr.M_PrincValsNum(o.PreCor[i-1])
			σ1new, _, σ3new, _ = tsr.M_PrincValsNum(o.PreCor[i])
			if math.Abs(σ3new-σ3) > 1e-7 || math.Abs(σ1new-σ1) > 1e-7 {
				plt.Arrow(-σ3*tsr.SQ2, -σ1, -σ3new*tsr.SQ2, -σ1new, io.Sf("sc=%d, fc='%s', ec='%s'", o.ArrWid, o.ClrPC, o.ClrPC))
			}
		}
	}
	// settings
	if last {
		plt.Equal()
		plt.Gll("$-\\sqrt{2}\\sigma_3$", "$-\\sigma_1$", "leg=1, leg_out=1, leg_ncol=4, leg_hlen=2")
		if lims, ok := o.Lims["s3,s1"]; ok {
			plt.AxisLims(lims)
		}
		if lims, ok := o.Lims["s3,s1,ys"]; ok {
			plt.AxisLims(lims)
		}
	}
}

// PlotRamp plots the ramp function (contour)
func (o *Plotter) DrawRamp(xmi, xma, ymi, yma float64) {
	if o.Rmpf == nil {
		o.set_empty()
		return
	}
	if o.NptsRmp < 2 {
		o.NptsRmp = 101
	}
	if math.Abs(xma-xmi) < 1e-5 {
		xmi, xma = -0.1, 0.1
	}
	if math.Abs(yma-ymi) < 1e-5 {
		ymi, yma = -0.1, 0.1
	}
	xx := la.MatAlloc(o.NptsRmp, o.NptsRmp)
	yy := la.MatAlloc(o.NptsRmp, o.NptsRmp)
	zz := la.MatAlloc(o.NptsRmp, o.NptsRmp)
	dx := (xma - xmi) / float64(o.NptsRmp-1)
	dy := (yma - ymi) / float64(o.NptsRmp-1)
	for i := 0; i < o.NptsRmp; i++ {
		for j := 0; j < o.NptsRmp; j++ {
			xx[i][j] = xmi + float64(i)*dx
			yy[i][j] = ymi + float64(j)*dy
			zz[i][j] = xx[i][j] - o.Rmpf(xx[i][j]+yy[i][j])
		}
	}
	plt.ContourSimple(xx, yy, zz, "colors=['blue'], linewidths=[2], levels=[0]")
}

// Save saves figure
func (o *Plotter) Save(typ, num string) {
	if o.PlotFcn != nil {
		o.PlotFcn()
	}
	ext := ".png"
	if o.UseEps {
		ext = ".eps"
	}
	if o.SaveFnk != "" {
		if o.SaveDir != "" {
			plt.SaveD(o.SaveDir, o.SaveFnk+typ+num+ext)
			return
		}
		plt.Save(o.SaveFnk + typ + num + ext)
	}
}

// subplot sets subplot
func (o *Plotter) Subplot() {
	if o.Split {
		plt.Clf()
		return
	}
	plt.Subplot(o.Nrow, o.Ncol, o.Pidx)
	o.Pidx += 1
}

// Auxiliary ////////////////////////////////////////////////////////////////////////////////////

// set_default_clr_mrk sets default colors and markers
func (o *Plotter) set_default_clr_mrk() {
	if o.ClrPC == "" {
		o.ClrPC = "#85b9ff"
	}
	if o.ArrWid == 0 {
		o.ArrWid = 10
	}
	if o.Clr == "" {
		o.Clr = "red"
	}
	if o.Mrk == "" {
		o.Mrk = ""
	}
	if o.Ls == "" {
		o.Ls = "-"
	}
	if o.LsAlt == "" {
		o.LsAlt = "--"
	}
	if o.SpMrk == "" {
		o.SpMrk = "o"
	}
	if o.EpMrk == "" {
		o.EpMrk = "s"
	}
	if o.SpClr == "" {
		o.SpClr = "black"
	}
	if o.EpClr == "" {
		o.EpClr = "black"
	}
	if o.SpMs == 0 {
		o.SpMs = 3
	}
	if o.EpMs == 0 {
		o.EpMs = 3
	}
	if o.YsClr0 == "" {
		o.YsClr0 = "green"
	}
	if o.YsClr1 == "" {
		o.YsClr1 = "cyan"
	}
	if o.YsLs0 == "" {
		o.YsLs0 = "-"
	}
	if o.YsLs1 == "" {
		o.YsLs1 = "--"
	}
	if o.YsLw0 < 0.1 {
		o.YsLw0 = 0.7
	}
	if o.YsLw1 < 0.1 {
		o.YsLw1 = 0.7
	}
}

// fix_range fixes range of scale to plot contours
func (o *Plotter) fix_range(middle, Xmi, Xma, Ymi, Yma float64) (xmi, xma, ymi, yma float64) {
	xmi, xma, ymi, yma = Xmi, Xma, Ymi, Yma
	if xma-xmi < 1e-7 {
		xmi, xma = -1+middle, 1+middle
	}
	if yma-ymi < 1e-7 {
		ymi, yma = -1+middle, 1+middle
	}
	m := o.YsRangeM
	if m < 1e-7 {
		m = 0.2
	}
	xmi -= m * (xma - xmi)
	xma += m * (xma - xmi)
	ymi -= m * (yma - ymi)
	yma += m * (yma - ymi)
	return
}

func (o Plotter) calc_x(p float64) float64 {
	return math.Log(1.0 + (p+o.Pt)/o.Pr)
}

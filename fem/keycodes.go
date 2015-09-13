// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gosl/io"
)

/*
func GetIntegrationPoints(nip, nipf int, cell *inp.Cell) (ipsElem, ipsFace []shp.Ipoint, err error) {

	// get integration points of element
	ipsElem, err = shp.GetIps(cell.Type, nip)
	if err != nil {
		err = chk.Err("cannot get integration points for element with shape type=%q and nip=%d\n%v", cellType, nip, err)
		return
	}

	// get integration points of face
	faceType := shp.GetFaceType(cellType)
	ipsFace, err = shp.GetIps(faceType, nipf)
	if err != nil {
		err = chk.Err("cannot get integration points for face with face-shape type=%q and nip=%d\n%v", faceType, nip, err)
	}
	return
}
*/

func GetSolidFlags(axisym, pstress bool, extra string) (useB, debug bool, thickness float64) {

	// defaults
	useB = false
	debug = false
	thickness = 1.0

	// flag: use B matrix
	if s_useB, found := io.Keycode(extra, "useB"); found {
		useB = io.Atob(s_useB)
	}

	// fix useB flag in case of axisymmetric simulation
	if axisym {
		useB = true
	}

	// flag: thickess => plane-stress
	if s_thick, found := io.Keycode(extra, "thick"); found {
		thickness = io.Atof(s_thick)
	}

	// fix thickness flag
	if !pstress {
		thickness = 1.0
	}

	// flag: debug
	if s_debug, found := io.Keycode(extra, "debug"); found {
		debug = io.Atob(s_debug)
	}
	return
}

func GetSeepFaceFlags(extra string) (Macaulay bool, BetRamp, Kappa float64) {

	// defaults
	Macaulay = false
	BetRamp = math.Ln2 / 0.01
	Kappa = 1.0

	// use macaulay function ?
	if s_mac, found := io.Keycode(extra, "mac"); found {
		Macaulay = io.Atob(s_mac)
	}

	// coefficient for smooth ramp function
	if s_bet, found := io.Keycode(extra, "bet"); found {
		BetRamp = io.Atof(s_bet)
	}

	// κ coefficient
	if s_kap, found := io.Keycode(extra, "kap"); found {
		Kappa = io.Atof(s_kap)
	}
	return
}

func GetContactFaceFlags(extra string) (Macaulay bool, BetRamp, Kappa float64) {

	// defaults
	Macaulay = false
	BetRamp = math.Ln2 / 0.01
	Kappa = 1.0

	// use macaulay function ?
	if s_mac, found := io.Keycode(extra, "mac"); found {
		Macaulay = io.Atob(s_mac)
	}

	// coefficient for smooth ramp function
	if s_bet, found := io.Keycode(extra, "bet"); found {
		BetRamp = io.Atof(s_bet)
	}

	// κ coefficient
	if s_kap, found := io.Keycode(extra, "kap"); found {
		Kappa = io.Atof(s_kap)
	}
	return
}

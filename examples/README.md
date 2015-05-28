# Gofem examples

Check gofem/inp/sim.go
Check where 2/3, 5/6, 8/9 came from? (Zienkiewicz)

## Pullout test: rjoing_ex06_pullout

3D beam with one steel rod being pulled out

The plot shows the displacement of rod nodes for a number of time outputs. 'y' is the length of the
rod starting from 0.20 because the rod is inside the beam.

## Horizontal permeameter: seep_simple_flux

Cylinder with porous media under seepage only. With gravity. To the left hand side hydrostatic
pressure is kept constant. To the right hand side, pressure is decrease by the following _shift_

_hst_ specifies a **shift** of liquid pressure according to:

$$
  pl(t,z) = pl_0(z) - shift(t)
$$

_hydrost_ must be true in _stages_

_dvgctrl_ divergence control will reduce the time step in case divergence occurs

The plot shows ...


## Plotting

Use package _out_

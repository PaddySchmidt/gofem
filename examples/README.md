# Gofem &ndash; examples

## Summary
1.  **dynamics_sgbook** -- Solid dynamics. Examples from Smith, Griffiths and Margetts book [1].
2.  **rjoint_ex01_curved** -- Rod-Joint model for reinforced concrete [2,3].
3.  **rjoint_ex06_pullout** -- Pullout test in reinforced concrete [2,3].
4.  **seep_ex01_freesurf** -- Transient simulation with seepage face. Rectangle domain. [4]
5.  **seep_ex02_freesurf** -- Transient simulation with seepage face. Earthen Slope. [4]
6.  **seep_simple_flux** -- Simple flux in a horizontal permeameter. 2D and 3D simulations with computation of water discharge.
7.  **spo751_pressurised_cylinder** -- Internally pressurised cylinder. Plasticity problem number 7.5.1 (page 244) from Souza-Peric-Owen (SPO) book [5].
8.  **spo754_strip_footing_collapse** -- Strip-footing collapse. Plasticity problem number 7.5.4 (page 252) from Souza-Peric-Owen (SPO) book [5].
9.  **up_3mcolumn_desiccation** -- Desiccation of porous column. Unsaturated porous media simulation from [6].
10. **up_indentation2d_unsat** -- Indentation of unsaturated porous medium. Similar to example in [6].

## References
1. Smith, Griffiths and Margetts (2013) Programming the Finite Element Method. 5th Edition. Wiley. ISBN: 978-1-119-97334-8. 682 pages.
2. Durand, Farias and Pedroso (2015) A rod-joint finite element for reinforced solids without the need for mesh compatibility.
3. Durand, Farias and Pedroso (2015) Computing intersections between non-compatible curves and finite elements. Computational Mechanics. [doi:10.1007/s00466-015-1181-y](http://dx.doi.org/10.1007/s00466-015-1181-y).
4. Pedroso (2015) A solution to transient seepage in unsaturated porous media. Computer Methods in Applied Mechanics and Engineering. 285:791-816. [doi:10.1016/j.cma.2014.12.009](http://dx.doi.org/10.1016/j.cma.2014.12.009).
5. de Souza Neto, E. A., Peric, D., and Owen, D. R. J. (2008) Computational Methods for Plasticity: Theory and Applications. Wiley. 791 pages. [doi:10.1002/9780470694626](http://dx.doi.org/10.1002/9780470694626).
6. Pedroso (2015) A consistent u-p formulation for porous media with hysteresis. International Journal for Numerical Methods in Engineering. 101:606-634. [doi:10.1002/nme.4808](http://dx.doi.org/10.1002/nme.4808).



# 1 dynamics_sgbook -- Solid dynamics

## 1.1 Forced vibrations. Single beam element with dynamic load

See page 487 of [1].

### 1.1.1 Geometry, mesh and boundary conditions

<div id="container">
<p><img src="dynamics_sgbook/figs/sg111msh.png" width="400"></p>
Finite element mesh.
</div>

### 1.1.2 Results

<div id="container">
<p><img src="dynamics_sgbook/figs/sg111.png" width="400"></p>
Deflection at the right hand side.
</div>



## 1.2 Forced vibrations. Solid element with dynamic load

See page 491 of [1].

### 1.2.1 Geometry, mesh and boundary conditions

<div id="container">
<p><img src="dynamics_sgbook/figs/sg114msh.png" width="400"></p>
Finite element mesh.
</div>

### 1.2.2 Results

<div id="container">
<p><img src="dynamics_sgbook/figs/sg114.png" width="400"></p>
Deflection.
</div>



## 1.3 Plastic slab with impacting distributed load

See page 515 of [1].

### 1.3.1 Geometry, mesh and boundary conditions

<div id="container">
<p><img src="dynamics_sgbook/figs/sg1121msh.png" width="400"></p>
Finite element mesh.
</div>

### 1.3.2 Results

<div id="container">
<p><img src="dynamics_sgbook/figs/sg1121.png" width="400"></p>
Deflection.
</div>



# 2 rjoint_ex01_curved -- Rod-Joint model for reinforced concrete

A curved rod is immersed in a parallelepiped and has a force applied to one of its end.

## 2.1 Mesh

<div id="container">
<p><img src="rjoint_ex01_curved/figs/mesh.png" width="400"></p>
Mesh with curved rod (3-node element shown by straight lines in red).
</div>

## 2.2 Results

<div id="container">
<p><img src="rjoint_ex01_curved/figs/rjoint01.png" width="400"></p>
Shear stress along joint element.
</div>



# 3 rjoint_ex06_pullout -- Pullout test in reinforced concrete

Pullout test in reinforced concrete block.

## 3.1 Mesh

<div id="container">
<p><img src="rjoint_ex06_pullout/figs/mesh.png" width="400"></p>
Finite element mesh.
</div>

## 3.2 Results

The plot shows the displacement of rod nodes for a number of time outputs. 'y' is the length of the rod starting from 0.20 because the rod is inside the beam.

<div id="container">
<p><img src="rjoint_ex06_pullout/figs/rjoint_ex06_pullout_o2elast.png" width="400"></p>
Displacements of nodes along rod (rebar).
</div>



# 4 seep_ex01_freesurf -- Transient simulation with seepage face

Modelling the seepage face using enriched element as proposed in [4]. The right hand side of a block of porous medium has the liquid pressure decreased thus inducing flow from left to right. [See video of similar simulation here](https://vimeo.com/113379824).

## 4.1 Mesh

<div id="container">
<p><img src="seep_ex01_freesurf/figs/mesh.png" width="400"></p>
Finite element mesh.
</div>

## 4.2 Results

<div id="container">
<p><img src="seep_ex01_freesurf/figs/coarse.png" width="400"></p>
Results: pressure and saturation distribution at the end of the simulation.
</div>


# 5 seep_ex02_freesurf -- Transient simulation with seepage face

Modelling the seepage face using enriched element as proposed in [4]. The right hand side of a earthen slope has the liquid pressure decreased thus inducing flow from left to right. [See video of similar simulation here](https://vimeo.com/113380931).

## 5.1 Mesh

<div id="container">
<p><img src="seep_ex02_freesurf/figs/mesh.png" width="400"></p>
Finite element mesh.
</div>

## 5.2 Results

<div id="container">
<p><img src="seep_ex02_freesurf/figs/struct.png" width="400"></p>
Results: pressure and saturation distribution at the end of the simulation.
</div>



# 6 seep_simple_flux -- Simple flux in a horizontal permeameter

Cylinder with porous media under seepage only. With gravity. To the left hand side hydrostatic
pressure is kept constant. To the right hand side, pressure is decrease by the following _shift_

_hst_ specifies a **shift** of liquid pressure according to:

$$
  pl(t,z) = pl_0(z) - shift(t)
$$

_hydrost_ must be true in _stages_

_dvgctrl_ divergence control will reduce the time step in case divergence occurs

## 6.1 2D: Mesh

<div id="container">
<p><img src="seep_simple_flux/figs/d2-coarse-msh.png" width="400"></p>
Finite element mesh. Rectangle with 10 x 3 units of length.
</div>

## 6.2 2D: Results

<div id="container">
<p><img src="seep_simple_flux/figs/seep_simple_flux_d2-simple-flux.png" width="400"></p>
Results.
 (a) pressure along vertical section at middle of rectangle (x=5.0);
 (b) pressure along horizontal section at lower part of the rectangle (y=0.0);
 (c) filter velocity at x=5 and y=3.
</div>

<div id="container">
<p><img src="seep_simple_flux/figs/d2-results.png" width="400"></p>
Seepage velocity.
</div>

## 6.3 3D: Results

<div id="container">
<p><img src="seep_simple_flux/figs/d3-results.png" width="400"></p>
Seepage velocity.
</div>



# 7 spo751_pressurised_cylinder -- Internally pressurised cylinder

See page 244 of [5]. Pressure is applied to the edge tagged with -10. The material is modelled using the von Mises model.

## 7.1 Mesh

<div id="container">
<p><img src="spo751_pressurised_cylinder/figs/mesh.png" width="400"></p>
Finite element mesh.
</div>

## 7.2 Results

<div id="container">
<p><img src="spo751_pressurised_cylinder/figs/gofem_spo751_disp.png" width="400"></p>
Displacement of the outer face; node # 20 with tag -202.
</div>

<div id="container">
<p><img src="spo751_pressurised_cylinder/figs/gofem_spo751_sr.png" width="400"></p>
Radial stress component of integration points along a line from the inner to the outer faces.
</div>

<div id="container">
<p><img src="spo751_pressurised_cylinder/figs/gofem_spo751_st.png" width="400"></p>
Tangential stress component of integration points along a line from the inner to the outer faces.
</div>



# 8 spo754_strip_footing_collapse -- Strip-footing collapse

See page 252 of [5]. Displacement is applied to the region corresponding to the footing at the left hand side. von Mises model is used.

## 8.1 Mesh

<div id="container">
<p><img src="spo754_strip_footing_collapse/figs/mesh.png" width="400"></p>
Finite element mesh.
</div>

## 8.2 Results

TODO.



# 9 up_3mcolumn_desiccation -- Desiccation of porous column

Desiccation of a 3 metres high column of porous medium. The bottom of the column has pressure specified. This pressure varies causing desiccation of column.

## 9.1 Mesh

<div id="container">
<p><img src="up_3mcolumn_desiccation/figs/mesh.png"></p>
Finite element mesh.
</div>

### 9.2 Decreasing and increasing the pressure.

<div id="container">
<p><img src="up_3mcolumn_desiccation/figs/up_3mcolumn_dessication_onepulse-qua9co.png"></p>
Results.
</div>

### 9.3 Linearly decreasing the pressure to reach unsaturated state.

<div id="container">
<p><img src="up_3mcolumn_desiccation/figs/up_3mcolumn_dessication_linear-qua9co.png"></p>
Results.
</div>

### 9.4 Wetting of initially unsaturated column.

<div id="container">
<p><img src="up_3mcolumn_desiccation/figs/up_3mcolumn_dessication_wet-linear-qua9co.png"></p>
Results.
</div>



# 10 up_indentation2d_unsat -- Indentation of unsaturated porous medium

This example shows the indentation of an unsaturated soil. First, the domain is made unsaturated by means of decreasing the pressure at the bottom. The next simulation the applies the indenter at the surface. Two ways of modelling this problem are considered:

* a and b simulations -- lower the pressure at bottom and apply the load
* c and d simulations -- dry the medium by means of flux prescribed at the surface then apply the load

## 10.1 Mesh

The size of the domain is 3 x 3 in units of length.

<div id="container">
<p><img src="up_indentation2d_unsat/figs/mesh.png"></p>
Finite element mesh.
</div>

## 10.2 Results -- first stage by lowering pressure

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_a-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. First stage: lowering pressure.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_b-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. Second stage: applying indenter.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_lrm_a-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. First stage: liquid retention behaviour under the indenter.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_lrm_b-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. Second stage: liquid retention behaviour under the indenter.
</div>

## 10.3 Results -- first stage by drying surface

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_c-coarse-elast-d2-q9.png"></p>
Results from first stage by drying surface. First stage: drying surface by applying flux.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_d-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. Second stage: applying indenter.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_lrm_c-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. First stage: liquid retention behaviour under the indenter.
</div>

<div id="container">
<p><img src="up_indentation2d_unsat/figs/up_indentation2d_unsat_lrm_d-coarse-elast-d2-q9.png"></p>
Results from first stage by lowering pressure. Second stage: liquid retention behaviour under the indenter.
</div>

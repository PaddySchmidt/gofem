from FEMmesh import FEMmesh
from FEMsolver import FEMsolver

# 1) input mesh
# mesh
verts = [[0, -101,  0.0, 0.0],
         [1,    0,  2.0, 0.0],
         [2, -102,  4.0, 0.0],
         [3, -103,  6.0, 0.0],
         [4,    0,  8.0, 0.0],
         [5, -104, 10.0, 0.0],
         [6,    0,  2.0, 4.0],
         [7,    0,  4.0, 4.0],
         [8,    0,  6.0, 4.0],
         [9,    0,  8.0, 4.0]]
cells = [# lower horizontal
         [ 0, -1, [0,1]],
         [ 1, -1, [1,2]],
         [ 2, -1, [2,3]],
         [ 3, -1, [3,4]],
         [ 4, -1, [4,5]],
         # diagonals
         [ 5, -1, [0,6]],
         [ 6, -1, [2,6]],
         [ 7, -1, [3,7]],
         [ 8, -1, [2,8]],
         [ 9, -1, [3,9]],
         [10, -1, [5,9]],
         # verticals
         [11, -1, [1,6]],
         [12, -1, [2,7]],
         [13, -1, [3,8]],
         [14, -1, [4,9]],
         # upper horizontal
         [15, -1, [6,7]],
         [16, -1, [7,8]],
         [17, -1, [8,9]]]
m = FEMmesh(verts, cells)
#m.draw()
#m.show()

# 2) create dictionary with all parameters
p = {-1:{'E':2.1e8,'A':0.003}}

# 3) allocate fem solver object
s = FEMsolver(m,'EelasticRod',p)

# 4) set boundary conditions
vb = {-101:{'ux':0.0,'uy':0.0},
      -102:{'fy':-320.0},
      -103:{'fy':-350.0},
      -104:{'uy':0.0}} # fy:[MN]
s.set_bcs(vb=vb) # vertex bry conds

# 5) solve equilibrium problem
s.solve_steady(True)
s.calc_secondary()

# results
all_ux = s.get_u()['ux']
all_uy = s.get_u()['uy']

# to compare with gofem
l = '[\n{\n  "Kmats":[\n'
for i, e in enumerate(s.elems):
    if i > 0: l += ',\n'
    l += '    [\n'
    for j in range(4):
        if j > 0: l += ',\n'
        l += '      ['
        for k in range(4):
            if k > 0: l += ','
            l += '%23.15e' % e.K[j,k]
        l += ']'
    l += '\n    ]'
l += '\n  ],\n  "disp":[\n'
for i, ux in enumerate(all_ux):
    if i > 0: l += ',\n'
    l += '    [%23.15e,%23.15e]' % (ux, all_uy[i])
l += '\n  ],\n  "sigmas":[\n'
for ic in range(s.m.nc):
    if ic > 0: l += ',\n'
    l += '    [[%23.15e]]' % s.esvs['sa'][ic][0]
l += '\n  ]\n}\n]'

print l

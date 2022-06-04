import gmsh
import math
import sys
import numpy as np
import mesh_con

gmsh.initialize()
gmsh.model.add("testsimp")
p = 1
k = 0.5
c1 = 0.5
c2 = -0.5
r = 0.2
lc = 2e-2
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(p, p, 0, lc, 2)
gmsh.model.geo.addPoint(-p, p, 0, lc, 3)
gmsh.model.geo.addPoint(-p, -p, 0, lc, 4)
gmsh.model.geo.addPoint(p, -p, 0, lc, 5)
gmsh.model.geo.addLine(2, 3, 1)
gmsh.model.geo.addLine(3, 4, 2)
gmsh.model.geo.addLine(4, 5, 3)
gmsh.model.geo.addLine(5, 2, 4)
gmsh.model.geo.addCurveLoop([1,2,3,4], 1)
gmsh.model.geo.addPoint(c1+r, c1, 0, lc, 101)
gmsh.model.geo.addPoint(c1, c1+r, 0, lc, 102)
gmsh.model.geo.addPoint(c1-r, c1, 0, lc, 103)
gmsh.model.geo.addPoint(c1, c1-r, 0, lc, 104)
gmsh.model.geo.addPoint(c1, c1, 0, lc, 100)
gmsh.model.geo.addCircleArc(101, 100, 102, 5)
gmsh.model.geo.addCircleArc(102, 100, 103, 6)
gmsh.model.geo.addCircleArc(103, 100, 104, 7)
gmsh.model.geo.addCircleArc(104, 100, 101, 8)
gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 2)
gmsh.model.geo.addPoint(c2, c2, 0, lc, 200)
gmsh.model.geo.addPoint(c2+r, c2, 0, lc, 201)
gmsh.model.geo.addPoint(c2, c2+r, 0, lc, 202)
gmsh.model.geo.addPoint(c2-r, c2, 0, lc, 203)
gmsh.model.geo.addPoint(c2, c2-r, 0, lc, 204)
gmsh.model.geo.addCircleArc(201, 200, 202, 9)
gmsh.model.geo.addCircleArc(202, 200, 203, 10)
gmsh.model.geo.addCircleArc(203, 200, 204, 11)
gmsh.model.geo.addCircleArc(204, 200, 201, 12)
gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 3)
gmsh.model.geo.addPlaneSurface([1,2,3], 1)
gmsh.model.geo.addPlaneSurface([2], 2)
gmsh.model.geo.addPlaneSurface([3], 3)
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(1, [1,2,3,4], 4)
gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.addPhysicalGroup(2, [2], 2)
gmsh.model.addPhysicalGroup(2, [3], 3)
gmsh.model.mesh.generate(2)
gmsh.write("Input/testsimp1.mesh")
mesh_con.convert_2d("Input/testsimp1.mesh", "Input/testsimp1.h5")
gmsh.fltk.run()
gmsh.finalize()
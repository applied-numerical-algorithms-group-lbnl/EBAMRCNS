
This is a code that solves the compressible Navier Stokes equations
in the context of complex geometry.    The boundary conditions and 
initial conditions live in the IBC classes.  

Shocktube:

The boundary conditions in this
particular example are no flow/no slip velocity/insulated and the intitial 
condition is  jump in state (think of it as a shock tube with complex geometries). 
There are other implemented in the src directory.

inflow/outflow
For this  one you specify a preshock state and the mach number for your calculation.

See src/GodunovGeom.cpp for the available geometries.

The interface is hard to use but the results are cool.

Chombo documentation:
http://chombo.lbl.gov

Here is the paper that describes the algorithm:
@ARTICLE{compress_ns,
author ={Daniel T. Graves and  Phillip Colella and David Modiano and Jeffrey Johnson  and  Bjorn Sjogreen and Xinfeng Gao},
title = {A Cartesian Grid Embedded Boundary Method for the Compressible {N}avier {S}tokes Equations},
journal = {Communications in Applied Mathematics and Compuational Science},
year = 2013,
volume = 8,
number = 1,
pages = "99-122"
}

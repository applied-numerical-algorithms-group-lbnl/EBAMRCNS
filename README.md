#EBAMRCNS

## Summary

This is a code that solves the compressible Navier Stokes equations
in the context of complex geometry.    The boundary conditions and 
initial conditions live in the IBC classes.   This is all written in the
Chombo infrastructure from LBNL.  See [the Chombo web page:](http://chombo.lbl.gov) for  Chombo documentation. 

## usage :

mpirun -np NPROC navier2d...ex blah.inputs

where NPROC is the number of processors

blah.inputs is the input file name 

###Shocktube  initial/boundary condition


The initial conditions here are zero velocity and  a sharp high/low pressure/density
interface. The boundary conditions in this
particular example are no flow/no slip velocity/insulated and the intitial 
condition is  jump in state (think of it as a shock tube with complex geometries). 
There are other implemented in the src directory.

###input params you should be ok with changing:

1.   restart file        --- starting from checkpoint thing
2.   logflag             --- whether to take logs of pressure and density for output
3.   cfl                 --- cfl number (needs to be less than 1, preferably 0.1-0.5
4.   max_step            --- maximum number of time steps
5.   tag_buffer_size     --- number of cells added around each tagged cell (should be at least 1)
6.   regrid_interval     --- how often to regrid
7.   gamma               --- ratio of specfic heats
8.   domain_length       --- length of the domain (only the x value used, I think)
9.    max_level           --- highest amr level number
10. n_cell              --- number of cells on level 0
11. ref_ratio           --- refinement ratios
12. max_grid_size       --- maximum size of any box in domain
13. checkpoint_interval --- how often to checkpoint
14. plot_interval       --- how often to write plot files
15.  which_geom          --- which geometric configurtion you are using (1 ramp, 2 slab, 4 cylinder, 5 sphere... see GodunovGeom.cpp for more)
16. ramp_normal         --- normal vector of the ramp (for which_geom == 1)
17. ramp_alpha          --- y intercept of the ramp
18. do_diffusion        --- whether to turn on or off diffusion terms ( false == inviscid Euler)
19. mu_viscosity        --- viscosity 
20. lambda_viscosity    --- set this = -2/3 mu_viscosity
21. explosion_p0        --- low pressure
22. explosion_r0        --- low density
23. explosion_p1        --- high pressure
24. explosion_r1        --- high density
25. explosion_size      --- where the interface lies


=-=-=-=-=-=-=-=-
###inflow/outflow initial/boundary condition
=-=-=-=-=-=-=-=-
This is an inflow/outflow problem.
For this  one you specify a preshock state and the mach number for your calculation
as well as where the shock starrts

##input params you should be ok with changing:

1. restart file        --- starting from checkpoint thing
2. logflag             --- whether to take logs of pressure and density for output
3. cfl                 --- cfl number (needs to be less than 1, preferably 0.1-0.5
4. max_step            --- maximum number of time steps
5. tag_buffer_size     --- number of cells added around each tagged cell (should be at least 1)
6. regrid_interval     --- how often to regrid
7. gamma               --- ratio of specfic heats
8. domain_length       --- length of the domain (only the x value used, I think)
9. max_level           --- highest amr level number
10. n_cell              --- number of cells on level 0
11. ref_ratio           --- refinement ratios
12. max_grid_size       --- maximum size of any box in domain
13. checkpoint_interval --- how often to checkpoint
14. plot_interval       --- how often to write plot files
15. which_geom          --- which geometric configurtion you are using (1 ramp, 2 slab, 4 cylinder, 5 sphere... see GodunovGeom.cpp for more)
16. ramp_normal         --- normal vector of the ramp (for which_geom == 1)
17. ramp_alpha          --- y intercept of the ramp
18. do_diffusion        --- whether to turn on or off diffusion terms ( false == inviscid Euler)
19. mu_viscosity        --- viscosity 
20. lambda_viscosity    --- set this = -2/3 mu_viscosity
21. preshockdense       --- pre shock density
22. preshockpress       --- pre shock pressure
23. shock_mach          --- mach number of the shcok
24. shock_center        --- where in the domain the shock lives 

##See src/GodunovGeom.cpp for the available geometries.

##The interface is hard to use but the results are cool.

##[Chombo documentation](http://chombo.lbl.gov)

##Here is the paper that describes the algorithm:
@ARTICLE{compress_ns,
author ={Daniel T. Graves and  Phillip Colella and David Modiano and Jeffrey Johnson  and  Bjorn Sjogreen and Xinfeng Gao},
title = {A Cartesian Grid Embedded Boundary Method for the Compressible {N}avier {S}tokes Equations},
journal = {Communications in Applied Mathematics and Compuational Science},
year = 2013,
volume = 8,
number = 1,
pages = "99-122"
}

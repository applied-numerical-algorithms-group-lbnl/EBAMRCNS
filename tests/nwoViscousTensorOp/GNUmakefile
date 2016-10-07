
#these you can change
DIM=3
#turns on or off symbol table
DEBUG=FALSE 
#optimization can be TRUE, FALSE or HIGH
OPT=HIGH
#parallel can be TRUE or FALSE
MPI=TRUE
#threading can be TRUE or FALSE
OPENMPCC =FALSE 


#leave this stuff alone
USE_EB=TRUE

ebase = nwoViscousDriver


CHOMBO_HOME = ../../chombo_lib

LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTools BoxTools


base_dir = .
src_dirs = ../src


INPUT = viscous.inputs

include $(CHOMBO_HOME)/mk/Make.example



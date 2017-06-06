#leave these alone
USE_EB = TRUE
ebase = navierDriver

#these you can change
DIM=2
#turns on or off symbol table
DEBUG=FALSE 
#optimization can be TRUE, FALSE or HIGH
OPT=HIGH
#parallel can be TRUE or FALSE
MPI=TRUE
#threading can be TRUE or FALSE
OPENMPCC =FALSE

# the location of Chombo lib dir
CHOMBO_HOME = ../chombo_lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
INPUT = ramp.shortrun.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets


## Note: requires USE_EB=TRUE
## Note: requires USE_NEW_EB=TRUE

USE_EB := TRUE
USE_PETSC=TRUE
# the base name(s) of the application(s) in this directory
ebase = marcExample1

# the location of Chombo lib dir
CHOMBO_HOME = ../../chombo_lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTools EBTools Workshop AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = 

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

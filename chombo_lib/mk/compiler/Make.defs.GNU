### This is a compiler-specific definitions file for GNU compilers.

## It sets values for the defaults variables ($def*), which will be used
## if the calling makefile doesn't set the corresponding variables.
## The corresponding variables are the ones without the "def" prefix.

## Defaults variables:
##  defcppcomflags = C-preprocessor options for both C++ and Fortran code
##  defcppdbgflags = C-preprocessor options for both C++ and Fortran code when DEBUG=TRUE
##  defcppoptflags = C-preprocessor options for both C++ and Fortran code when OPT=TRUE or HIGH
##  defcxxcppflags = C-preprocessor options for only C++ code
##  defcxxcomflags = C++ compiler/linker options
##  defcxxdbgflags : for C++ compiling and linking when DEBUG=TRUE
##  defcxxoptflags : for C++ compiling and linking when OPT=TRUE or HIGH
##  defcxxprofflags: for C++ compiling and linking when PROFILE=TRUE
##  deffcppflags   = C-preprocessor options for only Fortran code
##  deffcomflags   = Fortran compiler options
##  deffdbgflag    : for Fortran compiling when DEBUG=TRUE
##  deffoptflag    : for Fortran compiling when OPT=TRUE or HIGH
##  deffprofflags  : for Fortran compiling when PROFILE=TRUE
##  defflibflags   : library options needed to link the Fortran code
##  defldcomflags  = linker options
##  deflddbgflags  = linker options when DEBUG=TRUE
##  defldoptflags  = linker options when OPT=TRUE or HIGH
##  defldprofflags = linker options when PROFILE=TRUE
##  cppcallsfort   = preprocessor macro as used in FORT_PROTO.H to specify
##                    how the C++ code should construct names for calling Fortran
##
## Special object models:
##  USE_64 = TRUE   Compile/link with 64bit pointers and longs
##
##  Note: 64bit only works on AMD64 or PowerPC computers with gcc 3.3 or later.  <dbs,Mar04>
##
## Notes:
##  GNU C/C++ compilers set the preprocessor macros: (use -E -dM options to display)
##    __GNUC__       (int) defined for gcc or g++, == major version number
##    __GNUG__       (int) defined for g++, == major version number
##    __GNUC_MINOR__ (int) defined for gcc or g++, == minor version number

makefiles+=compiler/Make.defs.GNU

cxxname := $(notdir $(firstword $(subst -, ,$(CXX))))
fname   := $(notdir $(firstword $(subst -, ,$(FC))))


################################################################

ifeq ($(cxxname),g++)
  _gppversion := $(subst ., ,$(shell $(CXX) -dumpversion))
  _gppmajorver:= $(word 1,$(_gppversion))
  _gppminorver:= $(word 2,$(_gppversion))

  # `-ftemplate-depth' is needed to keep mangled function names from
  #  getting too long for some linkers to handle.
  # HDF5 uses `long long', so disable warnings about it.
  # `-Wno-sign-compare' turns off warnings about comparing signed and unsigned int -
  #   this is a meaningful warning, but it appears way too often in our code.
  # -Wno-unused-local-typedefs turns off a set of annoying "typedef locally 
  #defined but not used warnings in RefCountedPtr.H

  _cxxbaseflags := -std=c++11 -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-55  -Wno-unused-local-typedefs

  # g++ v3.4.0 started objecting to some uses of HOFFSET()
  # TJL - Changed moving from csh to sh...
  ifeq (00,$(shell test $(_gppmajorver) -ge 3 ; echo $$?)$(shell test $(_gppminorver) -ge 4 ; echo $$?))
    _cxxbaseflags += -Wno-invalid-offsetof
  endif

  ifeq (0,$(shell test $(_gppmajorver) -gt 3 ; echo $$?))
    _cxxbaseflags += -Wno-variadic-macros
  endif

  # NOTE: -malign-double causes MPI programs that use pout() to
  # segfault in ~ostream(), so dont use it

ifeq ($(OPENMPCC),TRUE)
  _cxxbaseflags += -fopenmp
else
  _cxxbaseflags += -Wno-unknown-pragmas
endif

  defcxxoptflags := -march=native -O3 $(_cxxbaseflags)
  defcxxdbgflags := -march=native -g -pedantic -Wall $(_cxxbaseflags)


  ##  # together so force the preprocessor to use the Fortran compiler
  ifeq ($(system),Darwin)
    CH_CPP = /usr/bin/cpp -E -P
  endif

  # the preprocessor in v3 or later strips "//" by default,
  # but this is a valid Fortran operator, so force it to leave comments
  ifeq (0,$(shell test $(_gppmajorver) -ge 3 ; echo $$?))
    CH_CPP += 
  endif

  # the preprocessor for v4.8+ adds some header gibberish that we 
  # suppress with -nostdinc.
  ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 8 ; echo $$?))
    CH_CPP += -nostdinc
  endif

endif


################################################################

# GNU Fortran95 compiler in GCC version 4
ifeq ($(fname),gfortran)
  _g77version := $(subst ., ,$(shell $(FC) -dumpspecs | awk '/^\*version:/{getline;print $$1}'))
  _g77majorver:= $(word 1,$(_g77version))
  _g77minorver:= $(word 2,$(_g77version))

  deffcomflags = -Wno-line-truncation -fno-second-underscore
  # these flags get rid of unused variables warnings. (-Wno-unused-parameter -Wno-unused-variable)
  deffcomflags += -Wno-unused-parameter -Wno-unused-variable
  deffoptflags = -march=native -O3 -funroll-loops
  deffdbgflags = -march=native -g -ggdb -Wall -W

  # Try this for gcc 4.2 and above...
  # newer versions of gcc (4.2.1 works) can use this
  ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 2 ; echo $$?))
     deffdbgflags += -fbounds-check
  endif

  # Removed because this flag doesn't seem to play well with "valgrind"
  # Try this for gcc 4.3 and above...
  # newer versions of gcc (4.3.2 works) can use this
  # ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 3 ; echo $$?))
  #    deffdbgflags += -fstack-check 
  # endif

  # Note this is a copy/paste from the g77 ifeq above
  # non-default object models are meaningful ONLY on some machines
  ifeq ($(USE_64),TRUE)
    ifeq ($(arch),x86_64)
      # works on an AMD64 machine that defaults to 32bit
      # may also work on a Power{PC,4,5} but I'm not positive
      deffcomflags += -m64
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m64
      ifneq ($(wildcard /usr/local/lib/x86_64), ) # This only helps if it exists!!! -JNJ
        defflibflags += -L/usr/local/lib/x86_64
      endif
    endif
    # on IBM Power
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix64
    endif
  else
    ifeq ($(arch),x86_64)
      # on AMD64 machines with default 64bit, you have to ask
      # for an emulation mode to get 32 bits
      deffcomflags += -m32
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m32
    endif
    # on IBM Power with default 64bit mode
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix32
    endif
  endif

ifeq ($(OPENMPCC),TRUE)
  defflibflags  += -lgfortran -lm -lgomp
else
  defflibflags  += -lgfortran -lm 
endif
endif


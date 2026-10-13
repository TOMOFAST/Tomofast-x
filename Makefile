#
#========================================================================
#
#                          T o m o f a s t - x
#                        -----------------------
#
#           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
#
#               (c) 2021 The University of Western Australia.
#
# The full text of the license is available in file "LICENSE".
#
#========================================================================

################################################################################
# Compiler and linker flags.
################################################################################

# Set to 1 for Windows, 0 for Linux
WINDOWS = 0

# Compiler choice: 1 - GCC, 2 - Intel.
COMPILER = 1

# Intel compiler is the default for Windows. For gfortran (MSYS2) use: make WINDOWS=1 COMPILER=1 USE_MPI=0
# (a value given on the command line overrides the assignments here).
ifeq ($(WINDOWS), 1)
  COMPILER = 2
endif

# 1 when make runs commands in the Windows cmd shell (Intel compiler), 0 when it uses a Unix shell (Linux, MSYS2).
WINSHELL = 0
ifeq ($(WINDOWS), 1)
  ifeq ($(COMPILER), 2)
    WINSHELL = 1
  endif
endif

# Parallel runtime: 1 - MPI (default), 0 - no MPI (one process, uses the MPI stub in src/utils/mpi_stub.F90).
# Example without MPI: make USE_MPI=0
USE_MPI = 1

# Use OpenMP threads: 1 - yes, 0 - no. By default it is on only when MPI is off.
ifeq ($(USE_MPI), 0)
  USE_OPENMP = 1
else
  USE_OPENMP = 0
endif

ifeq ($(COMPILER), 1)
  ifeq ($(USE_MPI), 1)
    # Use MPI Fortran compiler and linker wrappers.
    FC = mpif90
    #FC = /usr/bin/mpif90.openmpi
  else
    FC = gfortran
  endif
else
  ifeq ($(USE_MPI), 1)
    # Intel compiler (mpiifort is deprecated).
    FC = mpiifx
  else
    FC = ifx
  endif
endif

# obj directory
# Use a separate object directory per build configuration, so that changing USE_MPI or USE_OPENMP never reuses stale objects.
OBJDIR = obj
ifeq ($(USE_MPI), 0)
  OBJDIR := $(OBJDIR)_nompi
endif
ifeq ($(USE_OPENMP), 1)
  OBJDIR := $(OBJDIR)_omp
endif

# Executable name.
ifeq ($(WINDOWS), 0)
	EXEC = tomofastx
else
	EXEC = tomofastx.exe
endif

ifeq ($(COMPILER), 1)
  # GNU gfortran pseudo-optimized.
  FFLAGS = -std=f2008 -fconvert=big-endian -O3 -fimplicit-none -frange-check -fmax-errors=10 -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -DUSE_FLUSH6 -J $(OBJDIR)

  # GNU gfortran with debug symbols and full debugging.
  #FFLAGS = -std=f2008 -fconvert=big-endian -O0 -g -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbacktrace -Wunreachable-code -Wunused-label -Wunused-variable -Wimplicit-interface -Wall -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -DUSE_FLUSH6 -J $(OBJDIR)

  # Detect GCC version.
  GFORTRAN_VERSION := $(shell gfortran -dumpversion | cut -f1 -d.)

  # The folder creation in file_utils.F90 differs on Windows.
  ifeq ($(WINDOWS), 1)
    FFLAGS += -DWINDOWS
    # Link the runtime libraries statically, so the .exe runs outside the MSYS2 shell.
    FFLAGS += -static
  endif

  # Add the flag only for GCC 10+.
  # The MPI type mismatches such as in parallel_tools.f90 are intentional (using MPI_IN_PLACE) and safe.
  ifeq ($(shell test $(GFORTRAN_VERSION) -ge 10; echo $$?),0)
    FFLAGS += -fallow-argument-mismatch
  endif

  # GNU flags to output the vector optimizations report.
  OPT_INFO = -ftree-vectorize -fopt-info-vec-optimized=vec.info

  # Adding vectorisation report for GNU compiler.
  #FFLAGS := $(OPT_INFO) $(FFLAGS)
else
  ifeq ($(WINDOWS), 0)
    # Intel compiler optimized for speed for production runs.
    FFLAGS = -convert big_endian -implicitnone -assume buffered_io -assume byterecl -warn truncated_source -warn interfaces -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -DUSE_FLUSH6 -ftz -fpe0 -check nobounds -O3 -xHost -module $(OBJDIR)
  else
    # Windows flags.
    FFLAGS = -fpp -DWINDOWS -DUSE_FLUSH6 -convert:big_endian -warn:declarations -assume:buffered_io -assume:byterecl -warn:truncated_source -warn:interfaces -warn:unused -warn:declarations -warn:alignments -warn:ignore_loc -warn:usage -Qftz -fpe:0 -check:nobounds -O3 -QxHost -module:$(OBJDIR)
  endif

  # Intel compiler with full checking options to debug (slow but very useful to check everything).
  #FFLAGS = -convert big_endian -implicitnone -assume buffered_io -assume byterecl -warn truncated_source -warn interfaces -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -DUSE_FLUSH6 -ftz -fpe0 -check all -debug -g -O0 -traceback -ftrapuv -module $(OBJDIR)
endif

# Build without MPI: use the single-process MPI stub.
ifeq ($(USE_MPI), 0)
  FFLAGS += -DNO_MPI
endif

# OpenMP flags.
ifeq ($(USE_OPENMP), 1)
  ifeq ($(COMPILER), 1)
    FFLAGS += -fopenmp
  else ifeq ($(WINDOWS), 0)
    FFLAGS += -qopenmp
  else
    FFLAGS += -Qopenmp
  endif
endif

# To print a variable run: make print-VARIABLE
print-%  : ; @echo $* = $($*)

default: | $(OBJDIR) ${EXEC}

all: clean default

# EXECUTABLES:

LIBS_FULL =
# For the IPM linking.
#LIBS_FULL =-L/usr/local/lib -lipmf -lipm

# add paths to the search path for files
vpath %.f90 ./src
vpath %.F90 ./src
vpath %.f90 ./src/libs
vpath %.F90 ./src/libs
vpath %.c   ./src/libs
vpath %.f90 ./src/forward
vpath %.F90 ./src/forward
vpath %.f90 ./src/forward/gravmag
vpath %.F90 ./src/forward/gravmag
vpath %.f90 ./src/forward/gravmag/grav
vpath %.F90 ./src/forward/gravmag/grav
vpath %.f90 ./src/forward/gravmag/mag
vpath %.F90 ./src/forward/gravmag/mag
vpath %.f90 ./src/inversion
vpath %.F90 ./src/inversion
vpath %.f90 ./src/utils
vpath %.F90 ./src/utils
vpath %.f90 ./src/tests
vpath %.F90 ./src/tests

# External libraries.
SRC_LIST_LIBS = \
ftnunit.f90

SRC_LIST_UTILS = \
file_utils.F90 \
mpi_tools.F90 \
costs.f90 \
vector.f90 \
string.f90 \
noise.f90 \
paraview.f90 \
parallel_tools.f90 \
memory_tools.f90 \
wavelet_transform.F90 \
sort.f90

# Source files for inversion.
SRC_LIST_INVERSION = \
parameters_inversion.f90 \
sparse_matrix.F90 \
wavelet_utils.F90 \
grid.F90 \
model.F90 \
model_IO.F90 \
inversion_arrays.f90 \
lsqr_solver2.F90 \
damping.F90 \
gradient.F90 \
cross_gradient.F90 \
admm_method.F90 \
clustering.F90 \
damping_gradient.F90 \
joint_inverse_problem.F90

# Source files for gravity and magnetism problems.
SRC_LIST_GRAVMAG = \
parameters_gravmag.f90 \
parameters_grav.f90 \
parameters_mag.f90 \
data_gravmag.f90 \
gravity_field.f90 \
magnetic_field.f90 \
weights_gravmag.f90 \
sensitivity_gravmag.F90

# Main routines.
SRC_LIST_MAIN = \
parameters_init.f90 \
problem_joint_gravmag.F90

# Unit tests.
SRC_LIST_TESTS = \
tests_inversion.f90 \
tests_lsqr.f90 \
tests_parallel_tools.f90 \
tests_sparse_matrix.f90 \
tests_wavelet_compression.f90 \
unit_tests.f90

# Single-process replacement of MPI (must precede global_typedefs.F90).
SRC_LIST_STUB =
ifeq ($(USE_MPI), 0)
  SRC_LIST_STUB = mpi_stub.F90
endif

# Note: inversion files precede the forward problem files,
# to make sure there are no dependencies.
SRC_LIST = \
$(SRC_LIST_STUB) \
global_typedefs.F90 \
$(SRC_LIST_LIBS) \
$(SRC_LIST_UTILS) \
$(SRC_LIST_INVERSION) \
$(SRC_LIST_GRAVMAG) \
$(SRC_LIST_MAIN) \
$(SRC_LIST_TESTS) \
program_tomofastx.F90

# Create object file lists from all source files list.
# First replace .F90 extension.
OBJ_LIST_1 = $(SRC_LIST:%.F90=$(OBJDIR)/%.o)
# Then replace .f90 extension.
OBJ_LIST_2 = $(OBJ_LIST_1:%.f90=$(OBJDIR)/%.o)
# Then replace .c extension.
OBJ_LIST = $(OBJ_LIST_2:%.c=$(OBJDIR)/%.o)

# Targets to create the object directory if it does not exist.
$(OBJ_LIST): | $(OBJDIR)

$(OBJDIR):
ifeq ($(WINSHELL), 0)
	@test -d $(OBJDIR) || (rm -f $(OBJDIR); mkdir -p $(OBJDIR))
else
	@if not exist $(OBJDIR) mkdir $(OBJDIR)
endif

# Implicit rules for object files.
$(OBJDIR)/%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

# Target to build the actual executable.
${EXEC}: $(OBJ_LIST)
	 $(FC) $(FFLAGS) -o ${EXEC} $(OBJ_LIST) $(LIBS_FULL) $(LDFLAGS)

clean:
ifeq ($(WINSHELL), 0)
	rm -rf *.o *.mod ${EXEC} obj obj_nompi obj_nompi_omp obj_omp
else
	-del /Q *.o *.mod ${EXEC} 2>nul
	-rmdir /S /Q obj 2>nul
	-rmdir /S /Q obj_nompi 2>nul
	-rmdir /S /Q obj_nompi_omp 2>nul
	-rmdir /S /Q obj_omp 2>nul
endif

# targets to clean up all object directories
purge: clean
realclean: purge
deepclean: purge
mrproper:  purge


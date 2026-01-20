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

# Always use Intel compiler for Windows.
ifeq ($(WINDOWS), 1)
  COMPILER = 2
endif

# Use MPI Fortran compiler and linker wrappers.
ifeq ($(COMPILER), 1)
  FC = mpif90
  #FC = /usr/bin/mpif90.openmpi
else
  # Intel compiler (mpiifort is deprecated).
  FC = mpiifx
endif

# obj directory
OBJDIR = obj

# Executable name.
ifeq ($(WINDOWS), 0)
	EXEC = tomofastx
else
	EXEC = tomofastx.exe
endif

ifeq ($(COMPILER), 1)
  # GNU gfortran pseudo-optimized.
  # Note: Added -fallow-argument-mismatch and removed -pedantic as the MPI type mismatches in parallel_tools.f90 are intentional (using MPI_IN_PLACE) and safe.
  FFLAGS = -std=f2008 -fconvert=big-endian -O3 -fimplicit-none -frange-check -fmax-errors=10 -fallow-argument-mismatch -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -DUSE_FLUSH6 -J $(OBJDIR)

  # GNU gfortran with debug symbols and full debugging.
  #FFLAGS = -std=f2008 -fconvert=big-endian -O0 -g -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbacktrace -Wunreachable-code -Wunused-label -Wunused-variable -Wimplicit-interface -Wall -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -DUSE_FLUSH6 -J $(OBJDIR)

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
sparse_matrix.f90 \
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

# Note: inversion files precede the forward problem files,
# to make sure there are no dependencies.
SRC_LIST = \
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
ifeq ($(WINDOWS), 0)
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
	 $(FC) $(FFLAGS) -o ${EXEC} $(OBJ_LIST) $(LIBS_FULL)

clean:
ifeq ($(WINDOWS), 0)
	rm -rf *.o *.mod ${EXEC} $(OBJDIR)
else
	-del /Q *.o *.mod ${EXEC} 2>nul
	-rmdir /S /Q $(OBJDIR) 2>nul
endif

# targets to clean up all object directories
purge: clean
realclean: purge
deepclean: purge
mrproper:  purge


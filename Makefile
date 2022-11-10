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
# Configuration section.
################################################################################

# Disables writing output files for Paraview visualization, and some other output.
# Note: disable this for performance measurements since the output is performed by
# the master CPU only, and thus slows down the parallel performance.
SUPPRESS_OUTPUT = NO

# TODO: Move here compiler choice (GNU/Intel) and optimization level.

# By default, keep the output enabled.
ifeq ($(strip $(SUPPRESS_OUTPUT)),)
SUPPRESS_OUTPUT = NO
endif

################################################################################
# Compiler and linker flags.
################################################################################

# NOTE: To set explicitly gfortran-4.9 compiler, execute
# export OMPI_FC=gfortran-4.9
# or put the above line into ~/.bashrc file.

# NOTE: When mpiifort is not found (sometimes happens), compile using mpif90 with
# export OMPI_FC=ifort
# export OMPI_CC=icc

# NOTE: To check what compiler is used by mpif90, execute
# mpif90 -v

# Use MPI Fortran and C compiler and linker wrappers.
#FC = mpiifort
FC = mpif90
CC = mpicc
#CC = mpiicc

# obj directory
OBJDIR = obj

# Executable name.
EXEC = tomofastx

# Intel ifort with full checking options to debug (slow but very useful to check everything).
#FFLAGS = -convert big_endian -vec-report0  -implicitnone -assume buffered_io -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -DUSE_FLUSH6 -ftz -fpe0 -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -module $(OBJDIR)

# Intel ifort optimized for speed for production runs (add -mcmodel=medium -shared-intel to use more than 2 GB of memory).
#FFLAGS = -convert big_endian -vec-report0 -implicitnone -assume buffered_io -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -DUSE_FLUSH6 -ftz -fpe0 -check nobounds -O3 -xHost -module $(OBJDIR)

# GNU gfortran pseudo-optimized.
FFLAGS = -std=f2008 -fconvert=big-endian -O3 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -DUSE_FLUSH6 -J $(OBJDIR)

# GNU gfortran with debug symbols and full debugging.
#FFLAGS = -std=f2008 -fconvert=big-endian -O0 -g -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Warray-temporaries -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbacktrace -Wunreachable-code -Wunused-label -Wunused-variable -Wimplicit-interface -Wall -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -DUSE_FLUSH6 -J $(OBJDIR)

# (GNU) Flags to output the vector optimizations report.
# possible flags: -fopt-info-vec-optimized, -fopt-info-vec-missed, -fopt-info-vec-note, -fopt-info-vec-all
OPT_INFO = -ftree-vectorize -fopt-info-vec-optimized=vec.info

# Comment this for non GNU compiler.
#FFLAGS := $(OPT_INFO) $(FFLAGS)

CFLAGS =

ifeq ($(strip $(SUPPRESS_OUTPUT)),YES)
FFLAGS := $(FFLAGS) -DSUPPRESS_OUTPUT
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
mpi_tools.F90 \
costs.f90 \
vector.f90 \
string.f90 \
noise.f90 \
paraview.f90 \
parallel_tools.f90 \
wavelet_transform.F90 \
sort.f90

# Source files for inversion.
SRC_LIST_INVERSION = \
sparse_matrix.f90 \
grid.F90 \
model.F90 \
model_IO.F90 \
parameters_inversion.f90 \
inversion_arrays.f90 \
lsqr_solver2.F90 \
sca_solver.F90 \
damping.F90 \
sensitivity_matrix.F90 \
gradient.F90 \
cross_gradient.F90 \
method_of_weights.F90 \
admm_method.F90 \
clustering.F90 \
damping_gradient.F90 \
inverse_problem.F90 \
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
tests_method_of_weights.f90 \
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
	@test -d $(OBJDIR) || (rm -f $(OBJDIR); mkdir -p $(OBJDIR))

# Implicit rules for object files.
$(OBJDIR)/%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<


# Target to build the actual executable.
${EXEC}: $(OBJ_LIST)
	 $(FC) $(FFLAGS) -o ${EXEC} $(OBJ_LIST) $(LIBS_FULL)


clean:
	rm -rf *.o *.mod ${EXEC} $(OBJDIR)

# targets to clean up all object directories
purge: clean
realclean: purge
deepclean: purge
mrproper:  purge


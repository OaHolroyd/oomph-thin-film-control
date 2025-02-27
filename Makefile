# ============================================================================ #
#   VARIABLE DEFINITIONS                                                       #
# ============================================================================ #
# compiler/linker
CC=mpicc
CXX=mpic++
LD=$(CXX)

# if the environment variable OOMPHLIB is not set, uncomment this and set it
OOMPHLIB=/storage/maths/mauvvq/oomph-lib

# flags (for MacOS)
#WARNINGS=-Wall -Wextra -pedantic -Wno-implicit-function-declaration -Wno-unused-parameter
#CFLAGS=-O3 -Wall
#CXXFLAGS=-O3 -Wall -std=c++17 -DgFortran
#LDFLAGS=$(CXXFLAGS)
#LDLIBS=-L$(OOMPHLIB)/build/lib $(OOMPHLIB)/build/lib/*.a \
#       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/aarch64-apple-darwin23/14 \
#       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current \
#       -lemutls_w -lheapt_w -lgfortran -lquadmath -lm -llapacke -L/opt/homebrew/opt/lapack/lib
#INCLUDES=-isystem $(OOMPHLIB)/build/include \
#         -isystem $(OOMPHLIB) \
#         -I/opt/homebrew/opt/lapack/include
#OLFLAGS=-DHAVE_CONFIG_H -DOOMPH_HAS_STACKTRACE -DOOMPH_HAS_UNISTDH \
#        -DOOMPH_HAS_TRIANGLE_LIB -DOOMPH_HAS_TETGEN_LIB -DUSING_OOMPH_SUPERLU \
#        -DUSING_OOMPH_SUPERLU_DIST

# flags (for Taskfarm)
WARNINGS=-Wall -Wextra -pedantic -Wno-implicit-function-declaration -Wno-unused-parameter
CFLAGS=-O3 -Wall
CXXFLAGS=-O3 -Wall -std=c++17 -DgFortran -fopenmp
LDFLAGS=$(CXXFLAGS)

LDLIBS=-L$(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib \
       -L$(OOMPHLIB)/external_distributions/boost/boost_default_installation/lib \
       -L$(OOMPHLIB)/external_distributions/gmp/gmp_default_installation/lib \
       -L$(OOMPHLIB)/external_distributions/mpfr/mpfr_default_installation/lib \
       -L$(OOMPHLIB)/external_distributions/cgal/cgal_default_installation/lib \
       -L$(OOMPHLIB)/build/lib $(OOMPHLIB)/build/lib/*.so \
       -L/software/easybuild/software/UCC/1.0.0-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/PMIx/4.1.2-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/libfabric/1.15.1-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/UCX/1.12.1-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/libevent/2.1.12-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/hwloc/2.7.1-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/zlib/1.2.12-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/Perl/5.34.1-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/Perl/5.34.1-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/pkgconf/1.8.0-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/GCCcore/11.3.0/lib64 \
       -L/software/easybuild/software/GCCcore/11.3.0/lib \
       -L/software/easybuild/software/binutils/2.38-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/libpciaccess/0.16-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/numactl/2.0.14-GCCcore-11.3.0/lib64 \
       -L/software/easybuild/software/OpenSSL/1.1/lib64 \
       $(OOMPHLIB)/build/lib/libpoisson.so \
       $(OOMPHLIB)/build/lib/libgeneric.so \
       -ldmumps \
       -lmumps_common \
       -lboost_thread \
       -lboost_system \
       $(OOMPHLIB)/external_distributions/mpfr/mpfr_default_installation/lib/libmpfr.so \
       $(OOMPHLIB)/external_distributions/gmp/gmp_default_installation/lib/libgmp.so \
       -lCGAL_Core \
       -lCGAL \
       $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib/libscalapack.a \
       $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib/blacs.a \
       $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib/blacsF77.a \
       $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib/blacs_copy.a \
       $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib/libpord.a \
       $(OOMPHLIB)/build/lib/liboomph_hsl.so \
       $(OOMPHLIB)/build/lib/liboomph_crbond_bessel.so \
       $(OOMPHLIB)/build/lib/liboomph_triangle.a \
       $(OOMPHLIB)/build/lib/liboomph_tetgen.so \
       /software/easybuild/software/GCCcore/11.3.0/lib/../lib64/libstdc++.so \
       $(OOMPHLIB)/build/lib/liboomph_superlu_4.3.a \
       $(OOMPHLIB)/build/lib/liboomph_parmetis_3.1.1.a \
       $(OOMPHLIB)/build/lib/liboomph_superlu_dist_3.0.a \
       $(OOMPHLIB)/build/lib/liboomph_metis_from_parmetis_3.1.1.a \
       /software/easybuild/software/OpenBLAS/0.3.20-GCC-11.3.0/lib/libopenblas.a \
       -L/software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib \
       -L/software/easybuild/software/hwloc/2.7.1-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/libevent/2.1.12-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/libarchive/3.6.1-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/cURL/7.83.0-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/bzip2/1.0.8-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/ncurses/6.3-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/UCC/1.0.0-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/PMIx/4.1.2-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/libfabric/1.15.1-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/UCX/1.12.1-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/libevent/2.1.12-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/OpenSSL/1.1/lib/../lib64 \
       -L/software/easybuild/software/hwloc/2.7.1-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/libpciaccess/0.16-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/libxml2/2.9.13-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/XZ/5.2.5-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/numactl/2.0.14-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/OpenBLAS/0.3.20-GCC-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/binutils/2.38-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/zlib/1.2.12-GCCcore-11.3.0/lib/../lib64 \
       -L/software/easybuild/software/GCCcore/11.3.0/lib/gcc/x86_64-pc-linux-gnu/11.3.0 \
       -L/software/easybuild/software/GCCcore/11.3.0/lib/gcc/x86_64-pc-linux-gnu/11.3.0/../../../../lib64 \
       -L/lib/../lib64 \
       -L/usr/lib/../lib64 \
       -L/software/easybuild/software/libarchive/3.6.1-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/cURL/7.83.0-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/bzip2/1.0.8-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/ncurses/6.3-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/UCC/1.0.0-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/PMIx/4.1.2-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/libfabric/1.15.1-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/UCX/1.12.1-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/OpenSSL/1.1/lib \
       -L/software/easybuild/software/libpciaccess/0.16-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/libxml2/2.9.13-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/XZ/5.2.5-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/numactl/2.0.14-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/OpenBLAS/0.3.20-GCC-11.3.0/lib \
       -L/software/easybuild/software/binutils/2.38-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/zlib/1.2.12-GCCcore-11.3.0/lib \
       -L/software/easybuild/software/GCCcore/11.3.0/lib/gcc/x86_64-pc-linux-gnu/11.3.0/../../.. \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libmpi_usempif08.so \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libmpi_usempi_ignore_tkr.so \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libmpi_mpifh.so \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libmpi.so \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libopen-rte.so \
       /software/easybuild/software/OpenMPI/4.1.4-GCC-11.3.0/lib/libopen-pal.so \
       -lz \
       /software/easybuild/software/hwloc/2.7.1-GCCcore-11.3.0/lib/libhwloc.so \
       /software/easybuild/software/libpciaccess/0.16-GCCcore-11.3.0/lib/libpciaccess.so \
       /software/easybuild/software/libevent/2.1.12-GCCcore-11.3.0/lib/libevent_core.so \
       /software/easybuild/software/libevent/2.1.12-GCCcore-11.3.0/lib/libevent_pthreads.so \
       /software/easybuild/software/GCCcore/11.3.0/lib/../lib64/libgfortran.so \
       /software/easybuild/software/GCCcore/11.3.0/lib/../lib64/libgomp.so \
       -ldl \
       /software/easybuild/software/GCCcore/11.3.0/lib/../lib64/libquadmath.so \
       -lm \
       -lpthread \
       -fopenmp \
       -llapacke

INCLUDES=-isystem $(OOMPHLIB)/build/include \
         -isystem $(OOMPHLIB) \
         -isystem $(OOMPHLIB)/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/include \
         -DOOMPH_HAS_MUMPS \
         -isystem $(OOMPHLIB)/external_distributions/boost/boost_default_installation/include \
         -DOOMPH_HAS_BOOST \
         -isystem $(OOMPHLIB)/external_distributions/gmp/gmp_default_installation/include \
         -DOOMPH_HAS_GMP \
         -isystem $(OOMPHLIB)/external_distributions/mpfr/mpfr_default_installation/include \
         -DOOMPH_HAS_MPFR \
         -isystem $(OOMPHLIB)/external_distributions/cgal/cgal_default_installation/include \
         -DOOMPH_HAS_CGAL
OLFLAGS=-DHAVE_CONFIG_H -DOOMPH_HAS_STACKTRACE -DOOMPH_HAS_UNISTDH -DOOMPH_HAS_MPI \
        -DOOMPH_HAS_TRIANGLE_LIB -DOOMPH_HAS_TETGEN_LIB -DUSING_OOMPH_SUPERLU \
        -DUSING_OOMPH_SUPERLU_DIST -DOOMPH_HAS_FPUCONTROLH -DOOMPH_HAS_MALLOCH

# executable
EXE=main

# directories
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./output

# files
SRC=$(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.c)
TMP_OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.cpp=.o)))
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(TMP_OBJ:.c=.o)))
DEPS=$(patsubst %.o,%.d,$(OBJ)) # dependency files


# ============================================================================ #
#   RULES                                                                      #
# ============================================================================ #
# link objects into single binary
$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(LDLIBS)

# compile c++ to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OLFLAGS) $(INCLUDES) $(WARNINGS) -MMD -MP -c -o $@ $<

# compile c to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c Makefile | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(WARNINGS) -MMD -MP -c -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# include dependency information
-include $(DEPS)

# remove build files and executable
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(EXE) $(OUT_DIR)/*.dat $(OUT_DIR)/*.png $(OUT_DIR)/*.gif

.PHONY: all
all: $(EXE)

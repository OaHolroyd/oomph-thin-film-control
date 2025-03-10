# ============================================================================ #
#   VARIABLE DEFINITIONS                                                       #
# ============================================================================ #
# compiler/linker
CC=gcc
CXX=g++
LD=$(CXX)

# flags (this assumes that the environment variable OOMPHLIB contains the root directory of the oomph-lib library)
WARNINGS=-Wall -Wextra -pedantic -Wno-implicit-function-declaration -Wno-unused-parameter
CFLAGS=-O3 -Wall
CXXFLAGS=-O3 -Wall -std=c++17 -DgFortran
LDFLAGS=$(CXXFLAGS)
LDLIBS=-L$(OOMPHLIB)/build/lib $(OOMPHLIB)/build/lib/*.a \
       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/aarch64-apple-darwin23/14 \
       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current \
       -lemutls_w -lheapt_w -lgfortran -lquadmath -lm -llapacke -L/opt/homebrew/opt/lapack/lib
INCLUDES=-isystem $(OOMPHLIB)/build/include \
         -isystem $(OOMPHLIB) \
         -I/opt/homebrew/opt/lapack/include
OLFLAGS=-DHAVE_CONFIG_H -DOOMPH_HAS_STACKTRACE -DOOMPH_HAS_UNISTDH \
        -DOOMPH_HAS_TRIANGLE_LIB -DOOMPH_HAS_TETGEN_LIB -DUSING_OOMPH_SUPERLU \
        -DUSING_OOMPH_SUPERLU_DIST

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
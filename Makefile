# ============================================================================ #
#   VARIABLE DEFINITIONS                                                       #
# ============================================================================ #
# compiler/linker
CC=g++
LD=$(CC)

# flags (this assumes that the environment variable OOMPHLIB contains the root directory of the oomph-lib library)
WARNINGS=-Wall -Wextra -pedantic -Wno-implicit-function-declaration
CFLAGS=-O3 -Wall -std=c++17 -DgFortran
LDFLAGS=$(CFLAGS)
LDLIBS=-L$(OOMPHLIB)/build/lib $(OOMPHLIB)/build/lib/*.a \
       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/aarch64-apple-darwin23/14 \
       -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current \
       -lemutls_w -lheapt_w -lgfortran -lquadmath
INCLUDES=-I$(OOMPHLIB)/build/include \
         -I$(OOMPHLIB)
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
SRC=$(wildcard $(SRC_DIR)/*.cpp)
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.cpp=.o)))
DEPS=$(patsubst %.o,%.d,$(OBJ)) # dependency files


# ============================================================================ #
#   RULES                                                                      #
# ============================================================================ #
# link objects into single binary
$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(LDLIBS)

# compile object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(OLFLAGS) $(INCLUDES) $(WARNINGS) -MMD -MP -c -o $@ $<

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
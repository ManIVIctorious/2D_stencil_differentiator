# Compiler
  CC = gcc
# List of compiler flags
 #CFLAGS += -g                      # Enable debug symbols
  CFLAGS += -O2                     # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native           # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native           # Tune for current chipset, remain backwards compatible

  CFLAGS += -Werror                 # Treat warnings as errors
  CFLAGS += -Wall                   # Enable base set of warnings
  CFLAGS += -Wextra                 # Enable additional warnings
  CFLAGS += -Wstrict-prototypes     # Enable strict-prototypes warning
  CFLAGS += -Wmissing-prototypes    # Enable missing-prototypes warning
  CFLAGS += -Wno-sign-compare       # Disable sign-compare warning
 #CFLAGS = -g -w                    # Disable all warnings

# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`


# Resulting executable
  EXE = bin/2d_differenciation

# List of resulting object files
  OBJ += main.o
  OBJ += InvertMatrix.o
  OBJ += Factorial.o
  OBJ += FiniteDifferenceCoefficients.o
  OBJ += InputFunction.o

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ) Makefile
	$(CC) $(CFLAGS) $(LIB) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)

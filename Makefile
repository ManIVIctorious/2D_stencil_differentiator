# Compiler
  CC = gcc
# List of compiler flags
  CFLAGS += -O2 -Wall -Wextra -Werror -march=native

# Resulting executable
  EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE = $(EXEDIR)/2d-stencil-differentiator


# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`

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

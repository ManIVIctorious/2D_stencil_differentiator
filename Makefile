CC = gcc
CFLAGS = -g -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Werror -Wno-sign-compare -mtune=native

EXE = bin/2D_differenciation
OBJ = main.o InvertMatrix.o Factorial.o FiniteDifferenceCoefficients.o InputFunction.o
LIB = `pkg-config --cflags --libs gsl`

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(LIB) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)

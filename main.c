
#include <stdio.h>
#include <stdlib.h>


int FiniteDifferenceCoefficients(unsigned int derivative, unsigned int point_number, int* point_location, double* result_vector);
int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2);

int main(int argc, char* argv[]){

//-------------------------------------------------------------------------------------------
// Declaration  Declaration  Declaration  Declaration  Declaration  Declaration  Declaration
//-------------------------------------------------------------------------------------------

    int      i, j, k, l;
    long int sx, sy;
    int      nx, ny, control;
    double   dx, dy;
    double   spacing_threshold = 1.0E-12; // abs(q[i] - q[i+1])

    char   * endptr    = NULL;
    char   * inputfile = NULL;
    FILE   * fd        = NULL;

    int    * Location_x     = NULL; // freed
    int    * Location_y     = NULL; // freed
    double * first_deriv_x  = NULL; // freed
    double * first_deriv_y  = NULL; // freed

    double * second_deriv_x = NULL; // freed
    double * second_deriv_y = NULL; // freed
    double * xy_cross_deriv = NULL; // freed

    double * x = NULL; // freed
    double * y = NULL; // freed
    double * V = NULL; // freed
    double * d2dxdxV = NULL; // freed
    double * d2dydyV = NULL; // freed
    double * d2dxdyV = NULL; // freed


//-------------------------------------------------------------------------------------------
// Initialisation Initialisation Initialisation Initialisation Initialisation Initialisation
//-------------------------------------------------------------------------------------------
    fd = stdout;

    sx = strtol(argv[1], &endptr, 10);
    if(endptr == argv[1]){
        fprintf(stderr, "\n (-) Input Error:");
        fprintf(stderr, "\n     Stencil x-dimension could not be converted to integer");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    sy = strtol(argv[2], &endptr, 10);
    if(endptr == argv[2]){
        fprintf(stderr, "\n (-) Input Error:");
        fprintf(stderr, "\n     Stencil y-dimension could not be converted to integer");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    if(sx % 2 == 0 || sy % 2 == 0){
        fprintf(stderr, "\n (-) Input Error:");
        fprintf(stderr, "\n     Please use an odd number to specify the stencil size");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

    inputfile = argv[3];
    if(inputfile == NULL){
        fprintf(stderr, "\n (-) Input Error:");
        fprintf(stderr, "\n     Please specify a valid input file");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

//-------------------------------------------------------------------------------------------
//   Generate stencils    Generate stencils    Generate stencils    Generate stencils
//-------------------------------------------------------------------------------------------
// fill PointLocation array with subsequent point locations
    Location_x = calloc(sx, sizeof(int));
    Location_y = calloc(sy, sizeof(int));

    for(i = 0; i < sx; ++i) Location_x[i] = i - (sx/2);
    for(i = 0; i < sy; ++i) Location_y[i] = i - (sy/2);

// calculate finite difference coefficients for first order differenciations
    first_deriv_x  = calloc(sx, sizeof(double));
    first_deriv_y  = calloc(sy, sizeof(double));

    FiniteDifferenceCoefficients(1, sx, Location_x, first_deriv_x);
    FiniteDifferenceCoefficients(1, sy, Location_y, first_deriv_y);

// calculate finite difference coefficients for second order differenciations
    second_deriv_x = calloc(sx, sizeof(double));
    second_deriv_y = calloc(sy, sizeof(double));

    FiniteDifferenceCoefficients(2, sx, Location_x, second_deriv_x);
    FiniteDifferenceCoefficients(2, sy, Location_y, second_deriv_y);

// free memory
    free(Location_x); Location_x = NULL;
    free(Location_y); Location_y = NULL;

// generate stencil for 2D cross derivative out of 1D derivatives (tensor product)
    xy_cross_deriv = calloc(sx * sy, sizeof(double));

    for(i = 0; i < sy; ++i){
        for(j = 0; j < sx; ++j){
            xy_cross_deriv[i*sx + j] = first_deriv_y[i] * first_deriv_x[j];
        }
    }

//-------------------------------------------------------------------------------------------
//  Input   Input   Input   Input   Input   Input    Input   Input   Input   Input   Input
//-------------------------------------------------------------------------------------------
    x = malloc(sizeof(double));
    y = malloc(sizeof(double));
    V = malloc(sizeof(double));
    control = InputFunction(inputfile, &x, &y, &V, &nx, &ny);
    if(x == NULL || y == NULL || V  == NULL){
        fprintf(stderr, "\n (-) Error in reading data from input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n     Memory allocation for x, y or V");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

// check if the "N nx ny" line in input file matches the number of data points
    if(control != nx*ny){
        fprintf(stderr, "\n (-) Error in reading data from input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d*%d\"", nx*ny, nx, ny);
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(2);
    }

// there must be at least as many data points as the stencil size
    if(nx <= sx || ny <= sy){
        fprintf(stderr, "\n (-) Error in reading data from input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n     Insufficient number of data points %dx%d for stencil size %dx%d.", nx, ny, sx, sy);
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(2);
    }

//------------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------------
    dx = x[ny] - x[0];
    dy = y[1]  - y[0];
    for(i = 1; i < (nx*ny-ny); ++i){

    // spacing between x[i] and x[i-1]
        if( (dx - (x[i+ny] - x[i]))*(dx - (x[i+ny] - x[i])) > spacing_threshold*spacing_threshold ){
            fprintf(stderr, "\n (-) Error in input-file: \"%s\"", inputfile);
            fprintf(stderr, "\n     Coordinate spacing not equivalent.");
            fprintf(stderr, "\n     Aborting - please check your input...\n\n");
            exit(3);
        }
        dx = x[i+ny] - x[i];

    // spacing between x[i] and y[i-1]
        if( (i+1)%ny != 0 ){
            if( (dy - (y[i+1] - y[i]))*(dy - (y[i+1] - y[i])) > spacing_threshold*spacing_threshold ){
                fprintf(stderr, "\n (-) Error in input-file: \"%s\"", inputfile);
                fprintf(stderr, "\n     Coordinate spacing not equivalent.");
                fprintf(stderr, "\n     Aborting - please check your input...\n\n");
                exit(3);
            }
            dy = y[i+1] - y[i];
        }
    }

//-------------------------------------------------------------------------------------------
//  Stencil application    Stencil application    Stencil application    Stencil application
//-------------------------------------------------------------------------------------------
// d^2/dx^2
    d2dxdxV = calloc(nx*ny, sizeof(double));

    // start with the first block without boundary and end before the last
    for(i = (sx/2)*ny; i < nx*ny-(sx/2)*ny; ++i){
        // iteration over stencil, V is read block wise
        for(j = -(sx/2); j < (sx/2 + 1); ++j){
            d2dxdxV[i] += V[i + j*ny] * second_deriv_x[(sx/2) + j] / dx / dx;
        }
    }

// d^2/dy^2
    d2dydyV = calloc(nx*ny, sizeof(double));

    // calculate each nx block on its own
    for(i = 0; i < nx; ++i){
        // all lines inside of a block except of the first and last (sy/2)
        for(j = (sy/2); j < ny-(sy/2); ++j){
            // iteration over stencil
            for(k = -(sy/2); k < (sy/2 + 1); ++k){
                d2dydyV[i*ny + j] += V[(i*ny + j) + k] * second_deriv_y[(sy/2) + k] / dy / dy;
            }
        }
    }

// cross derivative
    d2dxdyV = calloc(nx*ny, sizeof(double));

    // i*ny + j selects the points with sufficient distance from boundaries
    for(i = (sx/2); i < nx - (sx/2); ++i){
        for(j = (sy/2); j < ny - (sy/2); ++j){
            // k and l apply the indices for potential and stencil
            for(k = -(sy/2); k < (sy/2 + 1); ++k){
                for(l = -(sx/2); l < (sx/2 + 1); ++l){
                    d2dxdyV[i*ny + j] += V[(i*ny + j) + k + l*ny] * xy_cross_deriv[(k + sy/2)*sx + (l + sx/2)] / dx / dy;
                }
            }
        }
    }


//-------------------------------------------------------------------------------------------
//  Output   Output   Output   Output   Output   Output   Output   Output   Output   Output
//-------------------------------------------------------------------------------------------
// print stencils
    fprintf(fd, "# Stencils for Hessian matrix with %dx%d points\n#", sx, sy);
//  positions
    fprintf(fd, "\n# --==> Position map  (x,y) <==--\n#");
    for(i=0; i<sy; ++i){
        for(j=0; j<sx; ++j){
            fprintf(fd, "\t(% d,% d)", j-(sy/2), i-(sx/2));
        }
        fprintf(fd, "\n#");
    }
    fprintf(fd, "\n");

//  stencil for d^2/dx^2
    fprintf(fd, "# --==> d^2/dx^2 <==--\n#");
    for(i = 0; i < sx; ++i){
        fprintf(fd, "\t% 4d", i);
        fprintf(fd, "\t% lf", second_deriv_x[i]);
    }
    fprintf(fd, "\n#\n");

// free memory
    free(second_deriv_x); second_deriv_x = NULL;

//  stencil for d^2/dy^2
    fprintf(fd, "# --==> d^2/dy^2 <==--\n#");
    for(i = 0; i < sy; ++i){
        fprintf(fd, "\t% 4d", i);
        fprintf(fd, "\t% lf\n#", second_deriv_y[i]);
    }
    fprintf(fd, "\n");

// free memory
    free(second_deriv_y); second_deriv_y = NULL;

//  stencil for d^2/(dxdy)
    fprintf(fd, "# --==> d^2/(dxdy) <==--\n#");

    fprintf(fd, "\td/dy / d/dx|");
    for(j=0; j<sx; ++j){
        fprintf(fd, "\t\t% lf", first_deriv_x[j]);
    }
    fprintf(fd, "\n#");
    fprintf(fd, "\t");
    for(j=0; j<sx; ++j){
        fprintf(fd, "--------------------------");
    }
    fprintf(fd, "\n#");

    for(i=0; i<sy; ++i){
        fprintf(fd, "\t% lf  |", first_deriv_y[i]);
        for(j=0; j<sx; ++j){
            fprintf(fd, "\t% 4d",  sx*i + j);
            fprintf(fd, "\t% lf", xy_cross_deriv[sx*i + j]);
        }
        fprintf(fd, "\n#");
    }
    fprintf(fd, "\n");

// free memory
    free(first_deriv_x);  first_deriv_x  = NULL;
    free(first_deriv_y);  first_deriv_y  = NULL;
    free(xy_cross_deriv); xy_cross_deriv = NULL;

// output x, y, V and the 3 derivatives
    fprintf(fd, "# nx = % 4d, ny = % 4d\n", nx, ny);
    for(i = 0; i < nx*ny; ++i){
        if(i % ny == 0) fprintf(fd, "\n");
        fprintf(fd, "\t% 4d", i         );
        fprintf(fd, "\t% lf", x[i]      );
        fprintf(fd, "\t% lf", y[i]      );
        fprintf(fd, "\t% lf", V[i]      );
        fprintf(fd, "\t% lf", d2dxdxV[i]);
        fprintf(fd, "\t% lf", d2dydyV[i]);
        fprintf(fd, "\t% lf", d2dxdyV[i]);
        fprintf(fd, "\n");
    }

// free memory
    free(x); x = NULL;
    free(y); y = NULL;
    free(V); V = NULL;
    free(d2dxdxV); d2dxdxV = NULL;
    free(d2dydyV); d2dydyV = NULL;
    free(d2dxdyV); d2dxdyV = NULL;

    return 0;
}

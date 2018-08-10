
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

// Dependencies
int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);
int factorial(int n);

// Offered prototypes
int FiniteDifferenceCoefficients(unsigned int derivative, unsigned int point_number, int* point_location, double* result_vector);


// get finite difference coefficients for arbitrary functions with equispaced grid:
//
// derivative:      0 for f(x), 1 for f'(x), 2 for f''(x), etc.
// point number:    Number of input function points taken into account for derivative,
//                    e.g. 3 for f(x-h) f(x) f(x+h)
// point location:  Integer array, containing positioning values for h, 
//                    e.g. -1 0 1     for f(x-h)  f(x) f(x+h),
//                    or   -2 0 1 3 5 for f(x-2h) f(x) f(x+h) f(x+3h) f(x+5h)
// result vector:   Array containing the desired finite difference coefficients
int FiniteDifferenceCoefficients(unsigned int derivative, unsigned int point_number, int* point_location, double* result_vector){

    unsigned int i, j;
    double     * constraint_vector;
    gsl_matrix *     coefficient_matrix;
    gsl_matrix * inv_coefficient_matrix;

// check if enough points for desired derivative are given
    if(derivative >= point_number){
        fprintf(stderr, "Please enter a derivative order that is less than the number of points in your stencil\n");
        fprintf(stderr, "derivative order:\t%u\n", derivative);
        fprintf(stderr, "sample point number:\t%u\n", point_number);
        exit(1);
    }

// allocate memory for coefficient_matrix and its inverse counterpart
    coefficient_matrix     = gsl_matrix_alloc(point_number, point_number);
    inv_coefficient_matrix = gsl_matrix_alloc(point_number, point_number);
    if(coefficient_matrix == NULL || inv_coefficient_matrix == NULL){
        fprintf(stderr, "Error in memory allocation for coefficient_matrix\n");
        exit(1);
    }

// allocate memory for constraint and set the entries to zero
    constraint_vector = calloc(point_number, sizeof(double));
    if(constraint_vector == NULL){
        fprintf(stderr, "Error in memory allocation for constraint_vector\n");
        exit(1);
    }

// fill the constraint_vector with data points
//  all derivatives of f(x) have to cancel out, except of the desired one
//  f^(n)(x) / n! = 1  -=>  f^(n)(x) = n!
    constraint_vector[derivative] = (double) factorial(derivative);

// fill the coefficient matrix with data points
//  first  line c1^0 c2^0 c3^0 c4^0 ... cn^0
//  second line c1^1 c2^1 c3^1 c4^1 ... cn^1
//      ...     ...  ...  ...  ...  ... ... 
//  n^th   line c1^n c2^n c3^n c4^n ... cn^n
    for(i = 0; i < point_number; ++i){
        for(j = 0; j < point_number; ++j){
            gsl_matrix_set(coefficient_matrix, i, j, pow(point_location[j],i));
        }
    }

// solve coefficient_matrix * b = constraint_vector
//  inversion and matrix multiplication
    InvertMatrix(coefficient_matrix, inv_coefficient_matrix, point_number);
    for(i = 0; i < point_number; ++i){
        for(j = 0;j < point_number; ++j){
            result_vector[i] += gsl_matrix_get(inv_coefficient_matrix, i, j) * constraint_vector[j];
        }
    }

// free unused memory
    gsl_matrix_free(coefficient_matrix);        coefficient_matrix      = NULL;
    gsl_matrix_free(inv_coefficient_matrix);    inv_coefficient_matrix  = NULL;
    free(constraint_vector);                    constraint_vector       = NULL;


// print finite difference coefficients in human readable form
//    for(i = 0; i < point_number; ++i){
//        fprintf(stdout, "#% 14.6lf f(x%+0dh) / (h**%u)\n", result_vector[i], point_location[i], derivative);
//    }
    return 0;
}

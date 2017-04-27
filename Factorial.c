
#include <stdio.h>
#include <stdlib.h>

// Offered prototypes
int factorial(int n);
 
// calculate factorial of arbitrary not negative integer
int factorial(int n){

    int i;
    int factorial = 1;

    if(n < 0){
        fprintf(stderr, "Factorial is only defined for non negative integers\n");
        fprintf(stderr, "Please check your input, aborting...\n");
        exit(1);
    }

    for(i = 2; i <= n; ++i){
        factorial *= i;
    }

    return factorial;
}

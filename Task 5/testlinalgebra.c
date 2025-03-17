#include <stdio.h>
#include "linalgebra.h"

int main()
{
    double a[9] = {1,1,1,1,2,2,2,3,4};
    double b[3] = {3,5,9};
    double x[3]; 

    printmat(a, 3);

    gausselim(a, 3, b);
    
    printmat(a, 3);
    for(int i =0; i<3; ++i) printf("%e", b[i]);
    printf("\n");

    bks(a,3,b,x);

    printf("x: ");
    for(int i =0; i<3; ++i) printf("%e ", x[i]);
    printf("\n");

    return 0;
}
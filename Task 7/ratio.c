#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gaussrnd() 
{
  double u1 = drand48(); // random number
  double u2 = drand48();
  return -sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

int main(int argc, char** argv)
{
    double mx = atof(argv[1]);
    double sx = atof(argv[2]);
    double my = atof(argv[3]);
    double sy = atof(argv[4]);
    int N = atof(argv[5]);

    for(int i = 0; i<N; ++i)
    {
        double x = mx + sx*gaussrnd();
        double y = my + sy*gaussrnd();
        double ratio = x/y;
        printf("%e \n",ratio);
    }
    return 0;
}
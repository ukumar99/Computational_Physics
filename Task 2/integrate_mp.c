#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define g 9.81
#define pi 3.14159265358

double e(double y)
{
  return exp(y);
}


double w(double theta, double thetaNaught, double length)
{
  return sqrt(2*g*(cos(theta) - cos(thetaNaught))/length);
}

double dw(double Theta, double length, int N)
{
  double eps = -(2. * Theta)/(double)N;
  double res = 0;

  for (int i=0; i<N; i++)
  {
  res += -eps/w((Theta + i*eps) + (0.5*eps), Theta, length);
  printf("%lf\t \n", res);
  }
  return res;
}

//Midpoint
double mintegrate(double x1, double x2, int N)
{
  double eps = (x2-x1)/N;
  double res = 0;
  double lim = (double)N - (1/(double)N)*2;

  for(int i=0; i<N; i++){
  res += eps*e(x1+i*eps+(0.5*eps));
  }
  return res;
}

//Left hand 
double lintegrate(double x1, double x2, int N)
{
  double eps = (x2-x1)/N;
  double res = 0;
  for(int i=0; i<N; i++){
  res += eps*e(x1+i*eps);
  }
  return res;
}

int main(int argc, char *argv[])
{
  int N = atoi(argv[4]);

  double res = dw(atof(argv[2]) * pi/180., atof(argv[3]), N);

  printf("%.25lf %d\n", res, N);
  return 0;
}

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define pi 3.14159265
#define n 10000

int main() 
{
  double Ax;
  printf("Enter a value for Ax: ");
  scanf("%lf",&Ax);
 
  double Ay;
  printf("Enter a value for Ay: ");
  scanf("%lf",&Ay);
  
  double Phix;
  printf("Enter a value for Phix: ");
  scanf("%lf",&Phix);
  
  double Phiy; 
  printf("Enter a value for Phiy: ");
  scanf("%lf",&Phiy);

  double w_x;
  printf("Enter a value for w_x: ");
  scanf("%lf",&w_x);

  double w_y;
  printf("Enter a value for w_y: ");
  scanf("%lf",&w_y);

  FILE *pairs;
  pairs = fopen("pairs","w");
  double t,x_t,y_t;
 
  for (t=0; t<=2*pi; t+=(2.*pi/n))
  {
	x_t = Ax*sin(w_x*t+Phix);
	y_t = Ay*sin(w_y*t+Phiy);
	fprintf(pairs,"%lf %lf\n", x_t, y_t);	  
  }

  fclose(pairs);
  return 0;
}

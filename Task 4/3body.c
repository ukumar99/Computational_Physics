#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NBODY 3 
#define DIM 2

double mass[NBODY] = {1. ,3e-6 ,3.7e-8};

void compute_force(double* ri, double *rj, double mi, double mj, double* fij)
{
  double G=4*M_PI*M_PI;
  double rij[DIM];
  for(int a=0; a<DIM; ++a) rij[a] = rj[a]-ri[a];
  double rsq=0;
  for(int a=0; a<DIM; ++a) rsq += rij[a]*rij[a];
  double r = sqrt(rsq);
  for(int a=0; a<DIM; ++a) fij[a] = (G*mi*mj/(r*r*r))*rij[a];
}

void f(double *x, double t, double* res)
{
  for(int i=0; i<NBODY; ++i)
  {
    for(int a=0; a<DIM; ++a) res[i*2*DIM+a] = x[i*2*DIM+a+DIM];

    double totalforce[DIM];
    for(int a=0; a<DIM; ++a) totalforce[a] = 0;
   
    for(int j=0; j<NBODY; ++j) if(i!=j)
    {
      double force[DIM];
      compute_force(&x[i*2*DIM], &x[j*2*DIM], mass[i], mass[j], force);
      for(int a=0; a<DIM; ++a) totalforce[a] += force[a];
    }

    for(int a=0; a<DIM; ++a) res[i*2*DIM+a+DIM] = totalforce[a]/mass[i];
  }
}

void rk2(double* xn, double t, double dt, double* res)
{
  double k1[2*NBODY*DIM]; f(xn, t, k1);
  double xhalf[2*NBODY*DIM]; 
  for(int i=0; i<2*NBODY*DIM; ++i) xhalf[i] = xn[i] + k1[i]*dt/2;
  double k2[2*NBODY*DIM]; f(xhalf, t+dt/2, k2);

  for(int i=0; i<2*NBODY*DIM; ++i) res[i] = xn[i] + dt*k2[i];
}

int main(int argc, char** argv)
{
  double G=4*M_PI*M_PI;

  double dt = atof(argv[1]);
  double t = atof(argv[2]);
  double x1 = atof(argv[3]);
  double vx1 = atof(argv[4]);
  double y1 = atof(argv[5]);
  double vy1 = atof(argv[6]);
  double x2 = atof(argv[7]);
  double vx2 = atof(argv[8]);
  double y2 = atof(argv[9]);
  double vy2 = atof(argv[10]);
  double x3 = atof(argv[11]);
  double vx3 = atof(argv[12]);
  double y3 = atof(argv[13]);
  double vy3 = atof(argv[14]);

  vx1 = 1;
  vy1 = -1;
  vx2 = 0;
  vy2 = 6;
  vx3 = 0;
  vy3 = -6; 

  // vy2 = sqrt(G*(mass[0] + mass[1]) / 1);
  // vy3 = sqrt(G*(mass[0] + mass[1]) / 1) + sqrt(G*(mass[1] + mass[2]) / (x3-x2));
  
  int nstep = t/dt;
  
  double COM_v0_x = (mass[0]*vx1 + mass[1]*vx2 + mass[2]*vx3) / (mass[0]+mass[1]+mass[2]);  //x-center of mass
  double COM_v0_y = (mass[0]*vy1 + mass[1]*vy2 + mass[2]*vy3) / (mass[0]+mass[1]+mass[2]);  //y-center of mass
  
  vx1 -= COM_v0_x; 
  vy1 -= COM_v0_y;
  vx2 -= COM_v0_x;  
  vy2 -= COM_v0_y;
  vx3 -= COM_v0_x;  
  vy3 -= COM_v0_y;

  double x[4*NBODY] = {x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3};

  for(int i=0; i<nstep; ++i)
  {
    double xf[2*NBODY*DIM];
    rk2(x, i*dt, dt, xf);
    printf("%e ",i*dt+dt);
    
    for(int j=0; j<2*NBODY*DIM; ++j) printf("%e ", xf[j]);
    printf("\n");
    for(int j=0; j<2*NBODY*DIM; ++j) x[j] = xf[j];
    
  }

  return 0;
}
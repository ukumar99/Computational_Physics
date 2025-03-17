#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define GM 4.*M_PI*M_PI

void f(double* x, double* res)
{
    res[0] = x[2]; // x -> vx
    res[1] = x[3]; // y -> vy
    res[2] = GM/(x[0]*x[0] + x[1]*x[1])*-x[0]/sqrt(x[0]*x[0]+x[1]*x[1]); // vx -> ax
    res[3] = GM/(x[0]*x[0] + x[1]*x[1])*-x[1]/sqrt(x[0]*x[0]+x[1]*x[1]); // vy -> ay
}

void rk2(double* init, double dt, double* res) // Runge-Kutta integrator
{
    double k1[4];
    f(init, k1); // create 4-element vector

    double half[4]; 
    for(int i=0; i<4; i++)
        half[i] = init[i] + k1[i]*dt/2; // Runge-Kutta
    
    double k2[4];
    f(half, k2); // perform function f
    for(int i=0; i<4; ++i)
        res[i] = init[i] + dt*k2[i]; // steps
}

int main(int argc, char* argv[])
{
    double x0 = atof(argv[1]); // x
    double y0 = 0; // y
    double vx0 = 0; // xv
    double vy0 = 6.28; // yv
    double dt = atof(argv[2]); // dt
    double T = atof(argv[3]); // time

    double dts[7] = {0.0001,0.0003,0.001,0.003,0.01,0.03,0.1};
    
    double result[4];

        double state[4] = {x0,y0,vx0,vy0};
        
        for (double t = 0; t<T; t+=dt)
        {
        
        rk2(state, dt, result); // call RK and generate output
        for(int j=0; j<4; ++j) state[j] = result[j];
        
        double KE = 0.5*(state[2]*state[2]+state[3]*state[3]);
        double PE = -GM/sqrt(state[0]*state[0]+state[1]*state[1]);
        double SE = 0.5*(state[2]*state[2]+state[3]*state[3]) - GM/sqrt(state[0]*state[0]+state[1]*state[1]);
        
        if (t/5.-floor(t/5.)<0.01)
            // printf("%e %e %e %e %e\n", t+dt, result[0], result[1], result[2], result[3]);
            printf("%e %e %e %e %e\n",t,KE,PE,KE+PE,SE);
            // printf("%e %e\n", t+dt, OE);
        
        }
      
    return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include "linalgebra.h"

struct resistor
{
  int v1, v2;
  double R;
};

double compute_resistance(int nv, int nr, struct resistor* lr, int beg, int end)
{
  int neq = nr + nv - 1;

  double *m = (double*)malloc(neq*neq*sizeof(double));
  double *b = (double*)malloc(neq*sizeof(double));
  double *x = (double*)malloc(neq*sizeof(double));

  for(int i = 0; i < neq*neq; ++i) m[i] = 0;

  for(int i = 0; i < nr; ++i)
  {
    int v = lr[i].v1;
    if(v < nv - 1) m[(nr + v)*neq + i] = +1;
    v = lr[i].v2;
    if(v < nv - 1) m[(nr + v)*neq + i] = -1;
  }

  for(int i = 0; i < nv - 1; ++i) b[nr + i] = 0;
  if(beg < nv - 1) b[nr + beg] = +1;
  if(end < nv - 1) b[nr + end] = -1;

  for(int i = 0; i < nr; ++i)
  {
    m[i*neq + i] = -lr[i].R;
    int v = lr[i].v1;
    if(v < nv - 1) m[i*neq + nr + v] = +1;
    v = lr[i].v2;
    if(v < nv - 1) m[i*neq + nr + v] = -1;
    b[i] = 0;
  }

  printmat(m, neq);
  solve_system(m, neq, b, x);
  printmat(m, neq);
  double vbeg = (beg < nv - 1)? x[nr + beg] : 0;
  double vend = (end < nv - 1)? x[nr + end] : 0;

  free(b);
  free(x);
  free(m);

  return vbeg - vend;
}

double do_grid(int n, int in_i, int in_j, int out_i, int out_j)
{
  struct resistor* lr = (struct resistor*)malloc(2*n*n*sizeof(struct resistor));

  for(int i=0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      lr[i*n + j].R = 1;
      lr[i*n + j].v1 = i*n + j;
      lr[i*n + j].v2 = i*n + (j+1)%n;

      lr[i*n + j + n*n].R = 1;
      lr[i*n + j + n*n].v1 = n*i + j;
      lr[i*n + j + n*n].v2 = ((i + 1)%n)*n + j;
    }
  }

  double res = compute_resistance(n*n, 2*n*n, lr, n*in_i + in_j, n*out_i + out_j);
  free(lr);

  return res;
}

int main(int argc, char** argv)
{
  int n     = atoi(argv[1]);
  int in_i  = atoi(argv[2]);
  int in_j  = atoi(argv[3]);
  int out_i = atoi(argv[4]);
  int out_j = atoi(argv[5]);

  double res = do_grid(n, in_i, in_j, out_i, out_j);
  printf("Resistance: %f \n", res);

  return 0;
}

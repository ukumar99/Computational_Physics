#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void swap(double* a, double* b)
{
  double c = *a;
  *a = *b;
  *b = c;
}

void swap_rows(double* m, int N, int r1, int r2)
{
  for(int c = 0; c < N; ++c) swap(m + r1*N + c, m + r2*N + c);
}

int find_pivot(double* m, int N, int r)
{
  int p = r;
  double val = fabs(m[r*N + r]);
  for(int rp = r + 1; rp < N; ++rp)
    if(fabs(m[rp*N + r]) > val)
    {
      p = rp;
      val = fabs(m[rp*N + r]);
    }
  return p;
}

void elimination(double* m, int N, double* b)
{
  for(int r = 0; r < N - 1; ++r)
  {
    int p = find_pivot(m, N, r);
    if(p != r)
    {
      swap_rows(m, N, r, p);
      swap(b + r, b + p);
    }

    double piv = m[r*N + r];
    if(fabs(piv) < 1e-20)
    {
      printf("Matrix nearly singular. Exiting ... \n");
      abort();
    }

    for(int rp = r + 1; rp < N; ++rp)
    {
      double fact = -m[rp*N + r]/piv;
      for(int c = r; c < N; ++c) m[rp*N + c] += fact*m[r*N + c];
      b[rp] += fact*b[r];
    }
  }
}

void bks(double* m, int N, double* b, double* x)
{
  for(int r = N - 1; r > -1; --r)
  {
    double temp = b[r];
    for(int c = r + 1; c < N; ++c) temp -= m[r*N + c]*x[c];
    x[r] = temp/m[r*N + r];
  }
}

void solve_system(double* m, int N, double* b, double* x)
{
  elimination(m, N, b);
  bks(m, N, b, x);
}
	
void printmat(double* mat, int n)// prints a matrix "mat" with size n*n
{
	for(int r= 0; r<n; ++r)
	{
		for(int c = 0; c<n; ++c) printf("% 5.2f ", mat[r*n+c]);
		printf("\n");
	}
	printf("\n");
}


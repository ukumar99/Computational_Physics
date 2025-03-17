#include <stdio.h>

void swap(double* a, double* b);

void swap_rows(double* m, int N, int r1, int r2);

int find_pivot(double* m, int N, int r);

void elimination(double* m, int N, double* b);

void bks(double* m, int N, double* b, double* x);

void solve_system(double* m, int N, double* b, double* x);

void printmat(double* mat, int n);// prints a matrix "mat" with size n*n
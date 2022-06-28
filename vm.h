#ifndef VM_VM_H
#define VM_VM_H

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

double **dynamic_array_alloc( size_t N, size_t M );
void dynamic_array_free( double **A, size_t N );
double sum( const double *x, int n );
double *zeros( int n);
double *ones( int n );
double *random_array( int n, double min, double max );
double *random_real_array( int n );
double *range( double start, double stop );
double *arange( double start, double stop, double step );
double *array_sin( const double *x, int n );
double *array_cos( const double *x, int n );
double *array_exp( const double *x, int n );
double *array_ln( const double *x, int n );
double *array_abs( const double *x, int n );
double *array_power2( const double *x, int n );
double *array_power3( const double *x, int n );
double *array_power_n( const double *x, int n, int power );
double *array_sqrt( const double *x, int n );
double dot( const double *x, double *y, int n );
double *mxv( double **m,  double *v, int rows, int cols );
double **vxv( double *x, double *y, int n, int m );
double *plus( double *x, double *y, int n );
double *minus( double *x, double *y, int n );
double *multiply( double *x, double a, int n );
double **transpose( double **m, int rows, int columns );
int max_ind( const double *x, int n );
double max( const double *x, int n );
int min_ind( const double *x, int n );
double min( const double *x, int n );
int contains( const double *array, int x , int n );
double search( int a, const double *array, int n );
void private_sort( double *a, int l, int r );
void sort(double *a, int n );
double *reverse( double *x, int n );

#endif

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

double **dynamic_array_alloc( size_t N, size_t M )
{
    double **A = (double **)malloc(N*sizeof(double *));
    int i;
    for(i = 0; i < N; i++) {
        A[i] = (double *)malloc(M*sizeof(double ));
    }
    return A;
}

void dynamic_array_free( double **A, size_t N )
{
    int i;
    for(i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
}

double sum( const double *x, int n )
{
    double res = 0;
    int i;
    for( i = 0; i < n; i++ ){
        res += x[i];
    }
    return res;
}

double *zeros( int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for( i = 0; i < n; i++ ){
        res[i] = 0;
    }
    return res;

}

double *ones( int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for( i = 0; i < n; i++ ){
        res[i] = 1;
    }
    return res;
}

double *random_array( int n, double min, double max )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for( i = 0; i < n; i++ ){
        res[i] = (((double)rand())/RAND_MAX)*(max - min) + min;
    }
    return res;
}

double *random_real_array( int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for( i = 0; i < n; i++ ){
        res[i] = ((double)rand())/RAND_MAX;
    }
    return res;
}

double *range( double start, double stop )
{
    int n = (int)(start+stop) + 1;
    double *res = (double*)malloc(n*sizeof(double));
    int i = 0;
    while (start + i < stop)
    {
        res[i] = start + i;
        i++;
    }
    return res;
}

double *arange( double start, double stop, double step )
{
    int n = (int)((start+stop)/step) + 1;
    double *res = (double*)malloc(n*sizeof(double));
    int i = 0;
    while (start+i*step < stop)
    {
        res[i] = start+i*step;
        i++;
    }

    return res;
}

double *array_sin( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        res[i] = sin(x[i]);
    return res;
}

double *array_cos( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        res[i] = cos(x[i]);
    return res;
}

double *array_exp( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        res[i] = exp(x[i]);
    return res;
}

double *array_ln( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        if (x[i] > 0)
            res[i] = log(x[i]);
        else
            res[i] = NAN;

    return res;
}

double *array_abs( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        if (x[i] > 0)
            res[i] = x[i];
        else
            res[i] = -1*x[i];

    return res;
}

double *array_power2( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        res[i] = x[i]*x[i];
    return res;
}

double *array_power3( const double *x, int n)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        res[i] = x[i]*x[i]*x[i];
    return res;
}

double *array_power_n( const double *x, int n, int power)
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        if (x[i] == 0 && power == -1)
            res[i] = NAN;
        else
            res[i] = pow(x[i], power);
    return res;
}

double *array_sqrt( const double *x, int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i ++)
        if (x[i] < 0)
            res[i] = NAN;
        else
            res[i] = sqrt(x[i]);
    return res;
}

double dot( const double *x, double *y, int n )
{
    double res = 0;
    int i;

    for (i = 0; i < n; i++) {
        res += x[i] * y[i];
    }
    return res;

}

double *mxv( double **m,  double *v, int rows, int cols )
{
    double *res = (double*)malloc(rows*sizeof(double));
    int i;
    for (i = 0; i < rows; i++)
    {
        res[i] = dot(m[i], v, cols);
    }
}

double **vxv( double *x, double *y, int n, int m )
{
    double **res = dynamic_array_alloc(m, n);
    int i, j;
    for (i = 0; i < n; i++)
        for(j = 0; j < m; j++)
            res[i][j] = x[i]*y[j];
}

double *plus( double *x, double *y, int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]+y[i];
    }
    return res;
}

double *minus( double *x, double *y, int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]-y[i];
    }
    return res;
}

double *multiply( double *x, double a, int n )
{
    double *res = (double*)malloc(n*sizeof(double));
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]*a;
    }
    return res;
}
double **transpose( double **m, int rows, int columns )
{
    double **res = dynamic_array_alloc(columns, rows);
    int i, j;
    for(i = 0; i < rows; i++)
        for(j = 0; j < columns; j++)
            res[j][i] = m[i][j];

    return res;
}

int max_ind( const double *x, int n )
{
    int i, j = 0;
    double max = DBL_MIN;
    for(i = 0; i < n; i++) {
        if (x[i] > max)
            j = i;
            max = x[i];
    }
    return j;
}

double max( const double *x, int n )
{
    int i;
    double max = DBL_MIN;
    for(i = 0; i < n; i++) {
        if (x[i] > max)
        max = x[i];
    }
    return max;
}

int min_ind( const double *x, int n )
{
    int i, j = 0;
    double min = DBL_MAX;
    for(i = 0; i < n; i++) {
        if (x[i] < min)
            j = i;
        min = x[i];
    }
    return j;
}

double min( const double *x, int n )
{
    int i;
    double min = DBL_MAX;
    for(i = 0; i < n; i++) {
        if (x[i] < min)
            min = x[i];
    }
    return min;
}

int contains( const double *array, int x , int n )
{
    int i;
    for(i = 0; i < n; i++)
        if (array[i] == x)
            return 1;
    return 0;
}

double search( int a, const double *array, int n )
{
    int low, high, middle;
    low = 0;
    high = n - 1;
    while (low <= high)
    {
        middle = (low + high) / 2;
        if (a < array[middle])
            high = middle - 1;
        else if (a > array[middle])
            low = middle + 1;
        else
            return middle;
    }
    return -1;
}

void private_sort( double *a, int l, int r )
{
    if (l == r) return;
    int mid = (l + r) / 2;
    private_sort(a, l, mid);
    private_sort(a, mid + 1, r);
    int i = l;
    int j = mid + 1;
    double *tmp = (double *)malloc(r * sizeof(double ));
    int step;
    for (step = 0; step < r - l + 1; step++)
    {
        if ((j > r) || ((i <= mid) && (a[i] < a[j]))) {
            tmp[step] = a[i];
            i++;
        } else {
            tmp[step] = a[j];
            j++;
        }
    }

    for (step = 0; step < r - l + 1; step++)
        a[l + step] = tmp[step];
    //free(tmp);
}

void sort( double *x, int n )
{
    private_sort(x, 0, n-1);
}

double *reverse( double *x, int n )
{
    double *res = (double*)malloc(n * sizeof(double));
    int i;
    for(i = 0; i < n; i++)
        res[i] = x[n-i-1];
    return res;
}


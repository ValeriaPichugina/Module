#include <limits.h>
#include <float.h>
#include <stdlib.h>

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

double dot( double *x, double *y, int n )
{
    double res = 0;
    int i;

    for (i = 0; i < n; i++) {
        res += x[i] * y[i];
    }
    return res;

}

void mxv( double **m,  double *v, double *res, int rows, int cols )
{
    int i;
    for (i = 0; i < rows; i++)
    {
        res[i] = dot(m[i], v, cols);
    }
}

void vxv( double *x, double *y, double **res, int n, int m )
{
    int i, j;
    for (i = 0; i < n; i++)
        for(j = 0; j < m; j++)
            res[i][j] = x[i]*y[j];
}

void plus( double *x, double *y, double *res, int n )
{
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]+y[i];
    }
}

void minus( double *x, double *y, double *res, int n )
{
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]-y[i];
    }
}

void multiply( double *x, double a, double *res, int n )
{
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = x[i]*a;
    }
}
double **transpose( double **m, int rows, int columns )
{
    double **res = dynamic_array_alloc(columns, rows);
    int i, j;
    for(i = 0; i < rows; i++)
        for(j = 0; j < columns; j++)
            res[i][j] = m[j][i];

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

int max( const double *x, int n )
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

int min( const double *x, int n )
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

double binarysearch( int a, const double *array, int n )
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

void sort(double *a, int l, int r)
{
    if (l == r) return;
    int mid = (l + r) / 2;
    sort(a, l, mid);
    sort(a, mid + 1, r);
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
    free(tmp);
}

double *reverse( double *x, int n )
{
    double *res = (double *)malloc(n * sizeof(double));
    int i;
    for(i = 0; i < n; i++)
        res[i] = x[n-i-1];
    return res;
}


/* Geometric calculations */

#include <stdio.h>
#include <math.h>

double max(const double *a, const double *b)
{
    if (*a >= *b) {
        return *a;
    } else {
        return *b;
    }
}

double growth_rate(const double x[], const double y[], int *npt)
{
    /* Computes max panel growth rate for an airfoil */

    int i;
    double len1, len2, growth1, growth2, thismaxgrowth, maxgrowth;

    maxgrowth = 0.;
    len1 = sqrt(pow(x[1]-x[0], 2.) + pow(y[1]-y[0], 2.));
    for ( i = 1; i < *npt; i++ )
    {
        len2 = sqrt(pow(x[i+1]-x[1], 2.) + pow(y[i+1]-y[i], 2.));
        growth1 = len2/len1;
        thismaxgrowth = max(&growth1, &growth2);
        if (thismaxgrowth > maxgrowth) { 
            maxgrowth = thismaxgrowth;
        }
    }

    return maxgrowth;
}

double le_angle(const double xt[], const double yt[],
                const double xb[], const double yb[])
{
    /* Calculates leading edge panel angle */

    double angt, angb;

    angt = atan((yt[1]-yt[0])/(xt[1]-xt[0]));
    angb = atan((yb[0]-yb[1])/(xb[1]-xb[0]));

    return (angt-angb)*180./M_PI;
}

int main()
{
    double a = 2.;
    double b = 1.;
    printf("%f\n", max(&a, &b));
    return 0;
}

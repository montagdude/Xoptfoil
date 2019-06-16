/* Geometric calculations */

#include <vector>
#include <cmath>

double max(const double & a, const double & b)
{
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}

double growth_rate(const std::vector<double> & x, const std::vector<double> & y)
{
    /* Computes max panel growth rate for an airfoil */

    unsigned int i, npt;
    double len1, len2, growth1, growth2, thismaxgrowth, maxgrowth;

    npt = x.size();
    maxgrowth = 0.;
    len1 = sqrt(pow(x[1]-x[0], 2.) + pow(y[1]-y[0], 2.));
    for ( i = 1; i < npt; i++ )
    {
        len2 = sqrt(std::pow(x[i+1]-x[1], 2.) + std::pow(y[i+1]-y[i], 2.));
        growth1 = len2/len1;
        thismaxgrowth = max(growth1, growth2);
        if (thismaxgrowth > maxgrowth) { 
            maxgrowth = thismaxgrowth;
        }
    }

    return maxgrowth;
}

double le_angle(const std::vector<double> & xt, const std::vector<double> & yt,
                const std::vector<double> & xb, const std::vector<double> & yb)
{
    /* Calculates leading edge panel angle */

    double angt, angb;

    angt = atan((yt[1]-yt[0])/(xt[1]-xt[0]));
    angb = atan((yb[0]-yb[1])/(xb[1]-xb[0]));

    return (angt-angb)*180./M_PI;
}

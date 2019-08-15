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

    return (angt+angb)*180./M_PI;
}

int check_te_angle(const std::vector<double> & x,
                   const std::vector<double> & thickness, const double & x0,
                   const double & min_allowable_angle)
{
    /* Checks whether TE wedge angle behind x0 is less than allowed. Returns
       thickness array index of violation, or 0 if there are no violations. */

    unsigned int i, npt;
    double angle, dx;

    npt = x.size();
    for ( i = 0; i < npt-1; i++ )
    {
        if (x[i] < x0)
            continue;
        dx = x[npt-1] - x[i];
        angle = 2.*atan(0.5*thickness[i]/dx)*180./M_PI;
        if (angle < min_allowable_angle) {
            return i;
        }
    }

    return 0;
}

double derv1f(const double & u_plus2, const double & u_plus1, const double & u,
              const double & h)
{
    return (-3.*u + 4.*u_plus1 - u_plus2) / (2.*h);
}

double derv1b(const double & u_minus2, const double & u_minus1, const double & u,
              const double & h)
{
    return (3.*u - 4.*u_minus1 + u_minus2) / (2.*h);
}

double derv1c(const double & u_plus, const double & u_minus, const double & h)
{
    return (u_plus - u_minus) / (2.*h);
}

double derv2f(const double & u_plus2, const double & u_plus, const double & u,
              const double & h)
{
    return (u - 2.*u_plus + u_plus2) / std::pow(h,2.);
}

double derv2b(const double & u_minus2, const double & u_minus, const double & u,
              const double & h)
{
    return (u - 2.*u_minus + u_minus2) / std::pow(h,2.);
}

double derv2c(const double & u_plus, const double & u, const double & u_minus,
              const double & h)
{
    return (u_plus - 2.*u + u_minus) / std::pow(h,2.);
}

std::vector<double> curvature(const std::vector<double> & x, const std::vector<double> & y)
{
    /* Computes curvature for a function gam(s) = x(s) + y(s) */

    std::vector<double> curv, svec;
    unsigned int i, npt;
    double se, se2, xe, ye, xe2, ye2, xs, ys, xs2, ys2;

    npt = x.size();
    curv.resize(npt);
    svec.resize(npt);

    // Length vector s
    svec[0] = 0.;
    for ( i = 1; i < npt; i++ )
    {
        svec[i] = svec[i-1] + sqrt(std::pow(x[i]-x[i-1],2.) + std::pow(y[i]-y[i-1],2.));
    }

    // Compute first and second derivatives and curvature vector
    for ( i = 0; i < npt; i++ )
    {
        if (i == 0) {
            // Grid metric ds/de and d2s/de2
            se = derv1f(svec[i+2], svec[i+1], svec[i], 1.);
            se2 = derv2f(svec[i+2], svec[i+1], svec[i], 1.);

            // Derivatives of x and y with respect to grid parameter e
            xe = derv1f(x[i+2], x[i+1], x[i], 1.);
            ye = derv1f(y[i+2], y[i+1], y[i], 1.);
            xe2 = derv2f(x[i+2], x[i+1], x[i], 1.);
            ye2 = derv2f(y[i+2], y[i+1], y[i], 1.);

        } else if (i == npt-1) {
            // Grid metric ds/de and d2s/de2
            se = derv1b(svec[i-2], svec[i-1], svec[i], 1.);
            se2 = derv2b(svec[i-2], svec[i-1], svec[i], 1.);

            // Derivatives of x and y with respect to grid parameter e
            xe = derv1b(x[i-2], x[i-1], x[i], 1.);
            ye = derv1b(y[i-2], y[i-1], y[i], 1.);
            xe2 = derv2b(x[i-2], x[i-1], x[i], 1.);
            ye2 = derv2b(y[i-2], y[i-1], y[i], 1.);

        } else {
            // Grid metric ds/de and d2s/de2
            se = derv1c(svec[i+1], svec[i-1], 1.);
            se2 = derv2c(svec[i+1], svec[i], svec[i-1], 1.);

            // Derivatives of x and y with respect to grid parameter e
            xe = derv1c(x[i+1], x[i-1], 1.);
            ye = derv1c(y[i+1], y[i-1], 1.);
            xe2 = derv2c(x[i+1], x[i], x[i-1], 1.);
            ye2 = derv2c(y[i+1], y[i], y[i-1], 1.);
        }

        // Derivatives of x and y with respect to surface length s
        xs = xe/se;
        ys = ye/se;
        xs2 = (xe2 - se2/se*xe)/std::pow(se,2.);
        ys2 = (ye2 - ye2/se*ye)/std::pow(se,2.);

        // Curvature
        curv[i] = (xs*ys2 - ys*xs2) / std::pow(std::pow(xs,2.) + std::pow(ys,2.), 1.5);
    }

    return curv;
}

std::vector<double> curvature_reversals(const std::vector<double> & x,
                                        const std::vector<double> & curvature,
                                        const double & threshold)
{
    /* Locates curvature reversals in x */

    unsigned int i, npt;
    double curv1, curv2;
    std::vector<double> reversalsx;

    npt = x.size();
    reversalsx.resize(0);
    curv1 = 0.;
    for ( i = 1; i < npt-1; i++ )
    {
        if (std::abs(curvature[i]) >= threshold) {
            curv2 = curvature[i];
            if (curv2*curv1 < 0.) { reversalsx.push_back(x[i]); }
            curv1 = curv2;
        }
    }

    return reversalsx;
}

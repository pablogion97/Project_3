#include "Curvebase.h"
using namespace std;

double Curvebase::adaptiveSimpsons(double a, double b, double tol) {
    // Compute the integral of a function f using the Adaptive Simpsons Integration

    double i1 = simpsons(a, b);
    double i2 = simpsons(a, (a + b) / 2) + simpsons((a + b) / 2, b);

    if (abs(i1 - i2) < 15 * tol)
        return i2;

    return adaptiveSimpsons(a, (a + b) / 2, tol / 2) + adaptiveSimpsons((a + b)/2, b, tol/2);
}


double Curvebase::simpsons(double a, double b) {
    // Helper function to make the code more readable
    return (b - a) / 6 * (iFunc(a) + 4 * iFunc((a + b) / 2) + iFunc(b));
}


double Curvebase::iFunc(double p) {
    return sqrt(dxp(p) * dxp(p) + dyp(p) * dyp(p));
}


double Curvebase::x(const double& s) {
    if (s == 0.0)
        return xp(pmin);

    if (s == 1.0)
        return xp(pmax);

    return xp(newtons(s));
}


double Curvebase::y(const double& s) {
    if (s == 0.0)
        return yp(pmin);

    if (s == 1.0)
        return yp(pmax);

    return yp(newtons(s));
}


double Curvebase::newtons(const double& s) {
    double p = pmin * (1 - s) + pmax * s;
    double p_new;
    double err = 10, tol = 1e-10;
    int it = 0, maxIter = 50;

    while (it < maxIter && err > tol) {
        p_new = p - (integrate(p) - s * length) / iFunc(p);
        err = abs(p_new - p);
        p = p_new;
        it++;
    }

    if (it == maxIter) {
        cerr << "ERROR: Newtons method failed to converge." << endl;
        exit(EXIT_FAILURE);
    }

    return p;
}


double Curvebase::integrate(double p) {
    if (p<pmin || p>pmax) {
        cerr << "ERROR: P out of bounds." << endl;
        exit(EXIT_FAILURE);
    }
    else if (p == pmin)
        return 0.0;

	return adaptiveSimpsons(pmin, p, 1e-10);
}
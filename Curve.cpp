#include "Curve.h"

Curve::Curve(double a, double b, double offSet, bool dir, double peakPos)
	: peak(peakPos) {
	offset = offSet;
	rev = dir;
	pmin = (a > b) ? b : a;
	pmax = (a < b) ? b : a;
	length = integrate(pmax);
}


Curve::Curve(const Curve& other)
	: peak(other.peak) {
	offset = other.offset;
	length = other.length;
	rev = other.rev;
	pmin = other.pmin;
	pmax = other.pmax;
}


double Curve::xp(const double& p) {
	return p;
}


double Curve::yp(const double& p) {
	if (p < peak)
		return offset + 0.5 / (1.0 + exp(-3.0 * (p - peak + 3.0) ) );

	else
		return offset + 0.5 / (1.0 + exp(3.0 * (p - peak - 3.0) ) );
}


double Curve::dxp(const double& p) {
	return 1.0;
}


double Curve::dyp(const double& p) {

	double res;

	if (p < peak) {
		res = exp(-3.0 * (p - peak + 3.0));
		res =  1.5 * res / ((res + 1.0) * (res + 1.0));
	}
	else {
		res = exp(3.0 * (p - peak - 3.0));
		res = -1.5 * res / ((res + 1.0) * (res + 1.0));
	}

	return res;
}
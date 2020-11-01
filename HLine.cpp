#include "HLine.h"

HLine::HLine(double a, double b, double offSet, bool dir) {
	offset = offSet;
	length = abs(a - b);
	rev = dir;
	pmin = (a > b) ? b : a;
	pmax = (a < b) ? b : a;
}


HLine::HLine(const HLine& other) {
	offset = other.offset;
	length = other.pmax - other.pmin;
	rev = other.rev;
	pmin = other.pmin;
	pmax = other.pmax;
}


double HLine::xp(const double& p) {
	return p;
}


double HLine::yp(const double& p) {
	return offset;
}


double HLine::dxp(const double& p) {
	return 1.0;
}


double HLine::dyp(const double& p) {
	return 0.0;
}


double HLine::x(const double& s) {
	if (s > 1.0 || s < 0.0)
		exit(EXIT_FAILURE);

	return pmin*(1-s) + pmax*s;
}


double HLine::y(const double& s) {
	return offset;
}
#include "VLine.h"

VLine::VLine(double a, double b, double offSet, bool dir) {
	offset = offSet;
	length = abs(a - b);
	rev = dir;
	pmin = (a > b) ? b : a;
	pmax = (a < b) ? b : a;
}


VLine::VLine(const VLine& other) {
	offset = other.offset;
	length = other.pmax - other.pmin;
	rev = other.rev;
	pmin = other.pmin;
	pmax = other.pmax;
}


double VLine::xp(const double& p) {
	return offset;
}


double VLine::yp(const double& p) {
	return p;
}


double VLine::dxp(const double& p) {
	return 0.0;
}


double VLine::dyp(const double& p) {
	return 1.0;
}


double VLine::y(const double& s) {
	if (s > 1.0 || s < 0.0)
		exit(EXIT_FAILURE);

	return pmin * (1 - s) + pmax * s;
}


double VLine::x(const double& s) {
	return offset;
}
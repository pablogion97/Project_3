#pragma once
#ifndef VLINE_H
#define VLINE_H
#include "Curvebase.h"
class VLine :
	public Curvebase
{
public:
	double xp(const double& p) override;
	double yp(const double& p) override;
	double dxp(const double& p) override;
	double dyp(const double& p) override;
	double x(const double& s) override;
	double y(const double& s) override;

	VLine(double a, double b, double offSet, bool dir);
	VLine(const VLine& other);
};
#endif // !VLINE_H



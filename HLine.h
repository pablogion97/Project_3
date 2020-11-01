#pragma once
#ifndef HLINE_H
#define HLINE_H
#include "Curvebase.h"
class HLine :
    public Curvebase
{
public:
	double xp(const double& p) override;
	double yp(const double& p) override;
	double dxp(const double& p) override;
	double dyp(const double& p) override;
	double x(const double& s) override;
	double y(const double& s) override;

	HLine(double a, double b, double offSet, bool dir);
	HLine(const HLine& other);
};
#endif


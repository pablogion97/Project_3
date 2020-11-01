#pragma once
#include "Curvebase.h"

class Curve :
    public Curvebase
{
private:
	double peak;

public:
	double xp(const double& p) override;
	double yp(const double& p) override;
	double dxp(const double& p) override;
	double dyp(const double& p) override;

	Curve(double a, double b, double offSet, bool dir, double peakPos);
	Curve(const Curve& other);
};


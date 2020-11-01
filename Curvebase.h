#pragma once
#ifndef CURVEBASE_H
#define CURVEBASE_H
#include <cmath>
#include <iostream>
class Curvebase
{
private:
	double iFunc(double p);
	double simpsons(double a, double b);
	double adaptiveSimpsons(double a, double b, double tol);
	double newtons(const double& s);

protected:
	double pmin = 0.0; // minimal value for p
	double pmax = 0.0; // maximal value for p
	bool   rev = false;  // orientation of the curve
	double length = 0.0; // Curve length
	double offset = 0.0; // Curve offset

	// Curve parametrization
	virtual double xp(const double& p) = 0;
	virtual double yp(const double& p) = 0;

	// Curve derivatives
	virtual double dxp(const double& p) = 0;
	virtual double dyp(const double& p) = 0;

	// Integrate curve length
	double integrate(double p);

public:
	Curvebase() { ; };   //constructor
	virtual double x(const double& s);  //arc length parametrization
	virtual double y(const double& s);  //arc length parametrization
};

#endif // !CURVEBASE_H
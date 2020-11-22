#pragma once
#ifndef DOMAIN_H
#define DOMAIN_H
#include "Curve.h"
#include "HLine.h"
#include "VLine.h"
#include <fstream>

struct Point {
	double x, y;
};

class Domain {
private:
	int nx, ny;
	double* x_, * y_;

	Curvebase* sides[4];
	bool check_consistency() const;

public:
	Domain(Curvebase& bottom, Curvebase& right, Curvebase& top, Curvebase& left);
	Domain(const Domain& other);
	~Domain();

	Domain& operator=(const Domain& other);

	void generate_grid(const int& nX, const int& nY);
	void print_grid() const;
};

#endif // !DOMAIN_H
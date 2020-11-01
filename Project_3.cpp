// Project_3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "VLine.h"
#include "HLine.h"
#include "Curve.h"
#include "Domain.h"

using namespace std;

int main()
{
    // Generate curves
    Curve bottomCurve(-10, 5, 0, true, -3);
    VLine rightCurve(0, 3, 5, true);
    HLine topCurve(-10, 5, 3, false);
    VLine leftCurve(0, 3, -10, false);

    // Set grid spacings
    const int NX = 49;
    const int NY = 19;

    // Generate domain
    Domain domain(bottomCurve, rightCurve, topCurve, leftCurve);

    // Generate grid
    domain.generate_grid(NX, NY);

    // Export grid to csv file
    domain.print_grid();
}
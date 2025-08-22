#pragma once
#include <cassert>
#include <cmath>


double W(double r, double h) {
    double q = r / h;
    double sigma_3d = 8 / (M_PI * h * h * h);
    assert(q >= 0);
    if (q <= 0.5) {
        return sigma_3d * (6 * (q * q * q - q * q) + 1);
    } else if (q <= 1) {
        double tmp = 1 - q;
        return sigma_3d * (2 * tmp * tmp * tmp);
    } else {
        return 0;
    }
}

double gradW(double r, double h) {
    double q = r / h;
    double sigma_3d = 8 / (M_PI * h * h * h);
    assert(q >= 0);
    if (q <= 0.5) {
        return sigma_3d * 6 * q * (3 * q - 2) / (h * r);
    } else if (q <= 1) {
        double tmp = 1 - q;
        return -sigma_3d * 6 * tmp * tmp / (h * r);
    } else {
        return 0;
    }
}
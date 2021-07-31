#pragma once

#include <vector>
#include <cmath>

namespace helpers {
    double radians(double angle) {
        const double c = M_PI / 180;

        return angle * c;
    }

    std::vector<double> linspace(double a, double b, int n) {
        double h = (b - a) / (n - 1);
        std::vector<double> xs(n);

        std::vector<double>::iterator x;
        double val;

        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
            *x = val;
        
        return xs;
    }
}
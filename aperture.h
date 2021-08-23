#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <exception>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "config.h"
#include "helpers.h"
#include "ray.h"

typedef std::pair<std::vector<double>, std::vector<double>> point_2d;

enum class Instrument {Telescope, Autoguider, Finderscope};

class Aperture {
    private:
        float radius;
        float secRadius;
        float sampleRate;

        point_2d sampleDisk(double minRadius = 0);
        
        Eigen::Transform<double, 3, Eigen::Affine> transform(double hourAngle, double declination);
        Eigen::MatrixXd sampleAperture(double hourAngle, double declination, std::vector<double> x, std::vector<double> z);
        Eigen::Vector3d apertureDirection(double hourAngle, double declination);
        
        bool isRayBlocked(Eigen::Vector3d& origin, double hourAngle, double declination, double domeAzimuth);
        std::vector<bool> isBlocked(Eigen::MatrixXd& origins, double hourAngle, double declination, double domeAzimuth);

    public:
        Aperture(double radius, double secondaryRadius = 0, int sampleRate = 3);
        ~Aperture();

        double obstruction(double hourAngle, double declination, double domeAzimuth);
};
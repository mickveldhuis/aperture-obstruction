#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <exception>
#include <vector>
#include <cmath>
#include <tuple>

const double L1 = 0.814;
const double L2 = 1.098;
const double L3 = 0.439;
const double LAT = 53.24;

double radians(double angle) {
    const double c = M_PI / 180;

    return angle * c;
}

typedef std::pair<std::vector<double>, std::vector<double>> point_2d;

enum class Instrument {Telescope, Autoguider, Finderscope};

class Aperture {
    private:
        float radius;
        float secRadius;
        float sampleRate;

        point_2d sample_disk(double minRadius = 0.0);
        
        Eigen::Transform<double, 3, Eigen::Affine> transform(double hourAngle, double declination);
        Eigen::MatrixXd sampleAperture(double hourAngle, double declination, Eigen::VectorXd x, Eigen::VectorXd z);
        Eigen::Vector3d apertureDirection(double hourAngle, double declination);
        
        bool isRayBlocked(Eigen::Vector3d origin, double hourAngle, double declination, double domeAzimuth);
        bool isBlocked(Eigen::MatrixXd origins, double hourAngle, double declination, double domeAzimuth);

    public:
        Aperture(double radius, double secondaryRadius = 0.0, int sampleRate = 3);
        ~Aperture();

        double obstruction(double hourAngle, double declination, double domeAzimuth);
};
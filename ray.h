#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

struct capsule_t {
    double radius;
    double extent;
};

class Ray {
    private:
        Eigen::Vector3d origin;
        Eigen::Vector3d direction;
        
        capsule_t dome;
        
        double t;
        double delta = 1e-20;

    public:
        Ray(Eigen::Vector3d origin, Eigen::Vector3d direction);
        ~Ray();

        bool findIntersection();
        Eigen::Vector3d getIntersection();
};
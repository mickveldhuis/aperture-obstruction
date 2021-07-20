#include "ray.h"

Ray::Ray(Eigen::Vector3d origin, Eigen::Vector3d direction) {
    this->origin = origin;
    this-> direction = direction;

    const double domeDiameter = 6.25;
    
    capsule_t hemispericalDome;
    hemispericalDome.radius = domeDiameter/2;
    hemispericalDome.extent = 1.66;

    this->dome = hemispericalDome;
}

Ray::~Ray() {}

bool Ray::findIntersection() {
    bool hasIntersection = false;

    if (abs(this->direction(0)) < this->delta && abs(this->direction(1)) < this->delta) {
        double z = this->dome.extent + sqrt(this->dome.radius - pow(this->origin(0), 2) - pow(this->origin(1), 2));
        this->t = z - this->origin(2);

        hasIntersection = true;
        return hasIntersection;
    }
    
    double a0 = pow(this->origin(0), 2) + pow(this->origin(1), 2) - pow(this->dome.radius, 2);
    double a1 = this->origin(0)*this->direction(0) + this->origin(1)*this->direction(1);
    double a2 = pow(this->direction(0), 2) + pow(this->direction(1), 2);

    double discriminant = pow(a1, 2) - a0 * a2;

    if (discriminant >= 0.0) {
        this->t = ( -a1 + sqrt( discriminant ) ) / a2;
        hasIntersection = true;
    } else {
        return hasIntersection;
    }

    if (this->origin(2) + this->t * this->direction(2) >= this->dome.extent) {
        a0 = pow(this->origin(0), 2) + pow(this->origin(1), 2) - pow(origin(2) - this->dome.extent, 2) - pow(this->dome.radius, 2);
        a1 = this->origin(0)*this->direction(0) + this->origin(1)*this->direction(1) + (origin(2) - this->dome.extent) * this->direction(2);

        this->t = -a1 + sqrt(pow(a1, 2) - a0);
        hasIntersection = true;
    }

    return hasIntersection;
}

Eigen::Vector3d Ray::getIntersection() {

    if (this->t < this->delta) {
        throw std::underflow_error("The parameter 't' is nearly zero, i.e., 't' is not set..!");
    }

    return this->origin + this->t * this->direction;
}
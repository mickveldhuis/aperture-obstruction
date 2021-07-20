#include "aperture.h"
        
Aperture::Aperture(double radius, double secondaryRadius = 0.0, int sampleRate = 3) {
    this->radius = radius;
    this->secRadius = secondaryRadius;
    this->sampleRate = sampleRate;
}

Aperture::~Aperture() {
    // TODO: remove/de-init stuff....
}

point_2d Aperture::sample_disk(double minRadius = 0) {
    double dr = 1/this->sampleRate;

    std::vector<double> x;
    std::vector<double> y;

    // TODO: add sampling logic; figure out HOW TO 'linspace'

    return std::make_pair(x, y);
}

Eigen::Transform<double, 3, Eigen::Affine> Aperture::transform(double hourAngle, double declination) {
    // Creating the base to primary aperture transformation matrix 
    Eigen::Transform <double, 3, Eigen::Affine> H = Eigen::Transform <double, 3, Eigen::Affine>::Identity();

    // Frame 0 -> 1: pier to RA axis
    H.translate(Eigen::Vector3d(0.0, 0.0, L1));
    
    // Frame 1 -> 2: RA axis to Dec axis
    H.rotate(Eigen::AngleAxisd(radians(90-LAT), Eigen::Vector3d::UnitX()));
    H.rotate(Eigen::AngleAxisd(radians(-hourAngle), Eigen::Vector3d::UnitZ()));
    H.translate(Eigen::Vector3d(0.0, 0.0, L2));
    
    // Frame 2 -> 3: Dec axis to optical axis
    H.rotate(Eigen::AngleAxisd(radians(declination), Eigen::Vector3d::UnitX()));
    H.translate(Eigen::Vector3d(-L3, 0.0, 0.0));

    return H;
}

Eigen::MatrixXd Aperture::sampleAperture(double hourAngle, double declination, Eigen::VectorXd x, Eigen::VectorXd z) {
    
}

Eigen::Vector3d Aperture::apertureDirection(double hourAngle, double declination)  {
    
}

bool Aperture::isRayBlocked(Eigen::Vector3d origin, double hourAngle, double declination, double domeAzimuth) {
    
}

bool Aperture::isBlocked(Eigen::MatrixXd origins, double hourAngle, double declination, double domeAzimuth) {
    
}

double Aperture::obstruction(double hourAngle, double declination, double domeAzimuth) {

}
#include<iostream>
#include<getopt.h>

#include "aperture.h"
#include "config.h"

int main(int argc, char** argv) {
    double hourAngle;
    double declination;
    double domeAzimuth;
    int rate;

    // Use GetOpt to define command line arguments
    int opt;
    
    while ((opt = getopt(argc, argv, ":h:d:a:r:")) != -1) 
    {
        switch (opt) {
            case 'h':
                hourAngle = atof(optarg);
                break;
            case 'd':
                declination = atof(optarg);
                break;
            case 'a':
                domeAzimuth = atof(optarg);
                break;
            case 'r':
                rate = atoi(optarg);
                std::cout << rate << std::endl;
                break;
            case '?':
                std::cout << "unknown option..." << static_cast<char>(optopt) << std::endl;
                break;
            case ':':
                std::cout << "Missing argument for " << static_cast<char>(optopt) << std::endl;
                break;
        }
    }

    Aperture telescope(cfg::APERTURE_DIAMETER/2);
    
    double blockage = telescope.obstruction(hourAngle, declination, domeAzimuth);

    std::cout << "The obstruction is " << blockage << std::endl;

    return 0;
}
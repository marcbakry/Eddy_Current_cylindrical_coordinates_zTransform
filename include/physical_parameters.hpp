#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

struct PhysicalParameters
{
    double sigma = 1000.0; // conductivity
    double nu    = 1.0;    // reluctance

    PhysicalParameters(const double _sigma, const double _nu);
};

#endif
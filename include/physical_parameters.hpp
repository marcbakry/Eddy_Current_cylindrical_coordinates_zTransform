#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

struct PhysicalParameters
{
    double sigma = 1000.0;
    double nu    = 1.0;

    PhysicalParameters(const double _sigma, const double _nu);
};

#endif
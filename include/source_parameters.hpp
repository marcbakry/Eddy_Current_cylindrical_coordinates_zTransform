#ifndef SOURCE_PARAMETERS_HPP
#define SOURCE_PARAMETERS_HPP

#include <iostream>
#include <vector>
#include <algorithm>

// This class implements a constant source which can take multiple
// values over time, starting at t = 0. The inputs should be understood
// as "takes value j[i] between _t[i] and _t[i+1]". The final time 
// should be excluded as it will be added automatically. A size check
// is performed at construction. The first time value must be lower or
// equal to 0 or the constructor will fail.
class SourceParameters {
public:
    SourceParameters(std::vector<double> const& _j, std::vector<double> const& _t);

    // compute the source for a given time _t
    double compute(double _t);

private:
    std::vector<double> m_j;
    std::vector<double> m_t;
};

#endif
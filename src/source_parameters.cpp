#include "source_parameters.hpp"

SourceParameters::SourceParameters(std::vector<double> const &_j, std::vector<double> const& _t) {
    if(_j.size() != _t.size()) {
        std::cout << "ERROR: SourceParameters::SourceParameters(): source and time sizes do not match. The program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if(!std::is_sorted(_t.begin(),_t.end())) {
        std::cout << "ERROR: SourceParameters::SourceParameters(): the time vector is not sorted. Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if(_t.front() > 0.0) {
        std::cout << "ERROR: SourceParameters::SourceParameters(): the first time value is greater strictly than 0.0. The program will exit." << std:endl;
        std::exit(EXIT_FAILURE);
    }
    m_j = _j; m_t = _t; m_t.push_back(1E50); // add 'impossible' value
}

double SourceParameters::compute(double _t) {
    for(const auto i=0; i<m_t.size()-1; ++i) {
        if(_t >= m_t[i] && _t < m_t[i+1]) return m_j[i];
    }
    return 0.0;
}
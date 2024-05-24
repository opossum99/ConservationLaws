#ifndef CONSERVATIONLAWSPROJECT_ANALYTICALSOLUTION_HPP
#define CONSERVATIONLAWSPROJECT_ANALYTICALSOLUTION_HPP

#include <vector>
#include "Aux.hpp"

std::vector<RUVP> analytical_test(const std::vector<std::pair<double, double>> &xy, double x_s, double time, const RUVP &parL,
                     const RUVP &parR);

#endif //CONSERVATIONLAWSPROJECT_ANALYTICALSOLUTION_HPP

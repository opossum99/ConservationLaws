#ifndef CONSERVATIONLAWSPROJECT_INITIAL_HPP
#define CONSERVATIONLAWSPROJECT_INITIAL_HPP

#include <vector>
#include "Solver.hpp"

void x_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy);

void y_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy);

void cyl_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy);

void diag_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy);

void und(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy);

#endif //CONSERVATIONLAWSPROJECT_INITIAL_HPP

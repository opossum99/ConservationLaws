#include <ceres/ceres.h>

#include "AnalyticalSolution.hpp"

std::vector<RUVP>
analytical_test(const std::vector<std::pair<double, double>> &xy, double x_s, double time, const RUVP &parL,
                const RUVP &parR) {
    std::vector<RUVP> anal_solution(xy.size());
    double p_star = 0.30313;
    double vel_star = 0.92745;
    double r_L_star = 0.42632;
    double r_R_star = 0.26557;
    double c_L = sqrt(GAMMA * parL.p / parL.r);
    double c_L_star = sqrt(GAMMA * p_star / r_L_star);
    double D = (p_star - parR.p) / (parR.r * vel_star);
    for (int i = 0; i < xy.size(); i++) {
        auto coor = xy[i].first - x_s;
        if (coor < -c_L * time) {
            anal_solution[i].p = parL.p;
            anal_solution[i].u = 0;
            anal_solution[i].v = 0;
            anal_solution[i].r = parL.r;
        } else if (coor < (vel_star - c_L_star) * time) {
            anal_solution[i].u = vel_star * (coor + c_L * time) / (vel_star - c_L_star + c_L) / time;
            anal_solution[i].p =
                    parL.p * std::pow(1 - (GAMMA - 1) / 2 * (anal_solution[i].u) / c_L, 2 * GAMMA / (GAMMA - 1));
            anal_solution[i].v = 0;
            anal_solution[i].r = parL.r * std::pow(1 - (GAMMA - 1) / 2 * (anal_solution[i].u) / c_L, 2 / (GAMMA - 1));
        } else if (coor < vel_star * time) {
            anal_solution[i].p = p_star;
            anal_solution[i].u = vel_star;
            anal_solution[i].v = 0;
            anal_solution[i].r = r_L_star;
        } else if (coor < D * time) {
            anal_solution[i].p = p_star;
            anal_solution[i].u = vel_star;
            anal_solution[i].v = 0;
            anal_solution[i].r = r_R_star;
        } else {
            anal_solution[i].p = parR.p;
            anal_solution[i].u = parR.u;
            anal_solution[i].v = parR.v;
            anal_solution[i].r = parR.r;
        }
    }
    return anal_solution;
}

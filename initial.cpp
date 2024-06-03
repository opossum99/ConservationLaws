#include "initial.hpp"
#include "Solver.hpp"


void x_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy) {
    for (int i = 0; i < xy.size(); i++) {
        if (xy[i].first < 0.5) {
//            initial[i].r = 5.99924;
//            initial[i].p = 460.894;
//            initial[i].u = 19.5975;
//            initial[i].v = 0;
//            initial[i].r = 1.;
//            initial[i].p = .4;
//            initial[i].u = -2.;
//            initial[i].v = 0;
            initial[i].r = 1.;
            initial[i].p = 1.;
            initial[i].u = 0.;
            initial[i].v = 0;
        } else {
//            initial[i].r = 1;
//            initial[i].p = 0.4;
//            initial[i].u = 2.;
//            initial[i].v = 0;
            initial[i].r = .125;
            initial[i].p = 0.1;
            initial[i].u = 0.;
            initial[i].v = 0;
//            initial[i].r = 5.99242;
//            initial[i].p = 46.0950;
//            initial[i].u = -6.19633;
//            initial[i].v = 0;
        }
    }
}

void y_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy) {
    for (int i = 0; i < xy.size(); i++) {
        if (xy[i].second < 0.5) {
            initial[i].r = 5.99924;
            initial[i].p = 460.894;
            initial[i].v = 19.5975;
            initial[i].u = 0;
//            initial[i].r = 1.;
//            initial[i].p = .4;
//            initial[i].v = -2.;
//            initial[i].u = 0;
//            initial[i].r = 1.;
//            initial[i].p = 1.;
//            initial[i].v = 0.;
//            initial[i].u = 0;
        } else {
//            initial[i].r = 1;
//            initial[i].p = 0.4;
//            initial[i].v = 2.;
//            initial[i].u = 0;
//            initial[i].r = .125;
//            initial[i].p = 0.1;
//            initial[i].v = 0.;
//            initial[i].u = 0;
            initial[i].r = 5.99242;
            initial[i].p = 46.0950;
            initial[i].v = -6.19633;
            initial[i].u = 0;
        }
    }
}

void cyl_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy) {
    const std::pair<double, double> center = {1., 1.};
    const double radius = 0.3;
    for (int i = 0; i < xy.size(); i++) {
        initial[i].u = 0;
        initial[i].v = 0;
        if ((center.second - xy[i].second) * (center.second - xy[i].second) +
            (center.first - xy[i].first) * (center.first - xy[i].first) < radius) {
            initial[i].r = 1.;
            initial[i].p = 1.;
        } else {
            initial[i].r = .125;
            initial[i].p = .1;
        }
    }
}

void diag_Riemann(std::vector<RUVP> &initial, const std::vector<std::pair<double, double>> &xy) {
    for (int i = 0; i < xy.size(); i++) {
        initial[i].u = 0;
        initial[i].v = 0;
        if (xy[i].first + xy[i].second < .5) {
            initial[i].r = 1.;
            initial[i].p = 1.;
        } else {
            initial[i].r = .125;
            initial[i].p = .1;
        }
    }
}


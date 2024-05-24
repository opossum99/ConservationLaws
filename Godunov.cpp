#include "Godunov.hpp"
#include "Aux.hpp"
#include <ceres/ceres.h>

Godunov::Godunov(const RUP &left, const RUP &right) {
    rup_['l'] = left;
    rup_['r'] = right;

    par_['l'].a = GAMMA * left.p / left.r;
    par_['l'].A = 2 / (GAMMA + 1) / left.r;
    par_['l'].B = (GAMMA - 1) / (GAMMA + 1) * left.p;
    par_['r'].a = GAMMA * right.p / right.r;
    par_['r'].A = 2 / (GAMMA + 1) / right.r;
    par_['r'].B = (GAMMA - 1) / (GAMMA + 1) * right.p;
};

double Godunov::flux(char l_or_r, double p) {
    if (p > rup_[l_or_r].p) {
        return (p - rup_[l_or_r].p) * std::sqrt(par_[l_or_r].A / (p + par_[l_or_r].B));
    } else {
        return 2 * par_[l_or_r].a / (GAMMA - 1) * (std::pow(p / rup_[l_or_r].p, (GAMMA - 1) / 2 / GAMMA) - 1);
    }
}

std::pair<double, double> Godunov::find_UP() {
    auto fun = [&](double p) {
        return (flux('r', p) + flux('l', p) + rup_['r'].u - rup_['l'].u) *
               (flux('r', p) + flux('l', p) + rup_['r'].u - rup_['l'].u);
    };

    return std::pair<double, double>();
}

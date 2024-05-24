#ifndef CONSERVATIONLAWSPROJECT_GODUNOV_HPP
#define CONSERVATIONLAWSPROJECT_GODUNOV_HPP

#include <functional>
#include <map>
#include "Aux.hpp"

struct RUP{
    double r;
    double u;
    double p;
};

struct aAB{
    double a;
    double A;
    double B;
};


class Godunov {
public:
    Godunov(const RUP& left, const RUP& right);

    double flux(char l_or_r, double p);

    std::pair<double, double> find_UP();

    template <typename T>
    bool operator()(const T* const p, T* residual) {
        residual[0] = flux('r', p) + flux('l', p) + rup_['r'].u - rup_['l'].u;
        return true;
    }

private:
    std::map<char, RUP> rup_;
    std::map<char, aAB> par_;
};


#endif //CONSERVATIONLAWSPROJECT_GODUNOV_HPP

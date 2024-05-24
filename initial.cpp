#include "initial.hpp"

void x_Riemann(std::vector<RUVP>& initial, const std::vector<std::pair<double, double>>& xy){
    for(int i = 0; i < xy.size(); i++){
        initial[i].u = 0;
        initial[i].v = 0;
        if(xy[i].first < 0.5){
            initial[i].r = 1.;
            initial[i].p = 1.;
        } else {
            initial[i].r = .125;
            initial[i].p = .1;
        }
    }
}

void y_Riemann(std::vector<RUVP>& initial, const std::vector<std::pair<double, double>>& xy){
    for(int i = 0; i < xy.size(); i++){
        initial[i].u = 0;
        initial[i].v = 0;
        if(xy[i].second < 0.5){
            initial[i].r = 1.;
            initial[i].p = 1.;
        } else {
            initial[i].r = .125;
            initial[i].p = .1;
        }
    }
}

void cyl_Riemann(std::vector<RUVP>& initial, const std::vector<std::pair<double, double>>& xy){
    for(int i = 0; i < xy.size(); i++){
        initial[i].u = 0;
        initial[i].v = 0;
        if((1 - xy[i].second)*(1 - xy[i].second) + (1 - xy[i].first)*(1 - xy[i].first) < 0.3){
            initial[i].r = 1.;
            initial[i].p = 1.;
        } else {
            initial[i].r = .125;
            initial[i].p = .1;
        }
    }
}


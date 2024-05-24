#ifndef CONSERVATIONLAWSPROJECT_SOLVER_HPP
#define CONSERVATIONLAWSPROJECT_SOLVER_HPP

#include <vector>
#include <Eigen/Dense>
#include <string>
#include "Geometry.hpp"
#include "Aux.hpp"

auto Q2var(Eigen::Vector4d Q);

auto Q2F(Eigen::Vector4d Q);

auto Q2G(Eigen::Vector4d Q);

auto var2Q(RUVP var);


class Solver {
public:
    Solver(const std::vector<std::pair<double, double>> &mesh, const std::vector<RUVP> &init_solution, double dt,
           int N_x, double courant = 1.);

    void LaxFriedrichs(double T_end);

    void Godunov();

    void MUSCL_PLM(double T_end);

    void boundaries();

    void CFL();

    template<class T>
    auto ij2k(T i, T j){
        return j * stride_ + i;
    }

    template<class T>
    std::pair<T, T> k2ij(T k){
        std::pair<T, T> ij;
        ij.second = k / stride_;
        ij.first = k % stride_;
        return ij;
    }

    double compare(const std::vector<RUVP> &expected);
    void VisualVTK(const std::string& filename);

private:
    std::vector<std::pair<double, double>> xy_;
    std::vector<RUVP> solution_;
    std::vector<Eigen::Vector4d> Q_;
    std::vector<Eigen::Vector4d> F_;
    std::vector<Eigen::Vector4d> G_;

    Geometry geometry_;
    double x_step_{};
    double y_step_{};
    double dt_;
    int stride_; // Number of cells in x-direction and stride in arrays.
    double time_;
    double cour_;
};


#endif //CONSERVATIONLAWSPROJECT_SOLVER_HPP

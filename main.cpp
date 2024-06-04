#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <chrono>

#include "Solver.hpp"
#include "initial.hpp"


int main() {
    std::cout << "Hello, World!" << std::endl;
    Eigen::Matrix2d a;
    a << 0, 2, -2, 0;
    std::cout << a << std::endl;
    const double length = 2.;
    const double height = 2.;
    const double fraction = 1.;
    const int N_y_inner = 1000;
    const int N_x_inner = N_y_inner * static_cast<int>(length / height * fraction);
    const int N_y = N_y_inner + 2;
    const int N_x = N_x_inner + 2;
    double x_step = length / N_x_inner;
    double y_step = height / N_y_inner;

    std::vector<std::pair<double, double>> init_mesh(N_x * N_y);
    for (int i = 0; i < init_mesh.size(); i++) {
        init_mesh[i].first = x_step * (i % N_x) - x_step / 2;
        init_mesh[i].second = y_step * static_cast<double>(i / N_x) - y_step / 2;
    }


    std::vector<RUVP> initial_conditions(N_x * N_y);
    cyl_Riemann(initial_conditions, init_mesh);
    const double T_end = 0.75;
    double dt = 0.01;
    Solver lax_solver(init_mesh, initial_conditions, dt, N_x, 0.5);
    lax_solver.LaxFriedrichs(T_end);
    lax_solver.VisualVTK("../data.vtk");

//    Solver mhm_solver(init_mesh, initial_conditions, dt, N_x, 0.5);
//    mhm_solver.MUSCL_MHM(T_end);
//    mhm_solver.VisualVTK("mhm.vtk");

//    auto anal = analytical_test_x(init_mesh, 0.5, T_end, {5.99924, 19.5975, 0, 460.894}, {5.99924, -6.19633, 0, 46.0950});
    //auto anal = analytical_test_x(init_mesh, 0.5, T_end, {1., 0, 0, 1.}, {.125, 0, 0, .1});
//    Solver anal_sol(init_mesh, anal, dt, N_x, 0.999);
//    anal_sol.VisualVTK("anal.vtk");
//    std::ofstream residual("res_mhm.txt", std::ios::app);
//    residual << x_step << "  " << mhm_solver.compare(anal) << std::endl;
//    residual.close();
    return 0;
}

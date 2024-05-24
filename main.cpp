#include <iostream>
#include <Eigen/Dense>
#include <fstream>

#include "Solver.hpp"
#include "initial.hpp"
#include "AnalyticalSolution.hpp"

int main() {
    std::cout << "Hello, World!" << std::endl;
    Eigen::Matrix2d a;
    a << 0, 2, -2, 0;
    std::cout << a << std::endl;
    const double length = 1.;
    const double height = 1.;
    const double fraction = 1.;
    const int N_y_inner = 500;
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
    x_Riemann(initial_conditions, init_mesh);

    const double T_end = 0.25;
    double dt = 0.01;
    Solver solver(init_mesh, initial_conditions, dt, N_x, 0.999);
    solver.LaxFriedrichs(T_end);
    solver.VisualVTK("data.vtk");

    auto anal = analytical_test(init_mesh, 0.5, T_end, {1., 0, 0, 1.}, {0.125, 0, 0, 0.1});
    Solver anal_sol(init_mesh, anal, dt, N_x, 0.999);
    anal_sol.VisualVTK("anal.vtk");
    std::ofstream residual("res.txt", std::ios::app);
    residual << x_step << "  " << solver.compare(anal) << std::endl;
    residual.close();
    return 0;
}

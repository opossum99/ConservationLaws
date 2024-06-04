#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>


__device__ void Q2F(double* Q, double* F, int idx, double gamma) {
    double rho = Q[idx * 4];
    double u = Q[idx * 4 + 1] / rho;
    double v = Q[idx * 4 + 2] / rho;
    double E = Q[idx * 4 + 3];
    double p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = u * (E + p);
}

__device__ void Q2G(double* Q, double* G, int idx, double gamma) {
    double rho = Q[idx * 4];
    double u = Q[idx * 4 + 1] / rho;
    double v = Q[idx * 4 + 2] / rho;
    double E = Q[idx * 4 + 3];
    double p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

    G[0] = rho * v;
    G[1] = rho * u * v;
    G[2] = rho * v * v + p;
    G[3] = v * (E + p);
}

// Kernel to update the solution
__global__ void update_solution(double* Q, double* Q_new, int nx, int ny, double dx, double dy, double dt, double gamma) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
        int idx = i + j * nx;
        int idx_ip1 = (i + 1) + j * nx;
        int idx_im1 = (i - 1) + j * nx;
        int idx_jp1 = i + (j + 1) * nx;
        int idx_jm1 = i + (j - 1) * nx;

        double F_ip1[4], F_im1[4], G[4], G_jp1[4], G_jm1[4];
        Q2F(Q, F_ip1, idx_ip1, gamma);
        Q2F(Q, F_im1, idx_im1, gamma);
        Q2G(Q, G_jp1, idx_jp1, gamma);
        Q2G(Q, G_jm1, idx_jm1, gamma);

        // Lax-Friedrichs update
        for (int k = 0; k < 4; k++) {
            Q_new[idx * 4 + k] = 0.25 * (Q[idx_ip1 * 4 + k] + Q[idx_im1 * 4 + k] + Q[idx_jp1 * 4 + k] + Q[idx_jm1 * 4 + k])
                - dt / (2 * dx) * (F_ip1[k] - F_im1[k])
                - dt / (2 * dy) * (G_jp1[k] - G_jm1[k]);
        }
    }
}

// Kernel to apply boundary conditions
__global__ void apply_boundary_conditions(double* Q, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Apply boundary conditions at the left and right boundaries
    if (i == 0) {
        for (int k = 0; k < ny; k++) {
            int idx = i + k * nx;
            Q[idx * 4 + 1] = 0.0; // Set horizontal velocity (rho*u) to 0
            // Other are equal
            Q[idx * 4 + 0] = Q[(i + 1 + k * nx) * 4 + 0];
            Q[idx * 4 + 2] = Q[(i + 1 + k * nx) * 4 + 2];
            Q[idx * 4 + 3] = Q[(i + 1 + k * nx) * 4 + 3];
        }
    }

    // Apply boundary conditions at the top and bottom boundaries
    if (j == 0) {
        for (int k = 0; k < nx; k++) {
            int idx = k + j * nx;
            Q[idx * 4 + 2] = 0.0; // Set vertical velocity (rho*v) to 0
            // Other are equal
            Q[idx * 4 + 0] = Q[(k + (j + 1) * nx) * 4 + 0];
            Q[idx * 4 + 1] = Q[(k + (j + 1) * nx) * 4 + 1];
            Q[idx * 4 + 3] = Q[(k + (j + 1) * nx) * 4 + 3];
        }
    }

    // Apply boundary conditions at the left and right boundaries
    if (i == nx - 1) {
        for (int k = 0; k < ny; k++) {
            int idx = i + k * nx;
            Q[idx * 4 + 1] = 0.0; // Set horizontal velocity (rho*u) to 0
            // Other are equal
            Q[idx * 4 + 0] = Q[(i - 1 + k * nx) * 4 + 0];
            Q[idx * 4 + 2] = Q[(i - 1 + k * nx) * 4 + 2];
            Q[idx * 4 + 3] = Q[(i - 1 + k * nx) * 4 + 3];
        }
    }

    // Apply boundary conditions at the top and bottom boundaries
    if (j == ny - 1) {
        for (int k = 0; k < nx; k++) {
            int idx = k + j * nx;
            Q[idx * 4 + 2] = 0.0; // Set vertical velocity (rho*v) to 0
            // Other are equal
            Q[idx * 4 + 0] = Q[(k + (j - 1) * nx) * 4 + 0];
            Q[idx * 4 + 1] = Q[(k + (j - 1) * nx) * 4 + 1];
            Q[idx * 4 + 3] = Q[(k + (j - 1) * nx) * 4 + 3];
        }
    }
}

// Kernel to compute the maximum wave speed
__global__ void compute_max_wave_speed(double* Q, double* max_wave_speed, int nx, int ny, double gamma) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nx && j < ny) {
        int idx = i + j * nx;

        double rho = Q[idx * 4 + 0];
        double u = Q[idx * 4 + 1] / rho;
        double v = Q[idx * 4 + 2] / rho;
        double E = Q[idx * 4 + 3];
        double p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
        double c = sqrt(gamma * p / rho);

        max_wave_speed[idx] = fabs(u) + c + fabs(v) + c;
    }
}


// Function to write data to a VTK file
void write_vtk(const std::vector<double>& Q, int nx, int ny, int timestep, double dx, double dy, double gamma) {
    std::ostringstream filename;
    filename << "output_" << std::setw(4) << std::setfill('0') << timestep << ".vtk";
    std::ofstream file(filename.str());
    std::vector<double> density(Q.size()/4);
    std::vector<double> velocity_x(Q.size() / 4);
    std::vector<double> velocity_y(Q.size() / 4);
    std::vector<double> pressure(Q.size() / 4);

    for (int idx = 0; idx < Q.size() / 4; idx++) {
        density[idx] = Q[idx * 4];
        velocity_x[idx] = Q[idx * 4 + 1]/density[idx];
        velocity_y[idx] = Q[idx * 4 + 2]/density[idx];
        pressure[idx] = (Q[idx * 4 + 3] - 0.5 * density[idx] * (velocity_x[idx] * velocity_x[idx] + velocity_y[idx] * velocity_y[idx])) * (gamma - 1.);
    }


    file << "# vtk DataFile Version 2.0\n";
    file << "Euler Equations\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx * ny << " double\n";

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << i * dx << " " << j * dy << " 0\n";
        }
    }

    file << "POINT_DATA " << nx * ny << "\n";

    file << "SCALARS density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = i + j * nx;
            file << std::fixed << density[idx + 0] << "\n";
        }
    }

    file << "SCALARS velocity_x double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = i + j * nx;
            file << std::fixed << velocity_x[idx] << "\n";
        }
    }

    file << "SCALARS velocity_y double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = i + j * nx;
            file << std::fixed << velocity_y[idx] << "\n";
        }
    }

    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = i + j * nx;
            file << std::fixed << pressure[idx] << "\n";
        }
    }

    file.close();
}

void cyl_Riemann(std::vector<double>& initial, int nx, int ny, double dx, double dy, double gamma) {
    const std::pair<double, double> center = { 1., 1. };
    const double radius = 0.3;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int idx = i + j * nx;
            if ((i * dx - center.first) * (i * dx - center.first) + (j * dy - center.second) * (j * dy - center.second) < radius) {
                initial[idx * 4 + 0] = 1.0;  // rho
                initial[idx * 4 + 1] = 0.0;  // rho*u
                initial[idx * 4 + 2] = 0.0;  // rho*v
                initial[idx * 4 + 3] = 1.0 / (gamma - 1.);  // E
            }
            else {
                initial[idx * 4 + 0] = 0.125;  // rho
                initial[idx * 4 + 1] = 0.0;  // rho*u
                initial[idx * 4 + 2] = 0.0;  // rho*v
                initial[idx * 4 + 3] = 0.1 / (gamma - 1.);  // E
            }
        }
    }
}

int main() {
    // Constants
    const int NX = 2001;  // Number of grid points in x direction
    const int NY = 2001;  // Number of grid points in y direction
    const double gamma = 1.4;  // Ratio of specific heats
    const double dx = 2.0 / (NX - 1);  // Grid spacing in x direction
    const double dy = 2.0 / (NY - 1);  // Grid spacing in y direction
    const double cfl_number = 0.75;  // CFL number

    // Allocate host memory
    std::vector<double> Q(NX * NY * 4);
    std::vector<double> Q_new(NX * NY * 4);
    std::vector<double> max_wave_speed(NX * NY);

    // Initialize Q (this would typically involve setting initial conditions)
    //for (int j = 0; j < NY; j++) {
    //    for (int i = 0; i < NX; i++) {
    //        int idx = i + j * NX;
    //        if (i * dx < 1.) {
    //            Q[idx * 4 + 0] = 1.0;  // rho
    //            Q[idx * 4 + 1] = 0.0;  // rho*u
    //            Q[idx * 4 + 2] = 0.0;  // rho*v
    //            Q[idx * 4 + 3] = 1.0 / (gamma - 1.);  // E
    //        }
    //        else {
    //            Q[idx * 4 + 0] = 0.125;  // rho
    //            Q[idx * 4 + 1] = 0.0;  // rho*u
    //            Q[idx * 4 + 2] = 0.0;  // rho*v
    //            Q[idx * 4 + 3] = 0.1 / (gamma - 1.);  // E
    //        }        
    //    }
    //}

    cyl_Riemann(Q, NX, NY, dx, dy, gamma);



    // Allocate device memory
    double* d_Q, * d_Q_new, * d_max_wave_speed;
    cudaMalloc(&d_Q, NX * NY * 4 * sizeof(double));
    cudaMalloc(&d_Q_new, NX * NY * 4 * sizeof(double));
    cudaMalloc(&d_max_wave_speed, NX * NY * sizeof(double));

    // Copy data from host to device
    cudaMemcpy(d_Q, Q.data(), NX * NY * 4 * sizeof(double), cudaMemcpyHostToDevice);

    // Define grid and block dimensions
    dim3 blockDim(16, 16);
    dim3 gridDim((NX + blockDim.x - 1) / blockDim.x, (NY + blockDim.y - 1) / blockDim.y);

    double Time = 0.0;       // Time
    double Time_end = 0.75;  // End
    int n = 0;               // Counter

    double another_time = 0.0;

    auto begin = std::chrono::steady_clock::now();

    while (Time < Time_end) {
        // Increase counter
        //if (Time >= another_time) {
        //    write_vtk(Q, NX, NY, n++, dx, dy, gamma);
        //    another_time += 0.01;
        //}

        // CFL Control
        compute_max_wave_speed <<<gridDim, blockDim >>> (d_Q, d_max_wave_speed, NX, NY, gamma);
        cudaMemcpy(max_wave_speed.data(), d_max_wave_speed, NX * NY * sizeof(double), cudaMemcpyDeviceToHost);
        double max_wave_speed_host = *std::max_element(max_wave_speed.begin(), max_wave_speed.end());
        double dt = cfl_number * std::min(dx, dy) / max_wave_speed_host;

        // Launch kernel to update solution
        update_solution <<<gridDim, blockDim >>> (d_Q, d_Q_new, NX, NY, dx, dy, dt, gamma);

        // Swap pointers
        std::swap(d_Q, d_Q_new);

        // Apply boundary conditions
        apply_boundary_conditions <<<gridDim, blockDim >>> (d_Q, NX, NY);

        // Copy data from device to host for output
        cudaMemcpy(Q.data(), d_Q, NX * NY * 4 * sizeof(double), cudaMemcpyDeviceToHost);

        // Write data to VTK file
        //write_vtk(Q, NX, NY, n, dx, dy, gamma);

        // Increase time
        Time += dt;

        std::cout << "Time = " << Time << std::endl;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

    std::cout << "All steps = " << n << std::endl;


    //write_vtk(Q, NX, NY, 0, dx, dy, gamma);


    // Copy data from device to host
    cudaMemcpy(Q.data(), d_Q, NX * NY * 4 * sizeof(double), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_Q);
    cudaFree(d_Q_new);
    cudaFree(d_max_wave_speed);

    return 0;
}
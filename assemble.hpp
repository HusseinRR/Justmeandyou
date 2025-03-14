#ifndef ASSEMBLE_HPP
#define ASSEMBLE_HPP

#include "Maillage.h"
#include "Matrice.hpp"
#include <cmath>

// Compute the area of a triangle defined by three points.
inline double aireTriangle(const Point &p1, const Point &p2, const Point &p3) {
    return 0.5 * std::fabs((p2.x - p1.x) * (p3.y - p1.y) -
                           (p3.x - p1.x) * (p2.y - p1.y));
}

// Compute the gradients of the local (P1) basis functions on a triangle.
inline void gradientLocal(const Point &p1, const Point &p2, const Point &p3, double grad[3][2]) {
    double A = aireTriangle(p1, p2, p3);
    double denom = 2.0 * A;
    // grad(L1) = [(y2 - y3), (x3 - x2)] / (2A)
    grad[0][0] = (p2.y - p3.y) / denom;
    grad[0][1] = (p3.x - p2.x) / denom;
    // grad(L2) = [(p3.y - p1.y), (p1.x - p3.x)] / (2A)
    grad[1][0] = (p3.y - p1.y) / denom;
    grad[1][1] = (p1.x - p3.x) / denom;
    // grad(L3) = [(p1.y - p2.y), (p2.x - p1.x)] / (2A)
    grad[2][0] = (p1.y - p2.y) / denom;
    grad[2][1] = (p2.x - p1.x) / denom;
}

// Evaluate A(x) and V(x) at a given point (x1, x2) using the covariance matrix Xi and rate r.
inline void calculer_A_V(const double Xi[2][2], double r, double x1, double x2,
                           double A_mat[2][2], double V_vec[2]) {
    A_mat[0][0] = 0.5 * Xi[0][0] * x1 * x1;
    A_mat[0][1] = 0.5 * Xi[0][1] * x1 * x2;
    A_mat[1][0] = 0.5 * Xi[1][0] * x1 * x2;
    A_mat[1][1] = 0.5 * Xi[1][1] * x2 * x2;
    
    V_vec[0] = (Xi[0][0] + 0.5 * Xi[1][0] - r) * x1;
    V_vec[1] = (Xi[1][1] + 0.5 * Xi[0][1] - r) * x2;
}

// Assemble the finite element matrices: M (mass), K (diffusion), and B (convection).
// These matrices are assembled into dense Matrice<double> objects.
inline void assemblerMatrice(const Maillage &mesh, Matrice<double> &M, 
                             Matrice<double> &K, Matrice<double> &B,
                             const double Xi[2][2], double r) {
    int N = static_cast<int>(mesh.points.size());
    M.resize(N, N);
    K.resize(N, N);
    B.resize(N, N);
    
    // Loop over all triangles in the mesh.
    for (const auto &T : mesh.triangles) {
        int i1 = T.s1, i2 = T.s2, i3 = T.s3;
        const Point &p1 = mesh.points[i1];
        const Point &p2 = mesh.points[i2];
        const Point &p3 = mesh.points[i3];
        
        double A = aireTriangle(p1, p2, p3);
        double gradL[3][2];
        gradientLocal(p1, p2, p3, gradL);
        
        // Evaluate coefficients at the triangle centroid.
        double xc = (p1.x + p2.x + p3.x) / 3.0;
        double yc = (p1.y + p2.y + p3.y) / 3.0;
        double A_mat[2][2], V_vec[2];
        calculer_A_V(Xi, r, xc, yc, A_mat, V_vec);
        
        // Local mass matrix for P1 elements: M_loc = (A/12)*[2 1 1; 1 2 1; 1 1 2]
        double facteurM = A / 12.0;
        double M_loc[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M_loc[i][j] = (i == j) ? 2.0 * facteurM : facteurM;
        
        // Local stiffness (diffusion) matrix: K_loc(i,j) = (grad L_i)^T * A_mat * (grad L_j) * A.
        double K_loc[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double gx_i = gradL[i][0], gy_i = gradL[i][1];
                double gx_j = gradL[j][0], gy_j = gradL[j][1];
                double Ax = A_mat[0][0] * gx_j + A_mat[0][1] * gy_j;
                double Ay = A_mat[1][0] * gx_j + A_mat[1][1] * gy_j;
                K_loc[i][j] = (gx_i * Ax + gy_i * Ay) * A;
            }
        }
        
        // Local convection matrix: B_loc(i,j) = (V_vec · grad L_j) * (A/3).
        double B_loc[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double gx_j = gradL[j][0], gy_j = gradL[j][1];
                double dotVg = V_vec[0] * gx_j + V_vec[1] * gy_j;
                B_loc[i][j] = dotVg * (A / 3.0);
            }
        }
        
        int loc[3] = { i1, i2, i3 };
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                M(loc[a], loc[b]) += M_loc[a][b];
                K(loc[a], loc[b]) += K_loc[a][b];
                B(loc[a], loc[b]) += B_loc[a][b];
            }
    }
}

#endif // ASSEMBLE_HPP

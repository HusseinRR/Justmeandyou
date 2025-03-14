#ifndef ASSEMBLE_HPP
#define ASSEMBLE_HPP

#include "Matrice.hpp" // Contient les fonctions aireTriangle, gradientLocal, calculer_A_V et assemblerMatrice
#include "Maillage.h"    // Définit les structures Point, Triangle et la classe Maillage
#include "vecteur.hpp"   // Autres définitions utiles
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

// Fonction pour calculer l'aire d'un triangle (déjà en français)
double aireTriangle(const Point &p1, const Point &p2, const Point &p3) {
    return 0.5 * std::fabs((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
}

// Calcule le gradient local sur l'élément triangulaire
void gradientLocal(const Point &p1, const Point &p2, const Point &p3, double grad[3][2]) {
    double A = aireTriangle(p1, p2, p3) * 2.0;
    grad[0][0] = (p2.y - p3.y) / A;  grad[0][1] = (p3.x - p2.x) / A;
    grad[1][0] = (p3.y - p1.y) / A;  grad[1][1] = (p1.x - p3.x) / A;
    grad[2][0] = (p1.y - p2.y) / A;  grad[2][1] = (p2.x - p1.x) / A;
}

// Calcule les contributions locales de la matrice A et du vecteur V
void calculer_A_V(const double Xi[2][2], double r, double x1, double x2, double A_mat[2][2], double V_vec[2]) {
    A_mat[0][0] = 0.5 * Xi[0][0] * x1 * x1;
    A_mat[0][1] = 0.5 * Xi[0][1] * x1 * x2;
    A_mat[1][0] = 0.5 * Xi[1][0] * x1 * x2;
    A_mat[1][1] = 0.5 * Xi[1][1] * x2 * x2;
    V_vec[0] = (Xi[0][0] + 0.5 * Xi[1][0] - r) * x1;
    V_vec[1] = (Xi[1][1] + 0.5 * Xi[0][1] - r) * x2;
}

// Fonction d'assemblage des matrices EF
void assemblerMatrice(const Maillage &maillage, Matrice<double> &M, Matrice<double> &K, Matrice<double> &B, const double Xi[2][2], double r) {
    size_t N = maillage.points.size();
    M = Matrice<double>(N, N);
    K = Matrice<double>(N, N);
    B = Matrice<double>(N, N);

    for (size_t it = 0; it < maillage.triangles.size(); it++) {
        const Triangle &T = maillage.triangles[it];
        const Point &p1 = maillage.points[T.s1];
        const Point &p2 = maillage.points[T.s2];
        const Point &p3 = maillage.points[T.s3];

        double A = aireTriangle(p1, p2, p3);
        double gradL[3][2];
        gradientLocal(p1, p2, p3, gradL);

        double xc = (p1.x + p2.x + p3.x) / 3.0;
        double yc = (p1.y + p2.y + p3.y) / 3.0;
        double A_mat[2][2], V_vec[2];
        calculer_A_V(Xi, r, xc, yc, A_mat, V_vec);

        Matrice<double> M_loc(3, 3), K_loc(3, 3), B_loc(3, 3);
        double facteurM = A / 12.0;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                M_loc(i, j) = (i == j) ? 2.0 * facteurM : facteurM;
                double gx_i = gradL[i][0], gy_i = gradL[i][1];
                double gx_j = gradL[j][0], gy_j = gradL[j][1];
                double Ax = A_mat[0][0] * gx_j + A_mat[0][1] * gy_j;
                double Ay = A_mat[1][0] * gx_j + A_mat[1][1] * gy_j;
                K_loc(i, j) = (gx_i * Ax + gy_i * Ay) * A;
                double produitVg = V_vec[0] * gx_j + V_vec[1] * gy_j;
                B_loc(i, j) = produitVg * (A / 3.0);
            }
        }

        int loc[3] = {T.s1, T.s2, T.s3};
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                M(loc[a], loc[b]) += M_loc(a, b);
                K(loc[a], loc[b]) += K_loc(a, b);
                B(loc[a], loc[b]) += B_loc(a, b);
            }
        }
    }
}
#endif // ASSEMBLE_HPP
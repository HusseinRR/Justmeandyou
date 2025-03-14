#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "Maillage.h"
#include "Matrice.hpp"
#include "vecteur.hpp"
#include "assemble.hpp"
#include <iomanip>

/*void afficherMatrice(const Matrice<double>& M) {
    std::cout << "🔎 Matrice " << " :" << std::endl;
    for (int i = 0; i < M.n; i++) {
        for (int j = 0; j < M.n; j++) {
            std::cout << std::setw(12) << M(i, j) << " "; // Alignement propre
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}*/

int main() {
    // 1) Définition des paramètres de l'EDP/options
    double Xi[2][2] = { {0.04, -0.024},
                        {-0.024, 0.04} };
    double r = 0.05;           // taux d'intéret
    double strike = 100.0;     // prix d'exercice (exemple)
    double T = 2.0;            // maturité en années
    double a = 200.0;          // domaine: [0, a] x [0, a]
    // Paramètres de discrétisation temporelle
    int Nt = (2*365)/5;        // nombre de pas de temps
    double dt = T / Nt;        // pas de temps

    //  2) Construction du maillage
    int Nx = 30, Ny = 30;
    Maillage mesh(0.0, a, 0.0, a, Nx, Ny);
    std::cout << "Mesh created: Npoints = " << mesh.points.size() 
              << ", Ntriangles = " << mesh.triangles.size() << std::endl;
    
    // 3) Assemblage des matrices EF
    int N = static_cast<int>(mesh.points.size());
    Matrice<double> M(N, N), Kmat(N, N), B(N, N);
    assemblerMatrice(mesh, M, Kmat, B, Xi, r);
   
    std::cout << "Matrices assembled, dimension = " << N << std::endl;
    
    // 4) Construction de D = K + B + r*M
    Matrice<double> D(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            D(i, j) = Kmat(i, j) + B(i, j) + r * M(i, j);
    
    // 5) Condition initiale Q0 = max(x1 + x2 - strike, 0)
    std::vector<double> Q0(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double x1 = mesh.points[i].x;
        double x2 = mesh.points[i].y;
        double payoff = x1 + x2 - strike;
        Q0[i] = (payoff > 0) ? payoff : 0;
    }
    
    // 6) Résolution de l'EDP en temps
    // Construction de E = M + dt*D
    Matrice<double> E(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            E(i, j) = M(i, j) + dt * D(i, j);
    
    // Factorisation de E.
    E.factorLU();
    
    // Initialisation de la solution Pk = Q0
    std::vector<double> Pk = Q0;
    
    // Boucle en temps
    for (int k = 0; k < Nt; ++k) {
        std::vector<double> Rk = M * Pk;
        // Résoudre E * P_next = Rk 
        E.solveLU(Rk);
        Pk = Rk;
        
        // Sauvegarder la solution chaque 10 pas de temps.
        if ((k + 1) % 10 == 0) {
            std::string filename = "Sol_time_" + std::to_string(k + 1) + ".dat";
            std::ofstream fout(filename);
            if (!fout) {
                std::cerr << "Error opening file " << filename << std::endl;
            } else {
                for (int i = 0; i < N; ++i)
                    fout << mesh.points[i].x << " " << mesh.points[i].y << " " << Pk[i] << "\n";
                fout.close();
                std::cout << "Wrote solution at time step " << (k + 1) << " to file " << filename << "\n";
            }
        }
    }
    
    std::cout << "Finished time-stepping." << std::endl;
    return 0;
}

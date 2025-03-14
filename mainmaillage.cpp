#include <iostream>
#include "maillage.h"

int main() {
    // Définir le domaine [a, b] x [c, d] et le nombre de subdivisions
    double a = 0.0, b = 1.0, c = 0.0, d = 1.0;
    int nx = 10, ny = 10;
    
    // Création du maillage
    Maillage mesh(a, b, c, d, nx, ny, 1);
    
    // Sauvegarde du maillage dans un fichier
    std::string filename = "maillage_output.txt";
    mesh.sauvegarderMaillage(filename);
    
    std::cout << "Maillage créé et sauvegardé dans " << filename << std::endl;
    return 0;
}
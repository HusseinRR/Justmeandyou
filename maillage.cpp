#include "Maillage.h"

Point::Point(double x_val, double y_val) : x(x_val), y(y_val) {}

Triangle::Triangle(int sommet1, int sommet2, int sommet3, int ordre_element)
    : s1(sommet1), s2(sommet2), s3(sommet3), ordre(ordre_element) {}

Maillage::Maillage(double a, double b, double c, double d, int nx_, int ny_) 
    : nx(nx_), ny(ny_) 
{
    genererPoints(a, b, c, d);
    genererTriangles();
}

void Maillage::genererPoints(double a, double b, double c, double d) {
    points.reserve((nx + 1) * (ny + 1));
    double dx = (b - a) / nx;
    double dy = (d - c) / ny;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            points.emplace_back(a + i * dx, c + j * dy);
        }
    }
}

void Maillage::genererTriangles() {
    triangles.reserve(nx * ny * 2);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int s1 = j * (nx + 1) + i;
            int s2 = s1 + 1;
            int s3 = s1 + (nx + 1);
            int s4 = s3 + 1;
            // Two triangles per cell
            triangles.emplace_back(s1, s2, s3, 1);
            triangles.emplace_back(s2, s3, s4, 1);
        }
    }
}

void Maillage::sauvegarderMaillage(const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }
    file << "Nombre de points: " << points.size() << "\n";
    file << "Nombre de triangles: " << triangles.size() << "\n\n";
    file << "Points:\n";
    for (size_t i = 0; i < points.size(); ++i)
        file << i << " " << points[i].x << " " << points[i].y << "\n";
    file << "\nTriangles:\n";
    for (size_t i = 0; i < triangles.size(); ++i)
        file << i << " " << triangles[i].s1 << " " 
             << triangles[i].s2 << " " << triangles[i].s3 
             << " " << triangles[i].ordre << "\n";
    file.close();
    std::cout << "Maillage sauvegardé dans " << filename << std::endl;
}

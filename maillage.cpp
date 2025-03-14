#include "maillage.h"

Point::Point(double x_val, double y_val) : x(x_val), y(y_val) {}

Triangle::Triangle(int sommet1, int sommet2, int sommet3, int ordre_element, int milieu1, int milieu2, int milieu3)
    : s1(sommet1), s2(sommet2), s3(sommet3), ordre(ordre_element), m1(milieu1), m2(milieu2), m3(milieu3) {}

Maillage::Maillage(double a, double b, double c, double d, int nx_, int ny_, int ordre_) 
    : nx(nx_), ny(ny_), ordre(ordre_)
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
    // Ajout des points intermédiaires si l'ordre est 2
    if (ordre == 2) {
        int base_size = points.size();
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int s1 = j * (nx + 1) + i;
                int s2 = s1 + 1;
                int s3 = s1 + (nx + 1);
                int s4 = s3 + 1;

                // Ajout des nœuds milieux
                points.emplace_back((points[s1].x + points[s2].x) / 2, (points[s1].y + points[s2].y) / 2);
                points.emplace_back((points[s2].x + points[s4].x) / 2, (points[s2].y + points[s4].y) / 2);
                points.emplace_back((points[s3].x + points[s4].x) / 2, (points[s3].y + points[s4].y) / 2);
            }
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
            // deux triangles par cellule
            if (ordre ==1) {
            triangles.emplace_back(s1, s2, s3, 1);
            triangles.emplace_back(s2, s3, s4, 1); }
        

        else {
            // Maillage d'ordre 2 (éléments quadratiques)
            int base_size = points.size();
            int m1 = base_size++;
            int m2 = base_size++;
            int m3 = base_size++;

            triangles.emplace_back(s1, s2, s3, 2, m1, m2, m3);
            triangles.emplace_back(s2, s3, s4, 2, m1, m2, m3);
        }
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
    for (size_t i = 0; i < triangles.size(); ++i) {
        file << i << " " << triangles[i].s1 << " " 
             << triangles[i].s2 << " " << triangles[i].s3 
             << " " << triangles[i].ordre << std::endl;
        if (ordre == 2){
            file << " " << triangles[i].m1 << " " << triangles[i].m2 << " " << triangles[i].m3;
        file << "\n";} }

    file.close();
    std::cout << "Maillage sauvegardé dans " << filename << std::endl;
}

#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

class Point {
public:
    double x, y;
    Point(double x_val, double y_val);
};

class Triangle {
public:
    int s1, s2, s3; // indices des sommets
    int ordre;      // ordre de l'élément 
    int m1, m2, m3; // indices des milieux
    Triangle(int sommet1, int sommet2, int sommet3, int ordre_element = 1, 
             int milieu1 = -1, int milieu2 = -1, int milieu3 = -1);
};

class Maillage {
public:
    std::vector<Point> points;      // Liste des points
    std::vector<Triangle> triangles; // Liste des triangles
    int nx, ny;                     // Nombre de subdivisions en x et en y
    int ordre; 

    // Construit un maillage sur le rectangle [a, b] x [c, d] avec nx et ny subdivisions.
    Maillage(double a, double b, double c, double d, int nx_, int ny_, int ordre_=1);
    void genererPoints(double a, double b, double c, double d);
    void genererTriangles();
    void sauvegarderMaillage(const std::string& filename);
};

#endif // MAILLAGE_H

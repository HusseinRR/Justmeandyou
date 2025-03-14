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
    int s1, s2, s3; // indices of vertices
    int ordre;      // element order (here, always 1 for linear)
    Triangle(int sommet1, int sommet2, int sommet3, int ordre_element = 1);
};

class Maillage {
public:
    std::vector<Point> points;      // List of nodes
    std::vector<Triangle> triangles; // List of triangles
    int nx, ny;                     // Number of subdivisions in x and y

    // Build a mesh over the rectangle [a, b] x [c, d] with nx, ny subdivisions.
    Maillage(double a, double b, double c, double d, int nx_, int ny_);
    void genererPoints(double a, double b, double c, double d);
    void genererTriangles();
    void sauvegarderMaillage(const std::string& filename);
};

#endif // MAILLAGE_H

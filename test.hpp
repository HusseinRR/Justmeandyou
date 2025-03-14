#ifndef SIMPLE_TESTS_HPP
#define SIMPLE_TESTS_HPP

#include "Maillage.h"
#include "Matrice.hpp"
#include "assemble.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// Utility: print a matrix.
void printMatrice(const Matrice<double>& M) {
    int n = M.nbLignes();
    int m = M.nbColonnes();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cout << std::setw(12) << M(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

// Utility: print a vector.
void printVector(const std::vector<double>& v) {
    for (size_t i = 0; i < v.size(); i++) {
        std::cout << "v[" << i << "] = " << v[i] << std::endl;
    }
}

// Test 1: K matrix multiplied by ones should be nearly zero (constant function is in the null space)
void testStiffnessWithOnes() {
    std::cout << "Test: K matrix * ones (constant function)" << std::endl;
    
    Maillage mesh;
    mesh.points.push_back(Point(0.0, 0.0));
    mesh.points.push_back(Point(1.0, 0.0));
    mesh.points.push_back(Point(0.0, 1.0));
    mesh.triangles.push_back(Triangle(0, 1, 2));
    
    // Identity for Xi and zero reaction for a pure Laplacian problem.
    double Xi[2][2] = { {1.0, 0.0}, {0.0, 1.0} };
    double r = 0.0;
    
    Matrice<double> M(3, 3), K(3, 3), B(3, 3);
    assemblerMatrice(mesh, M, K, B, Xi, r);
    
    std::vector<double> ones(3, 1.0);
    std::vector<double> result = K * ones;
    
    std::cout << "Computed values for K * ones:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "result[" << i << "] = " << result[i] << " (expected: 0)" << std::endl;
    }
    
    std::cout << "------------------------------" << std::endl;
}

// Test 2: M matrix integrated against a constant function should give area/3 per node
void testMassMatrixIntegration() {
    std::cout << "Test: M matrix * ones (mass matrix integration)" << std::endl;
    
    Maillage mesh;
    mesh.points.push_back(Point(0.0, 0.0));
    mesh.points.push_back(Point(1.0, 0.0));
    mesh.points.push_back(Point(0.0, 1.0));
    mesh.triangles.push_back(Triangle(0, 1, 2));
    
    double Xi[2][2] = { {1.0, 0.0}, {0.0, 1.0} };
    double r = 0.0;
    
    Matrice<double> M(3, 3), K(3, 3), B(3, 3);
    assemblerMatrice(mesh, M, K, B, Xi, r);
    
    std::vector<double> ones(3, 1.0);
    std::vector<double> result = M * ones;
    
    // For a triangle with vertices (0,0), (1,0), (0,1) the area is 0.5.
    // Integration of each linear basis function gives area/3.
    double expected = 0.5 / 3.0;
    
    std::cout << "Computed values for M * ones:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "result[" << i << "] = " << result[i]
                  << " (expected: " << expected << ")" << std::endl;
    }
    
    std::cout << "------------------------------" << std::endl;
}

// Test 3: Simple LU factorization with a 2x2 system
void testSimpleLU() {
    std::cout << "Test: Simple LU factorization (2x2 system)" << std::endl;
    
    Matrice<double> A(2, 2);
    A(0, 0) = 2.0;  A(0, 1) = 1.0;
    A(1, 0) = 1.0;  A(1, 1) = 3.0;
    
    std::vector<double> b = {3.0, 4.0};
    
    A.factoriserLU();
    A.resoudreLU(b);
    
    std::vector<double> expected = {1.0, 1.0};
    
    std::cout << "Computed solution:" << std::endl;
    for (size_t i = 0; i < b.size(); i++) {
        std::cout << "b[" << i << "] = " << b[i]
                  << " (expected: " << expected[i] << ")" << std::endl;
    }
    
    std::cout << "------------------------------" << std::endl;
}

// Test 4: K matrix applied to a nonconstant affine function (f(x,y)=x)
// For the triangle, the nodal values for f(x,y)=x are: f(0,0)=0, f(1,0)=1, f(0,1)=0.
void testAffineStiffness() {
    std::cout << "Test: K matrix * affine function (f(x,y)=x)" << std::endl;
    
    Maillage mesh;
    mesh.points.push_back(Point(0.0, 0.0)); // f=0
    mesh.points.push_back(Point(1.0, 0.0)); // f=1
    mesh.points.push_back(Point(0.0, 1.0)); // f=0
    mesh.triangles.push_back(Triangle(0, 1, 2));
    
    double Xi[2][2] = { {1.0, 0.0}, {0.0, 1.0} };
    double r = 0.0;
    
    Matrice<double> M(3, 3), K(3, 3), B(3, 3);
    assemblerMatrice(mesh, M, K, B, Xi, r);
    
    std::vector<double> affineFunction = {0.0, 1.0, 0.0};
    std::vector<double> result = K * affineFunction;
    
    std::cout << "Computed values for K * affine function:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "result[" << i << "] = " << result[i] << " (expected: 0)" << std::endl;
    }
    
    std::cout << "------------------------------" << std::endl;
}

#endif // SIMPLE_TESTS_HPP

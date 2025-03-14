#ifndef MATRICE_HPP
#define MATRICE_HPP

#include <vector>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iostream>

template<typename T>
class Matrice {
public:
    int m, n;           // number of rows and columns
    std::vector<T> data; // storage in row–major order
public:
    Matrice();
    Matrice(int m, int n);
    void resize(int m, int n);
    int rows() const;
    int cols() const;
    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    // Basic arithmetic operations
    Matrice<T> operator+(const Matrice<T>& other) const;
    Matrice<T>& operator+=(const Matrice<T>& other);
    Matrice<T> operator-(const Matrice<T>& other) const;
    Matrice<T>& operator-=(const Matrice<T>& other);
    Matrice<T> operator*(const T& scalar) const;
    Matrice<T>& operator*=(const T& scalar);
    Matrice<T> operator/(const T& scalar) const;
    Matrice<T>& operator/=(const T& scalar);

    // Matrix multiplication
    Matrice<T> operator*(const Matrice<T>& other) const;
    // Matrix-vector multiplication (with std::vector)
    std::vector<T> operator*(const std::vector<T>& vec) const;

    // LU factorization (in-place, no pivoting)
    void factorLU();
    // Solve LU * x = b; b is overwritten by the solution x.
    void solveLU(std::vector<T>& b) const;

    // For debugging: return a string representation
    std::string str() const;
};

template<typename T>
Matrice<T>::Matrice() : m(0), n(0) {}

template<typename T>
Matrice<T>::Matrice(int m_, int n_) : m(m_), n(n_), data(m_*n_, T(0)) {}

template<typename T>
void Matrice<T>::resize(int m_, int n_) {
    m = m_;
    n = n_;
    data.assign(m * n, T(0));
}

template<typename T>
int Matrice<T>::rows() const { return m; }

template<typename T>
int Matrice<T>::cols() const { return n; }

template<typename T>
T& Matrice<T>::operator()(int i, int j) {
    if (i < 0 || i >= m || j < 0 || j >= n)
        throw std::out_of_range("Index out of range in Matrice");
    return data[i * n + j];
}

template<typename T>
const T& Matrice<T>::operator()(int i, int j) const {
    if (i < 0 || i >= m || j < 0 || j >= n)
        throw std::out_of_range("Index out of range in Matrice");
    return data[i * n + j];
}

template<typename T>
Matrice<T> Matrice<T>::operator+(const Matrice<T>& other) const {
    if (m != other.m || n != other.n)
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    Matrice<T> result(m, n);
    for (int i = 0; i < m * n; ++i)
        result.data[i] = data[i] + other.data[i];
    return result;
}

template<typename T>
Matrice<T>& Matrice<T>::operator+=(const Matrice<T>& other) {
    if (m != other.m || n != other.n)
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    for (int i = 0; i < m * n; ++i)
        data[i] += other.data[i];
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator-(const Matrice<T>& other) const {
    if (m != other.m || n != other.n)
        throw std::invalid_argument("Matrix dimensions must match for subtraction.");
    Matrice<T> result(m, n);
    for (int i = 0; i < m * n; ++i)
        result.data[i] = data[i] - other.data[i];
    return result;
}

template<typename T>
Matrice<T>& Matrice<T>::operator-=(const Matrice<T>& other) {
    if (m != other.m || n != other.n)
        throw std::invalid_argument("Matrix dimensions must match for subtraction.");
    for (int i = 0; i < m * n; ++i)
        data[i] -= other.data[i];
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator*(const T& scalar) const {
    Matrice<T> result(m, n);
    for (int i = 0; i < m * n; ++i)
        result.data[i] = data[i] * scalar;
    return result;
}

template<typename T>
Matrice<T>& Matrice<T>::operator*=(const T& scalar) {
    for (int i = 0; i < m * n; ++i)
        data[i] *= scalar;
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator/(const T& scalar) const {
    Matrice<T> result(m, n);
    for (int i = 0; i < m * n; ++i)
        result.data[i] = data[i] / scalar;
    return result;
}

template<typename T>
Matrice<T>& Matrice<T>::operator/=(const T& scalar) {
    for (int i = 0; i < m * n; ++i)
        data[i] /= scalar;
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator*(const Matrice<T>& other) const {
    if (n != other.m)
        throw std::invalid_argument("Matrix multiplication dimension mismatch.");
    Matrice<T> result(m, other.n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < other.n; ++j) {
            T sum = T(0);
            for (int k = 0; k < n; ++k)
                sum += (*this)(i, k) * other(k, j);
            result(i, j) = sum;
        }
    }
    return result;
}

template<typename T>
std::vector<T> Matrice<T>::operator*(const std::vector<T>& vec) const {
    if (vec.size() != static_cast<size_t>(n))
        throw std::invalid_argument("Matrix-vector multiplication dimension mismatch.");
    std::vector<T> result(m, T(0));
    for (int i = 0; i < m; ++i) {
        T sum = T(0);
        for (int j = 0; j < n; ++j)
            sum += (*this)(i, j) * vec[j];
        result[i] = sum;
    }
    return result;
}

template<typename T>
void Matrice<T>::factorLU() {
    if (m != n)
        throw std::invalid_argument("LU factorization requires a square matrix.");
    for (int k = 0; k < m; ++k) {
        if (std::fabs((*this)(k, k)) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in LU factorization.");
        for (int i = k + 1; i < m; ++i) {
            (*this)(i, k) /= (*this)(k, k);
            for (int j = k + 1; j < m; ++j)
                (*this)(i, j) -= (*this)(i, k) * (*this)(k, j);
        }
    }
}

template<typename T>
void Matrice<T>::solveLU(std::vector<T>& b) const {
    if (m != n || b.size() != static_cast<size_t>(m))
        throw std::invalid_argument("LU solver requires square matrix and matching vector size.");
    std::vector<T> x(b); // copy b into x
    // Forward substitution for L (L has unit diagonal)
    for (int i = 1; i < m; ++i) {
        T sum = x[i];
        for (int j = 0; j < i; ++j)
            sum -= (*this)(i, j) * x[j];
        x[i] = sum; // no division since L(i,i)=1
    }
    // Back substitution for U
    for (int i = m - 1; i >= 0; --i) {
        T sum = x[i];
        for (int j = i + 1; j < m; ++j)
            sum -= (*this)(i, j) * x[j];
        if (std::fabs((*this)(i, i)) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in back substitution.");
        x[i] = sum / (*this)(i, i);
    }
    b = x;
}

template<typename T>
std::string Matrice<T>::str() const {
    std::ostringstream oss;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            oss << (*this)(i, j) << " ";
        oss << "\n";
    }
    return oss.str();
}

#endif // MATRICE_HPP

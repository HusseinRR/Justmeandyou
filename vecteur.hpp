#ifndef VECTEUR_HPP
#define VECTEUR_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

template<typename T>
std::vector<T> operator+(const std::vector<T>& u, const std::vector<T>& v) {
    if (u.size() != v.size())
        throw std::invalid_argument("Vector addition: sizes do not match.");
    std::vector<T> result(u.size());
    for (size_t i = 0; i < u.size(); ++i)
        result[i] = u[i] + v[i];
    return result;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& u, const std::vector<T>& v) {
    if (u.size() != v.size())
        throw std::invalid_argument("Vector subtraction: sizes do not match.");
    std::vector<T> result(u.size());
    for (size_t i = 0; i < u.size(); ++i)
        result[i] = u[i] - v[i];
    return result;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& u, const T& s) {
    std::vector<T> result(u.size());
    for (size_t i = 0; i < u.size(); ++i)
        result[i] = u[i] * s;
    return result;
}

template<typename T>
std::vector<T> operator*(const T& s, const std::vector<T>& u) {
    return u * s;
}

template<typename T>
std::vector<T> operator/(const std::vector<T>& u, const T& s) {
    std::vector<T> result(u.size());
    for (size_t i = 0; i < u.size(); ++i)
        result[i] = u[i] / s;
    return result;
}

template<typename T>
T dot(const std::vector<T>& u, const std::vector<T>& v) {
    if (u.size() != v.size())
        throw std::invalid_argument("Dot product: sizes do not match.");
    T sum = T(0);
    for (size_t i = 0; i < u.size(); ++i)
        sum += u[i] * v[i];
    return sum;
}

template<typename T>
T norme(const std::vector<T>& u) {
    return std::sqrt(dot(u, u));
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "(";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1)
            os << ", ";
    }
    os << ")";
    return os;
}

#endif // VECTEUR_HPP

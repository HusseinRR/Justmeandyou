#ifndef MATRICE_HPP
#define MATRICE_HPP

#include <vector>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>

// ============================================================================
// Classe de matrice dense
// ============================================================================
template<typename T>
class Matrice {
public:
    int lignes, colonnes;         // nombre de lignes et de colonnes
    std::vector<T> donnees;       // stockage en ordre ligne-major

public:
    Matrice();
    Matrice(int lignes, int colonnes);
    void redimensionner(int lignes, int colonnes);
    int nbLignes() const;
    int nbColonnes() const;
    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    // Opérations arithmétiques de base
    Matrice<T> operator+(const Matrice<T>& autre) const;
    Matrice<T>& operator+=(const Matrice<T>& autre);
    Matrice<T> operator-(const Matrice<T>& autre) const;
    Matrice<T>& operator-=(const Matrice<T>& autre);
    Matrice<T> operator*(const T& scalaire) const;
    Matrice<T>& operator*=(const T& scalaire);
    Matrice<T> operator/(const T& scalaire) const;
    Matrice<T>& operator/=(const T& scalaire);

    // Multiplication de matrices et multiplication matrice-vecteur
    Matrice<T> operator*(const Matrice<T>& autre) const;
    std::vector<T> operator*(const std::vector<T>& vect) const;

    // Factorisation LU (in-situ, sans pivotage)
    void factoriserLU();
    // Résolution de LU * x = b ; b est remplacé par la solution x.
    void resoudreLU(std::vector<T>& b) const;

    // Pour débogage : retourne une représentation sous forme de chaîne
    std::string toString() const;
};

template<typename T>
Matrice<T>::Matrice() : lignes(0), colonnes(0) {}

template<typename T>
Matrice<T>::Matrice(int lignes_, int colonnes_) : lignes(lignes_), colonnes(colonnes_), donnees(lignes_ * colonnes_, T(0)) {}

template<typename T>
void Matrice<T>::redimensionner(int lignes_, int colonnes_) {
    lignes = lignes_;
    colonnes = colonnes_;
    donnees.assign(lignes * colonnes, T(0));
}

template<typename T>
int Matrice<T>::nbLignes() const { return lignes; }

template<typename T>
int Matrice<T>::nbColonnes() const { return colonnes; }

template<typename T>
T& Matrice<T>::operator()(int i, int j) {
    return donnees[i * colonnes + j];
}

template<typename T>
const T& Matrice<T>::operator()(int i, int j) const {
    return donnees[i * colonnes + j];
}

template<typename T>
Matrice<T> Matrice<T>::operator+(const Matrice<T>& autre) const {
    Matrice<T> resultat(lignes, colonnes);
    for (int i = 0; i < lignes * colonnes; ++i)
        resultat.donnees[i] = donnees[i] + autre.donnees[i];
    return resultat;
}

template<typename T>
Matrice<T>& Matrice<T>::operator+=(const Matrice<T>& autre) {
    for (int i = 0; i < lignes * colonnes; ++i)
        donnees[i] += autre.donnees[i];
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator-(const Matrice<T>& autre) const {
    Matrice<T> resultat(lignes, colonnes);
    for (int i = 0; i < lignes * colonnes; ++i)
        resultat.donnees[i] = donnees[i] - autre.donnees[i];
    return resultat;
}

template<typename T>
Matrice<T>& Matrice<T>::operator-=(const Matrice<T>& autre) {
    for (int i = 0; i < lignes * colonnes; ++i)
        donnees[i] -= autre.donnees[i];
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator*(const T& scalaire) const {
    Matrice<T> resultat(lignes, colonnes);
    for (int i = 0; i < lignes * colonnes; ++i)
        resultat.donnees[i] = donnees[i] * scalaire;
    return resultat;
}

template<typename T>
Matrice<T>& Matrice<T>::operator*=(const T& scalaire) {
    for (int i = 0; i < lignes * colonnes; ++i)
        donnees[i] *= scalaire;
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator/(const T& scalaire) const {
    Matrice<T> resultat(lignes, colonnes);
    for (int i = 0; i < lignes * colonnes; ++i)
        resultat.donnees[i] = donnees[i] / scalaire;
    return resultat;
}

template<typename T>
Matrice<T>& Matrice<T>::operator/=(const T& scalaire) {
    for (int i = 0; i < lignes * colonnes; ++i)
        donnees[i] /= scalaire;
    return *this;
}

template<typename T>
Matrice<T> Matrice<T>::operator*(const Matrice<T>& autre) const {
    Matrice<T> resultat(lignes, autre.colonnes);
    for (int i = 0; i < lignes; ++i) {
        for (int j = 0; j < autre.colonnes; ++j) {
            T somme = T(0);
            for (int k = 0; k < colonnes; ++k)
                somme += (*this)(i, k) * autre(k, j);
            resultat(i, j) = somme;
        }
    }
    return resultat;
}

template<typename T>
std::vector<T> Matrice<T>::operator*(const std::vector<T>& vect) const {
    std::vector<T> resultat(lignes, T(0));
    for (int i = 0; i < lignes; ++i) {
        T somme = T(0);
        for (int j = 0; j < colonnes; ++j)
            somme += (*this)(i, j) * vect[j];
        resultat[i] = somme;
    }
    return resultat;
}

template<typename T>
void Matrice<T>::factoriserLU() {
    for (int k = 0; k < lignes; ++k) {
        for (int i = k + 1; i < lignes; ++i) {
            (*this)(i, k) /= (*this)(k, k);
            for (int j = k + 1; j < lignes; ++j)
                (*this)(i, j) -= (*this)(i, k) * (*this)(k, j);
        }
    }
}

template<typename T>
void Matrice<T>::resoudreLU(std::vector<T>& b) const {
    std::vector<T> x(b); // copie de b dans x
    // Substitution avant pour L (diagonale unitaire)
    for (int i = 1; i < lignes; ++i) {
        T somme = x[i];
        for (int j = 0; j < i; ++j)
            somme -= (*this)(i, j) * x[j];
        x[i] = somme;
    }
    // Substitution arrière pour U
    for (int i = lignes - 1; i >= 0; --i) {
        T somme = x[i];
        for (int j = i + 1; j < lignes; ++j)
            somme -= (*this)(i, j) * x[j];
        if (std::fabs((*this)(i, i)) < 1e-12)
            throw std::runtime_error("Pivot nul rencontré lors de la substitution arrière.");
        x[i] = somme / (*this)(i, i);
    }
    b = x;
}

template<typename T>
std::string Matrice<T>::toString() const {
    std::ostringstream oss;
    for (int i = 0; i < lignes; ++i) {
        for (int j = 0; j < colonnes; ++j)
            oss << (*this)(i, j) << " ";
        oss << "\n";
    }
    return oss.str();
}

// ============================================================================
// Classe pour matrice symétrique (stockage de la partie triangulaire inférieure)
// ============================================================================
template<typename T>
class MatriceSymetrique {
public:
    int taille;                    // dimension de la matrice carrée
    std::vector<T> donnees;        // stockage de la partie triangulaire inférieure

public:
    MatriceSymetrique();
    MatriceSymetrique(int taille);
    void redimensionner(int taille);
    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    // Opérations basiques
    MatriceSymetrique<T> operator+(const MatriceSymetrique<T>& autre) const;
    MatriceSymetrique<T>& operator+=(const MatriceSymetrique<T>& autre);
    MatriceSymetrique<T> operator*(const T& scalaire) const;
    MatriceSymetrique<T>& operator*=(const T& scalaire);

    // Conversion en matrice dense
    Matrice<T> toDense() const;

    // Factorisation de Cholesky (pour matrices symétriques définies positives)
    void factoriserCholesky();
    // Résolution de L*L^T * x = b, b est remplacé par la solution x.
    void resoudreCholesky(std::vector<T>& b) const;
};

template<typename T>
MatriceSymetrique<T>::MatriceSymetrique() : taille(0) {}

template<typename T>
MatriceSymetrique<T>::MatriceSymetrique(int taille_) : taille(taille_), donnees(taille_ * (taille_ + 1) / 2, T(0)) {}

template<typename T>
void MatriceSymetrique<T>::redimensionner(int taille_) {
    taille = taille_;
    donnees.assign(taille * (taille + 1) / 2, T(0));
}

template<typename T>
T& MatriceSymetrique<T>::operator()(int i, int j) {
    if (i < j) std::swap(i, j);
    return donnees[i * (i + 1) / 2 + j];
}

template<typename T>
const T& MatriceSymetrique<T>::operator()(int i, int j) const {
    if (i < j) std::swap(i, j);
    return donnees[i * (i + 1) / 2 + j];
}

template<typename T>
MatriceSymetrique<T> MatriceSymetrique<T>::operator+(const MatriceSymetrique<T>& autre) const {
    MatriceSymetrique<T> resultat(taille);
    int nb = taille * (taille + 1) / 2;
    for (int i = 0; i < nb; ++i)
        resultat.donnees[i] = donnees[i] + autre.donnees[i];
    return resultat;
}

template<typename T>
MatriceSymetrique<T>& MatriceSymetrique<T>::operator+=(const MatriceSymetrique<T>& autre) {
    int nb = taille * (taille + 1) / 2;
    for (int i = 0; i < nb; ++i)
        donnees[i] += autre.donnees[i];
    return *this;
}

template<typename T>
MatriceSymetrique<T> MatriceSymetrique<T>::operator*(const T& scalaire) const {
    MatriceSymetrique<T> resultat(taille);
    int nb = taille * (taille + 1) / 2;
    for (int i = 0; i < nb; ++i)
        resultat.donnees[i] = donnees[i] * scalaire;
    return resultat;
}

template<typename T>
MatriceSymetrique<T>& MatriceSymetrique<T>::operator*=(const T& scalaire) {
    int nb = taille * (taille + 1) / 2;
    for (int i = 0; i < nb; ++i)
        donnees[i] *= scalaire;
    return *this;
}

template<typename T>
Matrice<T> MatriceSymetrique<T>::toDense() const {
    Matrice<T> dense(taille, taille);
    for (int i = 0; i < taille; ++i)
        for (int j = 0; j < taille; ++j)
            dense(i, j) = (i >= j ? donnees[i*(i+1)/2 + j] : donnees[j*(j+1)/2 + i]);
    return dense;
}

template<typename T>
void MatriceSymetrique<T>::factoriserCholesky() {
    for (int i = 0; i < taille; ++i) {
        T somme = (*this)(i, i);
        for (int k = 0; k < i; ++k)
            somme -= (*this)(i, k) * (*this)(i, k);
        if (somme <= 0)
            throw std::runtime_error("Matrice non définie positive pour la factorisation de Cholesky.");
        (*this)(i, i) = std::sqrt(somme);
        for (int j = i + 1; j < taille; ++j) {
            T somme2 = (*this)(j, i);
            for (int k = 0; k < i; ++k)
                somme2 -= (*this)(j, k) * (*this)(i, k);
            (*this)(j, i) = somme2 / (*this)(i, i);
        }
    }
}

template<typename T>
void MatriceSymetrique<T>::resoudreCholesky(std::vector<T>& b) const {
    if (b.size() != static_cast<size_t>(taille))
        throw std::invalid_argument("La taille du vecteur ne correspond pas à la dimension de la matrice.");
    std::vector<T> y(taille, T(0));
    // Résolution de L * y = b
    for (int i = 0; i < taille; ++i) {
        T somme = b[i];
        for (int k = 0; k < i; ++k)
            somme -= (*this)(i, k) * y[k];
        y[i] = somme / (*this)(i, i);
    }
    std::vector<T> x(taille, T(0));
    // Résolution de L^T * x = y
    for (int i = taille - 1; i >= 0; --i) {
        T somme = y[i];
        for (int k = i + 1; k < taille; ++k)
            somme -= (*this)(k, i) * x[k];
        x[i] = somme / (*this)(i, i);
    }
    b = x;
}

// ============================================================================
// Classe pour matrice creuse (sparse) en format COO
// ============================================================================
template<typename T>
class MatriceCreuse {
public:
    int lignes, colonnes;
    std::vector<int> lignesIndices;    // indices de ligne
    std::vector<int> colonnesIndices;    // indices de colonne
    std::vector<T> valeurs;            // valeurs non nulles

public:
    MatriceCreuse();
    MatriceCreuse(int lignes, int colonnes);
    void redimensionner(int lignes, int colonnes);

    // Ajoute une valeur à la position (i, j) (la valeur est ajoutée à l'éventuelle valeur existante)
    void ajouter(int i, int j, T valeur);

    // Conversion en matrice dense
    Matrice<T> toDense() const;
};

template<typename T>
MatriceCreuse<T>::MatriceCreuse() : lignes(0), colonnes(0) {}

template<typename T>
MatriceCreuse<T>::MatriceCreuse(int lignes_, int colonnes_) : lignes(lignes_), colonnes(colonnes_) {}

template<typename T>
void MatriceCreuse<T>::redimensionner(int lignes_, int colonnes_) {
    lignes = lignes_;
    colonnes = colonnes_;
    lignesIndices.clear();
    colonnesIndices.clear();
    valeurs.clear();
}

template<typename T>
void MatriceCreuse<T>::ajouter(int i, int j, T valeur) {
    if (i < 0 || i >= lignes || j < 0 || j >= colonnes)
        throw std::out_of_range("Indice hors limites dans MatriceCreuse");
    lignesIndices.push_back(i);
    colonnesIndices.push_back(j);
    valeurs.push_back(valeur);
}

template<typename T>
Matrice<T> MatriceCreuse<T>::toDense() const {
    Matrice<T> dense(lignes, colonnes);
    for (size_t k = 0; k < valeurs.size(); ++k) {
        dense(lignesIndices[k], colonnesIndices[k]) += valeurs[k];
    }
    return dense;
}

#endif // MATRICE_HPP

//g++ main_gabriel.cpp vecteur.hpp Matrice.hpp -o main_gabriel.exe -Wall ; .\main_gabriel.exe
#include "vecteur.hpp"
#include "Matrice.hpp"

using namespace std;

int main(){
    cout << "Test en reel\n";
    /*
    //test de la classe vector
    vector V={1,0};
    vector W={0,1};
    cout<<"V = "<<V<<"\n";
    cout<<"W = "<<W<<"\n";
    cout<<"V+W = "<<V+W<<"\n";
    cout<<"V-W = "<<V-W<<"\n";
    double X=V|W;
    cout<<"V*W = "<<X<<"\n";*/

    //test de la classe matrice
/*
    Matrice<double> J(3,2);
    cout<<"Matrice J = "<<J;
    J.affichage();
    int n=10;
    Matrice<double> M(n,n);
    for(int i=0;i<n;i++)
    {M(i,i)=2.;
     if(i>0) M(i,i-1)=-1;
     if(i<n) M(i,i+1)=-1;
    }
    M.affichage();
    Matrice<double> I(n,n);
    for(int i=0;i<n;i++) I(i,i)=1.;
    Matrice<double> operation=(M*3.-I*2.)/2.;
    cout<<"(3.*M-2*I)/2"<<operation<<"\n";
    operation.affichage();
    cout<<"nouvelle matrice\n";
    Matrice<double> T(n,n+1);
    for(int i=0;i<n;i++) T(i,n)=i+1;// on accède avec i de 0 a n
    T.affichage();
    T.supprime(2,10);// colonnes entre 0 et m-1 puis lignes entre 0 et n-1
    T.affichage();
    Matrice<double> MT=M+T;
    MT.affichage();*/

    //tester LU

    Matrice<double> A(3, 3);
    A(0,0) = 2; A(0,1) = -1; A(0,2) = 1;
    A(1,0) = 3; A(1,1) = 3;  A(1,2) = 9;
    A(2,0) = 3; A(2,1) = 3;  A(2,2) = 5;    
    cout<<A;
    cout << "Matrice A :\n";
    //cout << "Nombre de coefficients stockes : " << A.coeffs.size() << endl;
   
    A.affichage();

    // Décomposition LU
    Matrices_L_U<double> LU = A.LUdecomposition();

    cout << "\nMatrice L :\n";
    LU.L.affichage();

    cout << "\nMatrice U :\n";
    LU.U.affichage();
    Matrice<double> P=LU.L*LU.U;
    P.affichage();

    cout<<" Solveur LUX=b :\n";
    vector<double> b(3, 1);
    b[0] = 8;
    b[1] = 171;
    b[2] = 12;
    
    // Résolution du système AX = b
    Matrice<double> X = A.solveurLU(b);

    

    // Affichage de la solution X
    cout << "Solution X : \n";
    X.affichage();
    //cout<<" comparaison : \n"<<X<<A*X<<b;
    
    vector<double> x(3, 1);
    x[0] = 13;
    x[1] = 17;
    x[2] = 24;

    cout<<"Nouveau test : b prime ="<<A*x;

    vector<double> bprime(3, 1);
    bprime[0] = 33;
    bprime[1] = 306;
    bprime[2] = 210;
    
    // Résolution du système AX = b
    Matrice<double> Xprime = A.solveurLU(bprime);
    cout << "Solution Xprime: \n";
    Xprime.affichage();
    

}
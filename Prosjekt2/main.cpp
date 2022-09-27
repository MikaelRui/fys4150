//g++ main.cpp -o main.exe -O2 -larmadillo -llapack -lblas

#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>


using namespace arma;
using namespace std;

int k = 0;
int l = 1;
double Eps = pow(10,-8);
double tau;
double t;
double s;
double c;
mat R;


///funksjon som tar inn en matrise og returnerer høyeste absoluttverdien.
double max_offdiag_symmetric(mat A, int n, int k, int l);

double max_offdiag_symmetric(mat A, int n, int k, int l) {
    double maxValue;
    for (int i = 0; i < n; ++i){
        for (int j = i+1; j < n; ++j){
            double aij = fabs(A(i,j));
            if (aij > maxValue){
                maxValue = aij; k=i; l=j;
            }
        }
    }
    cout << "maxverdi: " << maxValue << "index" << "[" << k << "," << l << "]" << endl;
    return maxValue;
}
///Jacobi's rotation algorithm
mat JacobiRotation(mat A, mat R, int n, int k, int l);
mat JacobiRotation(mat A, mat R, int n, int k, int l){

    mat Q;
    while (abs(A(k,l)) > Eps) {
        ///3.1
        double tau = (A(l,l) - A(k,k))/(2*A(k,k));
        ///3.2
        if (tau > 0){
            t = - tau + sqrt(1 + pow(tau,2));
        }
        if (tau < 0){
            t = (-1) / (- tau + sqrt(1 + pow(tau,2)));
        }
        c = 1/(sqrt(1 + pow(tau,2)));
        s = c*t;
        ///3.3 keeping track of elements, so i don't use the new value
        double a_kk = A(k,k);
        double a_ll = A(l,l);

        A(k,k) = A(k,k) *c*c -2 * A(k,l)*c*s +a_ll * s*s;
        A(l,l) = A(l,l)*c*c + 2*A(k,l)*c*s + a_kk*s*s;
        A(k,l) = 0;
        A(l,k) = 0;


        for (int i =0; i <n; i++){
            if (i !=k && i!= l) {
                double a_ik = A(i,k);
                double a_il = A(i,l);
                A(i,k) = a_ik * c - a_il * s;
                A(k,i) = A(i,k); /// symmetric matrix
                A(i,l) = a_il * c + a_ik * s;
                A(l,i) = A(i,l); /// symmetric matrix
            }
            ///3.4 update the overall rotation matrix for all i, also keeping track of elements here
            double r_ik = R(i,k);
            double r_il = R(i,l);

            R(i,k) = r_ik *c - r_il * s;
            R(i,l) = r_il * c - r_ik *s;

        }


    }
    ///3.5 find indices k and l for the new max off-diag element
    double maxV;
    for (int i = 0; i < n; ++i){
        for (int j = i+1; j < n; ++j){
            double aij = fabs(A(i,j));
            if (aij > maxV){
                maxV = aij; k=i; l=j;
            }
        }
    }
    if (abs(A(k,l))< Eps){
        cout << "A matrix:" << endl;
        cout << A << endl;
        cout << "R matrix:" << endl;
        cout << R << endl;
        }
    return JacobiRotation(A,R,n,k,l);
}

int main() {
    ///Problem 2
    ///making the tridiag matrix
    std::cout << "Hello, World!" << std::endl;
    int h= 6+1;
    mat A(6,6); A.zeros();
    A.diag(0) += (2/ pow(h,2));
    A.diag(1) += (-1/ pow(h,2));
    A.diag(-1) += (-1/ pow(h,2));

    cout << A;

vec eigval;
mat eigvec;
mat B;
B = A.t()*A;

eig_sym(eigval, eigvec,B);

cout << "eigenvalues:" << endl << eigval;
eigval.save("eigenvalues.dat", raw_ascii);
cout << "eigenvectors:" << endl << eigvec;
eigvec.save("eigenvector.dat", raw_ascii);

/// problem 3
///a)

max_offdiag_symmetric(eigvec, 6, k, l);



///3b)


mat b;
b << 1 << 0 << 0 << 0.5 << endr
  << 0 << 1 << -0.7 << 0 << endr
  << 0 << -0.7 << 1 << 0 << endr
  << 0.5 << 0 << 0 << 1 << endr;
cout << b << endl;


///setter matrise b inn i funksjonen for å finne høyeste off-diag verdien
max_offdiag_symmetric(b,4, k, l);

///Problem 4: Jacobi's rotation algorithm
///a, b)
mat R(6,6); R.zeros(); R.diag(0) += (1.);

max_offdiag_symmetric(eigvec, 6, 0, 1);

JacobiRotation(eigvec, R, 6, k, l);
cout << "A" << eigvec << "R" << R << endl;

///problem 5
///a)
///new R-matrix called G
mat G(10,10); G.zeros(); G.diag(0) += (1.0);
mat C = mat(10,10).randn();
C = symmatu(C);

auto t1 = std::chrono::high_resolution_clock::now();
JacobiRotation(C,G,10,0,1);
auto t2 = std::chrono::high_resolution_clock::now();

return 0;
}


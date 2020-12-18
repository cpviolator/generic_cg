
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <sys/time.h>

#include <Eigen/Eigenvalues>

using namespace std;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;

typedef std::complex<double> Complex;

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define I Complex(0,1.0)
#define cUnit Complex(1.0,0)

#include <cg.h>

int main(int argc, char **argv) {

  // Problem parameters 
  int N = 64;
  int maxiter = 100000;
  double tol = 1e-10;
  double diag = N/2;
  bool laplace_mat = false;
  
  // Starting guess
  std::vector<Complex> x(N, 0.0);
  // vector to solve
  std::vector<Complex> b(N, 1.0);  

  // The row major matrix
  std::vector<std::vector<Complex>> mat(N, std::vector<Complex> (N, 0.0));  
  
  // Use Eigen to generate a random matrix
  MatrixXcd ref = MatrixXcd::Random(N, N);
  
  // Copy the Eigen matrix into the array
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++) {
      mat[i][j] = (laplace_mat ? 0.0 : ref(i,j) + conj(ref(j, i)));
      if(i == j) mat[i][j] += (laplace_mat ? 2 : diag);
      if((i == j+1 || j == i+1) && laplace_mat) mat[i][j] = -1.0;
    }
  
  // Pass the row major matrix, guess, source, and constraints to CG
  normalise(b);
  cg(mat, x, b, tol, maxiter);
  
  // Test for correcteness
  std::vector<Complex> test(N, 0.0);
  matVec(mat, test, x);

  // The 'test' vect should be the same as b
  cout << "Source norm = " << norm(b) << endl;
  axpy(-1.0, test, b);
  cout << "Residual norm = " << norm(b) << endl;
  cout << "Solution norm = " << norm(x) << endl;  
}

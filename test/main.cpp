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

#include <cg.h>

int main(int argc, char **argv) {

  // Problem parameters
  int N = 1024; // Will use an NxN matrix
  
  int maxiter = 100000; // The maximul allowed CG iterations
  
  double tol = 1e-10; // The tolerance on the quality of the solution vector
  
  double diag = 80; // We add this value to the diagonal elements
  // of the matrix to ensure it is positive semi-definite, i.e.
  // all eigenvalues are real and greater than zero. If you change
  // N and the CG fails to converge, you need to increase this
  // value.
  
  bool laplace_mat = false; // If this is true, we use the 2D
  // Laplace operator rather than a random matrix. The 2D laplace
  // op is `sparse` in that most of the entries are zero. We always
  // take advantage of matrix sparsity.
  
  // Starting guess
  std::vector<Complex> x(N, 0.0);
  
  // Vector to solve (for the matrix M, find x = M^{-1} * b)
  std::vector<Complex> b(N, 0.0);
  b[0] = 1.0;

  // The row major matrix
  std::vector<std::vector<Complex>> mat(N, std::vector<Complex> (N, 0.0));  
  
  // Use Eigen to generate a random matrix
  MatrixXcd ref = MatrixXcd::Random(N, N);
  
  // Copy the Eigen matrix into the array
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {

      // If we are using the laplace matrix, construct it
      // else populate the matrix with the random elements
      // from eigen
      mat[i][j] = ref(i,j) + conj(ref(j, i));
      if(i == j) mat[i][j] += diag;
    }
  }
  
  // Pass the row major matrix, guess, source, and constraints to CG
  normalise(b);
  cg(laplace_mat, mat, x, b, tol, maxiter);
  
  // Test for correcteness
  std::vector<Complex> test(N, 0.0);
  matVec(laplace_mat, mat, test, x);

  // The 'test' vector should be the same as b
  cout << "test routine source norm = " << norm(b) << endl;
  axpy(-1.0, test, b);
  cout << "test routine residual norm = " << norm(b) << endl;
  cout << "test routine solution norm = " << norm(x) << endl;  
}

#include <cg.h>

// Basic Linear Algebra Subroutines

// Zero vector
void zero(std::vector<Complex> &x) {
#pragma omp parallel for
  for(int i=0; i<(int)x.size(); i++) x[i] = 0.0;
}

// Copy vector 
void copy(std::vector<Complex> &x, const std::vector<Complex> &y) {
#pragma omp parallel for
  for(int i=0; i<(int)x.size(); i++) x[i] = y[i];
}

// Scale vector 
void ax(double a, std::vector<Complex> &x) {
#pragma omp parallel for
  for(int i=0; i<(int)x.size(); i++) x[i] = a*x[i];
}

// Complex inner product
Complex cDotProd(const std::vector<Complex> &x, const std::vector<Complex> &y) {
  Complex prod = 0.0;
#pragma omp parallel for reduction(+:prod)
  for(int i=0; i<(int)x.size(); i++) prod += conj(x[i]) * y[i];
  return prod;
}

// Inner product
Complex dotProd(const std::vector<Complex> &x, const std::vector<Complex> &y) {
  Complex prod = 0.0;
#pragma omp parallel for reduction(+:prod)
  for(int i=0; i<(int)x.size(); i++) prod += x[i] * y[i];
  return prod;
}

// Norm squared
double norm2(const std::vector<Complex> &x) { 
  double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
  for(int i=0; i<(int)x.size(); i++) sum += (conj(x[i]) * x[i]).real();
  return sum;
}
  
// Norm 
double norm(const std::vector<Complex> &a) { 
  return sqrt(real(norm2(a)));
}

// Normalise the vector
void normalise(std::vector<Complex> &a) {
  double nrm = sqrt(real(norm2(a)));
  ax(1/nrm, a);
}

// axpy
void axpy(const double a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z) {
#pragma omp parallel for
  for(int i=0; i<(int)x.size(); i++) {    
    z[i]  = y[i];
    z[i] += a*x[i];
  }
}

// axpy (in place)
void axpy(const double a, const std::vector<Complex> &x, std::vector<Complex> &y) {
#pragma omp parallel for
  for(int i=0; i<(int)x.size(); i++) {    
    y[i] += a*x[i];
  }
}

void matVec(const bool laplace_mat, const std::vector<std::vector<Complex>> &mat, std::vector<Complex> &result, const std::vector<Complex> &vec) {

  // Sanity checks
  int N = vec.size();
  if(result.size() != vec.size()) {
    cout << "ERROR: Vectors passed to mat vec of unequal size" << endl;
    exit(0);
  }

  if(mat.size() != vec.size()) {
    cout << "ERROR: matrix passed to mat vec of unequal size to vector" << endl;
    exit(0);
  }
  
  // Apply matvec
  for(unsigned int i=0; i<vec.size(); i++) {
    if(laplace_mat) {
      // Take advantage of the sparse matrix pattern
      result[i] = 2.0*vec[i];
      if(i > 0) result[i] += vec[(i-1+N)%N];
      if(i < N-1) result[i] += vec[(i+1)%N];
    }
    else {
      // The matrix is dense so we must do N
      // dot products
      result[i] = dotProd(mat[i], vec);
    }
  }
}

// The Conjugate gradient routine
int cg(const bool laplace_mat, const std::vector<std::vector<Complex>> &mat, std::vector<Complex> &x, const std::vector<Complex> &b, const double tol, const int maxiter) {

  // Temp objects
  std::vector<Complex> temp(x.size(), 0.0);
  std::vector<Complex> res(x.size(), 0.0);
  std::vector<Complex>  p(x.size(), 0.0);
  std::vector<Complex> Ap(x.size(), 0.0);

  bool use_init_guess = false;
  bool verbose = false;
  int success = 0;
  int iter = 0;
  double eps = tol*tol;
  double rsq_new, rsq;
  double alpha, beta;
  double denom;

  // Find norm of rhs.
  double bnorm = norm2(b);
  double bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Error in inverterCG: inverting on zero source or (Nan!" << endl;
    exit(0);
  }

  zero(temp);
  zero(res);
  
  // compute initial residual
  //---------------------------------------  
  if (norm2(x) > 0.0) {        
    // Initial guess supplied: res = b - A*x0
    use_init_guess = true;
    matVec(laplace_mat, mat, temp, x);    
    axpy(-1.0, temp, b, res);
    
    // Update bnorm
    bnorm = norm2(res);
    bsqrt = sqrt(bnorm);

    // temp contains the original guess
    copy(temp, x);
    
    cout << "using initial guess, |x0| = " << norm(temp)
	 << ", |b| = " << bsqrt
	 << ", |res| = " << norm(res) << endl;
    
  } else {
    // No initial guess supplied. Initial residual is the source.    
    copy(res, b);
    zero(temp);
  }
  
  zero(x);
  copy(p, res);  
  rsq = norm2(res);
  //---------------------------------------
  
  // Iterate until convergence
  //---------------------------------------
  for (iter = 0; iter < maxiter; iter++) {
    
    // Compute Ap.
    // This part is usually the bottle neck in any program
    // and we work very hard to optimise it!
    matVec(laplace_mat, mat, Ap, p);

    denom = (cDotProd(p, Ap)).real(); // This is a reduction
    alpha = rsq/denom;
    
    axpy( alpha, p,    x);
    axpy(-alpha, Ap, res);
    
    // Exit if new residual is small enough
    rsq_new = norm2(res); // This is a reduction
    printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);
    if (rsq_new < eps*bnorm) {
      rsq = rsq_new;
      break;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    axpy(beta, p, res, p);
    
  } // End loop over iter
  //---------------------------------------

  if(iter == maxiter) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", iter+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = iter+1; 
  }

  // x contains the solution to the deflated system b - A*x0.
  // We must add back the exact part if using an initial guess too
  axpy(1.0, temp, x);  
  // x now contains the solution to the RHS b.
  
  //Sanity
  cout << "source norm = " << norm(b) << endl;
  cout << "sol norm = " << norm(x) << endl;
  
  matVec(laplace_mat, mat, temp, x);
  axpy(-1.0, temp, b, res);
  double truersq = norm2(res);
  printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter+1, rsq, truersq/(bsqrt*bsqrt));
  
  return success;    
}

void python_test() {
  cout << "Nice! You linked 'void python_test()' to your python script!" << endl;
}

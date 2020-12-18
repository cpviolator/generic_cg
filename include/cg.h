#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <sys/time.h>

using namespace std;

typedef std::complex<double> Complex;

// Zero vector
void zero(std::vector<Complex> &x);

// Copy vector 
void copy(std::vector<Complex> &x, const std::vector<Complex> &y);

// Inner product
Complex cDotProd(const std::vector<Complex> &x, const std::vector<Complex> &y);
  
// Norm squared
double norm2(const std::vector<Complex> &x);

// Norm 
double norm(const std::vector<Complex> &a);

// Normalise
void normalise(std::vector<Complex> &a);

// caxpby
void caxpby(const Complex a, const std::vector<Complex> &x, const Complex b, std::vector<Complex> &y);

// axpby
void axpby(const double a, const std::vector<Complex> &x, const double b, std::vector<Complex> &y);

// axpy in result
void axpby(const double a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z);

// axpy in place
void axpy(const double a, const std::vector<Complex> &x, std::vector<Complex> &y);

// matvec
void matVec(const std::vector<std::vector<Complex>> &mat, std::vector<Complex> &result, const std::vector<Complex> &vec);

// CG op
int cg(const std::vector<std::vector<Complex>> &mat, std::vector<Complex> &x, const std::vector<Complex> &b, const double eps, const int maxiter);

// Simple function to test pythin linkage
void python_test();

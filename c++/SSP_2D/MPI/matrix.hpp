/* Daniel R. Reynolds
   SMU Mathematics
   19 June 2015

   This file defines the Matrix class, as well a variety of linear 
   algebra functions defined based on matrices and vectors.  Here, 
   "vectors" are considered as std::vector<double> objects. */

#ifndef MATRIX_DEFINED__
#define MATRIX_DEFINED__

// Inclusions
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

// for accessing pi via M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#else
#include <cmath>
#endif


//-------------------------------------------------------------

// Arithmetic matrix class
class Matrix {

public:

  // components
  size_t nrows;
  size_t ncols;
  std::vector< std::vector<double> > data;       // matrix data (stored by column)

  // constructors
  Matrix() : nrows(0), ncols(0) {};              // default constructor (empty)
  Matrix(size_t m, size_t n);                    // general constructor (initializes values to 0.0)
  Matrix(size_t m);                              // column-vector constructor (initializes values to 0.0)
  Matrix(size_t m, size_t n, double* vals);      // constructors that copy input data
  Matrix(size_t m, size_t n, std::vector<double> vals);
  Matrix(std::vector<double> vals);
  Matrix(std::vector< std::vector<double> > vals);
  Matrix(std::string mat_spec);                  // string-based matrix constructor
  Matrix(const Matrix& A);                       // copy constructors
  Matrix& operator=(const Matrix& A);

  // dimension accessor routines
  size_t Size() const;
  size_t Columns() const;
  size_t Rows() const;

  // column accessor routines
  std::vector<double>& Column(size_t i);
  std::vector<double>& operator[](size_t i);

  // row accessor (copy) routine
  std::vector<double> Row(size_t i);

  // Matlab/Fortran-style accessors
  double& operator()(size_t i, size_t j);
  double operator()(size_t i, size_t j) const;
  double& operator()(size_t idx);
  double operator()(size_t idx) const;

  // output routines
  int Write() const;                                                      // write to stdout
  int Write(const char *outfile) const;                                   // write to a file
  friend std::ostream& operator<<(std::ostream& os, const Matrix& A);     // streaming output routine

  // in-place operations (C is the matrix calling the routine) -- 0=success, 1=fail
  int Product(const Matrix& A, const Matrix& X);                          // C = A*X
  int LinearSum(double a, const Matrix& A, double b, const Matrix& B);    // C = a*A + b*B
  int Add(const Matrix& A) { return LinearSum(1.0,*this,1.0,A); };        // C = C+A
  int Add(double a);                                                      // C = C+a
  int Subtract(const Matrix& A) { return LinearSum(1.0,*this,-1.0,A); };  // C = C-A
  int Subtract(double a) { return Add(-a); };                             // C = C-a
  int Multiply(const Matrix& A);                                          // C = C.*A
  int Multiply(double a);                                                 // C = a*C
  int Divide(const Matrix& A);                                            // C = C./A
  int Copy(const Matrix& A);                                              // C = A
  int Insert(const Matrix& A, long int is, long int ie,                   // C(is:ie,js:je) = A
	     long int js, long int je);
  int Constant(double a);                                                 // C = a
  int Power(double p);                                                    // C = C.^p
  int Abs();                                                              // Cij = |Cij|
  int Trans();                                                            // C = C^T
  int Inverse();                                                          // C = C^{-1}

  // in-place arithmetic shortcut operators -- 0=success, 1=fail
  int operator+=(const Matrix& A) { return Add(A); };      // C = C+A
  int operator+=(double a)        { return Add(a); };      // C = C+a (add a to all entries)
  int operator-=(const Matrix& A) { return Subtract(A); }; // C = C-A 
  int operator-=(double a)        { return Subtract(a); }; // C = C-a (add -a to all entries)
  int operator*=(const Matrix& A) { return Multiply(A); }; // C = C.*A (componentwise multiply)
  int operator*=(double a)        { return Multiply(a); }; // C = a*C
  int operator/=(const Matrix& A) { return Divide(A); };   // C = C./A (componentwise divide)
  int operator/=(double a)        { return Multiply(1.0/a); };  // C = (1/a)*C
  int operator^=(double p)        { return Power(p); };    // C = C.^p (componentwise power)
  int operator=(double a)         { return Constant(a); }; // C = a (all entries equal a)

  // derived matrix creation operations (C is the output, A calls the operation)
  Matrix T();                                              // C = A^T
  Matrix Extract(long int is, long int ie,                 // C = A(is:ie,js:je)
		 long int js, long int je);

  // Scalar output operators on matrices
  double Min() const;                             // min_ij Cij
  double Max() const;                             // min_ij Cij
  bool operator==(const Matrix& A) const;         // check for Matrix equality

};



//--- supplementary matrix arithmetic routines ---

double Dot(const Matrix& A, const Matrix& B);    // sum_ij (Cij * Aij)
double Norm(const Matrix& A);                    // sqrt(sum_ij Cij^2)
double InfNorm(const Matrix& A);                 // max_i sum_j |Cij|
double OneNorm(const Matrix& A);                 // max_j sum_i |Cij|



//--- new matrix creation routines (C is the output, A and B are the operands) ---

Matrix operator+(const Matrix& A, const Matrix &B);        // C = A+B
Matrix operator-(const Matrix& A, const Matrix &B);        // C = A-B
Matrix operator*(const Matrix& A, const Matrix &B);        // C = A*B
Matrix operator*(const Matrix& A, const double b);         // C = A*b
Matrix operator*(const double a, const Matrix& B);         // C = a*B
Matrix Linspace(double a, double b, size_t m, size_t n);   // linear span
Matrix Logspace(double a, double b, size_t m, size_t n);   // logarithmic span
Matrix Random(size_t m, size_t n);                         // Cij random in [0,1]
Matrix Eye(size_t n);                                      // Cij = delta_i(j)
Matrix MatrixRead(const char *infile);                     // creates from input file
Matrix Inverse(const Matrix& A);                           // C = A^{-1}



//--- supplementary matrix-vector arithmetic routines ---

// standard matrix-vector product -> new vector (function form)
std::vector<double> MatVec(const Matrix& A, const std::vector<double>& v);

// standard matrix-vector product -> new column vector
std::vector<double> operator*(const Matrix& A, const std::vector<double>& v);


//--- supplementary vector<double> vector-arithmetic routines ---

// output routine
int VecWrite(const std::vector<double>&v, const char *outfile);

// inner product between two vectors
double Dot(const std::vector<double>& v1, const std::vector<double>& v2);

// norms
double Norm(const std::vector<double>& v);      // sqrt(sum_i vi^2) (vector 2-norm)
double InfNorm(const std::vector<double>& v);   // max_i |vi|       (vector inf-norm)
double OneNorm(const std::vector<double>& v);   // sum_i |vi|       (vector 1-norm)
double Norm(const std::vector<std::vector<double>>& v);      // sqrt(sum_i vi^2) (vector 2-norm)
double InfNorm(const std::vector<std::vector<double>>& v);   // max_i |vi|       (vector inf-norm)
double OneNorm(const std::vector<std::vector<double>>& v);   // sum_i |vi|       (vector 1-norm)

std::vector<double> abs(std::vector<double>& u);

std::vector<std::vector<double>> abs(std::vector<std::vector<double>> &u);

double max(std::vector<double> &u);

double max(std::vector<std::vector<double>> &u);

// creates a vector of n linearly spaced values from a through b
std::vector<double> Linspace(double a, double b, size_t n);

// creates a vector of n logarithmically spaced values from 10^a through 10^b
std::vector<double> Logspace(double a, double b, size_t n);
	   
// creates a vector of n uniformly-distributed random values
std::vector<double> Random(size_t n);

// output routines 
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);

// extract/insert routines for portions of vectors
std::vector<double> VecExtract(std::vector<double>& x,       // y = x(is:ie)
                               long int is, long int ie);
int VecInsert(std::vector<double>& x, long int is,           // x(is:ie) = y
              long int ie, std::vector<double>& y);

// arithmetic operators for vector<double>
std::vector<double>& operator+=(std::vector<double>& v, const double c);
std::vector<double>& operator+=(std::vector<double>& v, const std::vector<double>& w);
std::vector<double>& operator-=(std::vector<double>& v, const double c);
std::vector<double>& operator-=(std::vector<double>& v, const std::vector<double>& w);
std::vector<double>& operator*=(std::vector<double>& v, const double c);
std::vector<double>& operator*=(std::vector<double>& v, const std::vector<double>& w);
std::vector<double>& operator/=(std::vector<double>& v, const double c);
std::vector<double>& operator/=(std::vector<double>& v, const std::vector<double>& w);
std::vector<double>& operator^=(std::vector<double>& v, const double c);
std::vector<double>& operator^=(std::vector<double>& v, const std::vector<double>& w);

std::vector<double> operator+(const std::vector<double>& v, const double c);
std::vector<double> operator+(const double c, const std::vector<double>& v);
std::vector<double> operator+(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator-(const std::vector<double>& v, const double c);
std::vector<double> operator-(const double c, const std::vector<double>& v);
std::vector<double> operator-(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator*(const std::vector<double>& v, const double c);
std::vector<double> operator*(const double c, const std::vector<double>& v);
std::vector<double> operator*(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator/(const std::vector<double>& v, const double c);
std::vector<double> operator/(const double c, const std::vector<double>& v);
std::vector<double> operator/(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator^(const std::vector<double>& v, const double c);
std::vector<double> operator^(const double c, const std::vector<double>& v);
std::vector<double> operator^(const std::vector<double>& v, const std::vector<double>& w);




//--- linear algebra routines ---

// backward substitution on the linear system U*X = B, filling in an existing Matrix X
//    U and B remain unchanged in this operation; X holds the result
//    B and X may have multiple columns
int BackSubstitution(const Matrix& U, Matrix& X, const Matrix& B);

// backward substitution on the linear system U*X = B, returning X as a new Matrix
//    U and B remain unchanged in this operation
Matrix BackSubstitution(const Matrix& U, const Matrix& B);

// backward substitution on U*x = b, filling in an existing vector<double> x
//    U and b remain unchanged in this operation; x holds the result
int BackSubstitution(const Matrix& U, std::vector<double>& x, const std::vector<double>& b);

// backward substitution on U*x = b, returning x as a new vector<double>
//    U and b remain unchanged in this operation
std::vector<double> BackSubstitution(const Matrix& U, const std::vector<double>& b);

// forward substitution on the linear system L*X = B, filling in the input Matrix X
//    L and B remain unchanged in this operation; X holds the result
//    B and X may have multiple columns
int ForwardSubstitution(const Matrix& L, Matrix& X, const Matrix& B);

// forward substitution on the linear system L*X = B, returning X as a new Matrix
//    L and B remain unchanged in this operation
Matrix ForwardSubstitution(const Matrix& L, const Matrix& B);

// forward substitution on L*x = b, filling in an existing vector<double) x
//    L and b remain unchanged in this operation; x holds the result
int ForwardSubstitution(const Matrix& L, std::vector<double>& x, const std::vector<double>& b);

// forward substitution on L*x = b, returning x as a new vector<double>
//    L and b remain unchanged in this operation
std::vector<double> ForwardSubstitution(const Matrix& L, const std::vector<double>& b);

// solves a linear system A*X = B, filling in the input Matrix X
//    A and B are modified in this operation; X holds the result
int LinearSolve(Matrix& A, Matrix& X, Matrix& B);

// solves a linear system A*X = B, returning X as a new Matrix
//    A and B are modified in this operation; X holds the result
Matrix LinearSolve(Matrix& A, Matrix& B);

// solves a linear system A*x = b, filling in the input vector<double> x
//    A and b are modified in this operation; x holds the result
int LinearSolve(Matrix& A, std::vector<double>& x, std::vector<double>& b);

// solves a linear system A*x = b, returning x as a new vector<double>
//    A and b are modified in this operation; x holds the result
std::vector<double> LinearSolve(Matrix& A, std::vector<double>& b);


#endif  // MATRIX_DEFINED__

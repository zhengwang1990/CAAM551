#include "matrix.hpp"

template<>
void matrix<float>::symeig(matrix <float> &W, matrix <float> &V){
  V = *this;
  W.resize(Nrows,1);

  char JOBZ = 'V';
  char UPLO = 'U';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <float> WORK(LWORK,1);

  int INFO = -999;

  ssyev_ (&JOBZ, &UPLO, &N, V.c_array(), &LDA, W.c_array(), WORK.c_array(), &LWORK, &INFO );
}

template<>
void matrix<double>::symeig(matrix <double> &W, matrix <double> &V){
  V = *this;
  W.resize(Nrows,1);

  char JOBZ = 'V';
  char UPLO = 'U';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <double> WORK(LWORK,1);

  int INFO = -999;

  dsyev_ (&JOBZ, &UPLO, &N, V.c_array(), &LDA, W.c_array(), WORK.c_array(), &LWORK, &INFO );
}

template<>
void matrix<float>::eig(matrix <float> &WR, matrix <float> &WI, matrix <float> &VL, matrix <float> &VR){
  matrix <float> A = *this;
  VL = *this;
  VR = *this;
  WR.resize(Nrows,1);
  WI.resize(Nrows,1);

  char JOBVL = 'V';
  char JOBVR = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <float> WORK(LWORK,1);

  int INFO = -999;

  sgeev_ (&JOBVL, &JOBVR, &N, A.c_array(), &LDA, WR.c_array(), WI.c_array(),
	  VL.c_array(), &LDA, VR.c_array(), &LDA, WORK.c_array(), &LWORK, &INFO);
}

template<>
void matrix<double>::eig(matrix <double> &WR, matrix <double> &WI, matrix <double> &VL, matrix <double> &VR){
  matrix <double> A = *this;
  VL = *this;
  VR = *this;
  WR.resize(Nrows,1);
  WI.resize(Nrows,1);

  char JOBVL = 'V';
  char JOBVR = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <double> WORK(LWORK,1);

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A.c_array(), &LDA, WR.c_array(), WI.c_array(),
	  VL.c_array(), &LDA, VR.c_array(), &LDA, WORK.c_array(), &LWORK, &INFO);
}

// General left matrix inverse not implemented
matrix <double> operator | (const matrix <double> & A, const matrix <double> &B){

  matrix <double> C = B;
  matrix <double> Acopy = A;

  int N    = A.nrows();
  int NRHS = B.ncolumns();
  int LDA  = N;
  int LDB  = B.nrows();
  int INFO;

  matrix <int> IPIV(N,1);

  // Solve A*X = B for X
  dgesv_( &N,
	  &NRHS,
	  Acopy.c_array(),
	  &LDA,
          IPIV.c_array(),
	  C.c_array(),
	  &LDB,
	  &INFO );

  return C;
}

// General left matrix inverse not implemented

matrix <float> operator | (const matrix <float> & A, const matrix <float> &B){

  matrix <float> C = B;
  matrix <float> Acopy = A;

  int N    = A.nrows();
  int NRHS = B.ncolumns();
  int LDA  = N;
  int LDB  = B.nrows();
  int INFO;

  matrix <int> IPIV(N,1);

  // Solve A*X = B for X

  sgesv_( &N,
	  &NRHS,
	  Acopy.c_array(),
	  &LDA,
          IPIV.c_array(),
	  C.c_array(),
	  &LDB,
	  &INFO );
  
  if(INFO)
    cout << "sgesv: INFO=" << INFO << endl;

  return C;
}

template<>
void matrix<double>:: lu (matrix <int> &IPIV){

  int M = Nrows;
  int N = Nrows;
  int LDA = N;
  IPIV.resize(N,1);
  int INFO;

  dgetrf_( &M,
	   &N,
	   data,
	   &LDA,
	   IPIV.c_array(),
	   &INFO);
}

template<>
void matrix<float>:: lu (matrix<int> &IPIV){
  
  int M = Nrows;
  int N = Nrows;
  int LDA = N;
  IPIV.resize(N,1);
  int INFO;
  
  sgetrf_( &M,
	   &N,
	   data,
	   &LDA,
	   IPIV.c_array(),
	   &INFO);

}

template<>
void matrix<double>:: solve(matrix<int> &IPIV, matrix<double> &B){

  char TRANS = 'N';
  int N = Nrows;
  int NRHS = B.ncolumns();
  int LDA = N;
  int LDB = B.nrows();
  int INFO;

  dgetrs_( &TRANS,
	   &N,
	   &NRHS,
	   data,
	   &LDA,
	   IPIV.c_array(),
	   B.c_array(),
	   &LDB,
	   &INFO);
}

template<>
void matrix<float>:: solve(matrix<int> &IPIV, matrix<float> &B){

  char TRANS = 'N';
  int N = Nrows;
  int NRHS = B.ncolumns();
  int LDA = N;
  int LDB = B.nrows();
  int INFO;

  sgetrs_( &TRANS,
	   &N,
	   &NRHS,
	   data,
	   &LDA,
	   IPIV.c_array(),
	   B.c_array(),
	   &LDB,
	   &INFO);
}

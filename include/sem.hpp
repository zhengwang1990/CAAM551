#include "projectHeader.hpp"
#include "matrix.hpp"
#include "setupAide.hpp"
#include "omp.h"

class sem{

private:

  /* basci settings */
  int K, K_c;
  int Ne1, Ne2; // Number of elements
  int N, N_c; // Polynomial order
  int Np1, Np2; 
  datafloat h; // Mesh size
  datafloat J; // Jacobi
  
  /* coordinate */
  fmatrix xLoc, w; // local x and weights
  imatrix xi,yi; // global xi, yi;
  fmatrix r; // global coordinate in 1 demension

  /* nodeId, list and offset */
  imatrix nodeId, list, offset;

  /* interior node numer */
  int Ni;

  /* local matrices */
  fmatrix MLoc, M, MJ; // Mass matrix
  fmatrix S; // Stiffness matrix

  /* right hand side */  
  fmatrix b;

  /* solutions */
  fmatrix uSoln, uExact;

  /* auxiliary array */
  fmatrix aAul;

  /* global stiffness */
  fmatrix globalS;

  /* No. of iteration */
  int iter;
  int iterMax;
  datafloat tol;

  /* additive Schwarz */
  fmatrix *A, *uLoc;
  imatrix *IPIV;
  vector<int> *id;

  /* time */
  double tStart, tEnd;
  double elapsed;

public:
  sem();
  sem(int _K, int _N);
  void init();
  void initSchwarz(int level);
  void gNumbering();
  void assembleMatrices();
  void setupRHS();
  void gather(fmatrix &uAul, fmatrix &u);
  void matVec(fmatrix &u, fmatrix &Au);
  void pcg(int level);
  void cg();
  void setSchwarzId(int level);
  void setSchwarzMatrix(int level);
  void additiveSchwarz(fmatrix &u, fmatrix &Mu, int &start, int &end);
  void postProcessing();
  void matrixGen();
  datafloat partialD(int n, datafloat x);
  datafloat source(datafloat x, datafloat y);

};

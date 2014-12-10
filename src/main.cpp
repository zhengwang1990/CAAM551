#include "sem.hpp"

int main(){

  int K = 10; // No. of elements in one dimension; 2/K is the mesh size h
  int N = 5; // degree of polynomials

  sem semDriver(K,N);

  semDriver.init();

  //semDriver.matrixGen();

  semDriver.pcg(1);

  semDriver.postProcessing();

  semDriver.pcg(2);

  semDriver.postProcessing();

  semDriver.cg();

  semDriver.postProcessing();

  return 0;
}

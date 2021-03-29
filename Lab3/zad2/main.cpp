#include"Functions.h"

int main() {
//   --------metoda sprzezonego gradientu -------------------- 
  double matrixA[N][N] = {0.0};
  double b[N] = {0.0};
  prepareMatrix(matrixA);
  prepareB(b);
  double x[N] = {0.0};
  prepareX(x);
  conjugateGradientMethod(matrixA, b, x);
//   printArray(x);
}


#include"Functions.h"

int main() {
  double matrixA[N][N] = {0.0};
  double b[N] = {0.0};
  prepareMatrix(matrixA);
  prepareB(b);
  double x[N] = {0.0};
  prepareX(x, 1.0);
  methodOfSteepestDescent(matrixA, b, x);

}


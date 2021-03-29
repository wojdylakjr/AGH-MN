#include "MatrixFunctions.h"

int main() {
  double c[N] = { 0,   0, 0, 0.5, 0, -5}; // ustalamy wspolczynniki c w kolejnosci od x^(N-1) do x^0

  double x[N] = {3, 2, 4, -2, -3, 5};
  double x_matrix[N][N];

  prepareMatrix(x, x_matrix);
  double originalMatrix[N][N] = {0};
  prepareMatrix(x, originalMatrix);

  std::cout << "Macierz A: " << std::endl;
  printMatrix(x_matrix);

  //---------------inverse-------------------

  double inverseMatrix[N][N] = {0};
  findInverseMatrix(x_matrix, inverseMatrix);
  double identityMatrix[N][N] = {0};
  matrixMultiply(originalMatrix, inverseMatrix, identityMatrix);
  std::cout << std::endl;
  std::cout << "Macierz jednostkowa uzyskana z mnożenia macierzy A i macierzy "
               "odwrotnej"
            << std::endl;
  printMatrix(identityMatrix);

  //   std::cout <<"Wyznacznik wylicznoej macierzy jednostkowej: ";//
  //   <<std::endl; double det_IdentityMatrix = det(identityMatrix);
  //   printf("%2.25f    ", det_IdentityMatrix);

  //-------------------draw------------------
  drawPlotTXT(c, -5, 5);

  //------------------det--------------------

  double detA = det(originalMatrix);
  std::cout << std::endl << "Wyznacznik macierzy U i macierzy x wynosi: ";
  printf("%2.1f    ", detA);

  //-----wskaznik uwarunkowania macierzy------
  prepareMatrix(x, originalMatrix);
  matrixIndicator(originalMatrix, inverseMatrix);
  return 0;
}
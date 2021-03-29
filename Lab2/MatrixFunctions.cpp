#include "MatrixFunctions.h"

void vectorMultiply(double x[N], double c[N], double result[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i] += pow(x[i], j) * c[j];
    }
  }
}


void printArray(double arr[N]) {
  std::cout << "[ ";
  for (int i = 0; i < N; i++) {
    printf("%2.12f    ", arr[i]);
  }
  std::cout << "]" << std::endl;
}


void prepareMatrix(double x[N], double x_matrix[][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      x_matrix[i][j] = pow(x[i], j);
    }
  }
}


void printMatrix(double matrix[][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%10.4f    ", matrix[i][j]);
    }
    std::cout << std::endl;
  }
}


void GaussMethodLU(double matrix[][N], double L[][N]) {
  for (int columns = 0; columns < N; columns++) {
    L[columns][columns] = 1;
    for (int rows = columns + 1; rows < N; rows++) {
      double l = matrix[rows][columns] / matrix[columns][columns];
      L[rows][columns] = l;
      for (int i = columns; i < N; i++) {
        matrix[rows][i] = matrix[rows][i] - (matrix[columns][i] * l);
      }
    }
  }
}


void solve_U(double U_matrix[][N], double c[N], double z[N]) {
  c[N - 1] = z[N - 1] / U_matrix[N - 1][N - 1];
  for (int rows = N - 2; rows >= 0; rows--) {
    double sum = 0;
    for (int columns = rows + 1; columns < N; columns++) {
      sum += U_matrix[rows][columns] * c[columns];
    }
    c[rows] = (z[rows] - sum) / U_matrix[rows][rows];
  }
}


void solve_L(double L_matrix[][N], double z[N], double y[N]) {
  z[0] = y[0];
  for (int rows = 1; rows < N; rows++) {
    double sum = 0;
    for (int columns = 0; columns < rows; columns++) {
      sum += L_matrix[rows][columns] * z[columns];
    }
    z[rows] = y[rows] - sum;
  }
}


double horner(double wsp[], int stopien, double x) {
  double wynik = wsp[0];

  for (int i = 1; i <= stopien; i++)
    wynik = wynik * x + wsp[i];

  return wynik;
}


double det(double x_matrix[][N]) {
  double L_matrix[N][N] = {0};
  GaussMethodLU(x_matrix, L_matrix);
  double det = 1;
  for (int i = 0; i < N; i++) {
    det *= x_matrix[i][i];
  }
  return det;
}


void drawPlotTXT(double c[N], double minRange, double maxRange) {
  FILE *fp = fopen("lab2_1.txt", "w");

  double x[N] = {0};

  while (minRange < maxRange) {

    for (int i = 0; i < N; i++) {
      if (minRange <= maxRange) {
        x[i] = minRange;
        minRange += 0.1;
      } else{
        return;
      }
    }

    double y[N] = {0};
    double z[N] = {0};
    double x_matrix[N][N];
    double L_matrix[N][N] = {0};
    vectorMultiply(x, c, y);
    prepareMatrix(x, x_matrix);
    GaussMethodLU(x_matrix, L_matrix);
    solve_L(L_matrix, z, y);
    solve_U(x_matrix, c, z);

    for (int i = 0; i < N; i++) {
      fprintf(fp, "%2.12f %2.12f\n", x[i], horner(c, N - 1, x[i]));
    }
  }
}


void findInverseMatrix_L(double L_matrix[][N], double result[][N]) {

  for (int i = 0; i < N; i++) {
    result[i][i] = 1;
    for (int j = 0; j < N; j++) {
      double sum = 0;
      for (int k = j + 1; k < i; k++) {
        sum += L_matrix[i][k] * result[k][j];
      }
      result[i][j] = (L_matrix[i][j] - sum);
    }
  }

  for (int i = 1; i < N; i++) {
    for (int j = 0; j < i; j++) {
      result[i][j] = -result[i][j];
    }
  }
}


void findInverseMatrix_U(double U_matrix[][N], double result[][N]) {
  for (int i = 0; i < N; i++) {
    result[i][i] = 1 / U_matrix[i][i];
    for (int j = i + 1; j < N; j++) {
      double sum = 0;
      for (int k = i; k < j; k++) {
        sum += (result[i][k] * U_matrix[k][j]);
      }
      result[i][j] = -sum / U_matrix[j][j];
    }
  }
  result[N - 1][N - 1] = 1 / U_matrix[N - 1][N - 1];
}


void matrixMultiply(double m1[][N], double m2[][N], double result[][N]) {
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      double sum = 0;
      for (int k = 0; k < N; k++) {
        result[i][j] += m1[i][k] * m2[k][j];
      }
    }
}


void findInverseMatrix(double x_matrix[][N], double result[][N]) {
  double L_matrix[N][N] = {0};
  double L_matrixInverse[N][N] = {0};
  double U_matrixInverse[N][N] = {};

  GaussMethodLU(x_matrix, L_matrix);

  std::cout << std::endl;

  std::cout << "Macierz U" << std::endl;
  printMatrix(x_matrix);
  std::cout << std::endl;
  std::cout << "Macierz L" << std::endl;
  printMatrix(L_matrix);

  // std::cout <<"Macierz odwrtona macierzy U" <<std::endl;
  findInverseMatrix_U(x_matrix, U_matrixInverse);
  //   printMatrix(U_matrixInverse);
  //   std::cout << std::endl;
  //   std::cout <<"Macierz odwrtona macierzy L" <<std::endl;
  findInverseMatrix_L(L_matrix, L_matrixInverse);
  //   printMatrix(L_matrixInverse);

  matrixMultiply(U_matrixInverse, L_matrixInverse, result);
  std::cout << std::endl;
  std::cout << "Macierz odwrtona macierzy A" << std::endl;
  printMatrix(result);
}


double findMaxRow(double matrix[][N]) {
  double max = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (matrix[i][j] > max) {
        max = matrix[i][j];
      }
    }
  }
  return max;
}


void matrixIndicator(double matrix[][N], double inverseMatrix[][N]) {
  double matrix1Max = findMaxRow(matrix);
  double matrix2Max = findMaxRow(inverseMatrix);
  std::cout << std::endl
            << "Norma pierwszej macierzy: " << matrix1Max << std::endl
            << "Norma drugiej macierzy: " << matrix2Max << std::endl
            << "Wskaźnik uwarunkowania macierzy wynosi: "
            << matrix1Max * matrix2Max << std::endl;
}
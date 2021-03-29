#include <cmath>
#include <iostream>
#define N 50
#define L 5
#define X 5 // liczba wartosci wlasnych, ktorych szukamy

void prepareMatrix(double matrix[][N]);
void printMatrix(double matrix[][N]);
void printArray(double arr[N]);
double eigenValues(double matrix[][N], double left, double right, double value, double iteration);
void printPlotTXT(double matrix[][N], double eigenvalue);
double scalarMultiply(double vec1[N], double vec2[N]);
double vectorNorm(double vec1[N], double vec2[N]);
int main() {
  double matrixA[N][N] = {0.0};
  double xVec[N] = {0.0};
  double lambdaVec[X] = {0.0};
  prepareMatrix(matrixA);
  //   printMatrix(matrixA);
  double right = -(matrixA[0][1] + matrixA[2][1]) + matrixA[N - 1][N - 1];
  double left = (matrixA[0][1] + matrixA[2][1]) - matrixA[N - 1][N - 1];

  for (int i = 0; i < X; i++) {
    lambdaVec[i] = eigenValues(matrixA, left, right, i, 50);
    printf("%10.10f  ", lambdaVec[i]);
  }
  printPlotTXT(matrixA, lambdaVec[4]);
}

void prepareMatrix(double matrix[][N]) {
  double delta = 2.0 * L / N;
  for (int i = 1; i < N; i++) {
    matrix[i][i - 1] = -1.0 / (2.0 * pow(delta, 2));
    matrix[i - 1][i] = matrix[i][i - 1];
  }
  for (int i = 0; i < N; i++) {
    double x = -L + (i + 1) * delta;
    matrix[i][i] = pow(delta, -2) + pow(x, 2) / 2.0;
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

void printArray(double arr[N]) {
  std::cout << "[ ";
  for (int i = 0; i < N; i++) {
    printf("%2.4f    ", arr[i]);
  }
  std::cout << "]" << std::endl;
}

double eigenValues(double matrix[][N], double left, double right, double value, double iteration) {
  double lambda = (left + right) / 2;
  if (iteration == 0) {
    return lambda;
  }
  double w[N] = {0};
  int counter = 0;
  w[0] = 1;
  w[1] = matrix[0][0] - lambda;
  if (w[0] * w[1] < 0) {
    counter++;
  }
  for (int i = 2; i < N; i++) {
    w[i] = (matrix[i - 1][i - 1] - lambda) * w[i - 1] -
           pow(matrix[i - 1][i], 2) * w[i - 2];
    if (w[i] * w[i - 1] < 0) {
      counter++;
    }
  }

  if (counter > value) {
    return eigenValues(matrix, left, lambda, value, iteration - 1);
  } else {
    return eigenValues(matrix, lambda, right, value, iteration - 1);
  }
}

void printPlotTXT(double matrix[][N], double eigenvalue) {
  FILE *fp = fopen("lab4_x.txt", "w");
  double xVec[N] = {0};
  xVec[0] = 1;
  xVec[1] = (eigenvalue - matrix[0][0]) / matrix[0][1];
  for (int i = 2; i < N; i++) {
    xVec[i] = (((eigenvalue - matrix[i - 1][i - 1]) * xVec[i - 1] -
                matrix[i - 2][i - 1] * xVec[i - 2]) /
               matrix[i - 1][i]);
  }
  double norm = vectorNorm(xVec, xVec);

  for (int i = 0; i < N; i++) {
    // fprintf(fp, "%d  %f \n", i,xVec[i]/norm);
    fprintf(fp, "%f  \n", xVec[i] / norm);
  }
}

double scalarMultiply(double vec1[N], double vec2[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    sum += vec1[i] * vec2[i];
  }
  return sum;
}

double vectorNorm(double vec1[N], double vec2[N]) {
  return sqrt(scalarMultiply(vec1, vec2));
}
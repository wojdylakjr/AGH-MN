#include "Functions.h"
void prepareMatrix(double matrix[][N]) {
  double m = 5.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (abs(i - j) <= m) {
        matrix[i][j] = static_cast<double>(1.0 / (1.0 + abs(i - j)));
      } else {
        matrix[i][j] = 0.0;
      }
    }
  }
}

void prepareB(double arr[N]) {
  for (int i = 0; i < N; i++) {
    arr[i] = static_cast<double>(i);
  }
}
void prepareX(double arr[N], double number) {
  for (int i = 0; i < N; i++) {
    arr[i] = number;
  }
}

void printArray(double arr[N]) {
  std::cout << "[ ";
  for (int i = 0; i < N; i++) {
    printf("%2.4f    ", arr[i]);
  }
  std::cout << "]" << std::endl;
}

void printMatrix(double matrix[][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%10.4f    ", matrix[i][j]);
    }
    std::cout << std::endl;
  }
}

void methodOfSteepestDescent(double A[][N], double b[N], double x[N]) {
    double r[N] = {0.0};
  FILE *fp = fopen("lab3_x.txt", "w");
  int counter = 0;
  do {
    double multiplied[N] = {0};
    matrixMultiplyWithVector(A, x, multiplied); // A * x
    vectorSubtracktion(b, multiplied, r); // b - A * x 
    double temp[N] = {0};
    matrixMultiplyWithVector(A, r, temp); // A * r 
    double alfa = scalarMultiply(r, r) / scalarMultiply(r, temp); // r*r / r*Ar
    double temp2[N] = {0};
    vectorMultiplyWithConst(alfa, r, temp2); // alfa * r

    fprintf(fp, "i: %d |r| =%5.10f alfa =%5.10f |x| =%5.5f\n", counter,
            vectorNorm(r, r), alfa, vectorNorm(x, x));

    vectorAdding(x, temp2); // x = x + alfa * r
    counter++;
  } while (vectorNorm(r, r) > 0.00001);
}

void matrixMultiplyWithVector(double matrix1[][N], double x[N],
                              double result[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i] += matrix1[i][j] * x[j];
    }
  }
}

void vectorSubtracktion(double vec1[N], double vec2[N], double result[N]) { // odejmowanie
  for (int i = 0; i < N; i++) {
    result[i] = vec1[i] - vec2[i];
  }
}

double scalarMultiply(double vec1[N], double vec2[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    sum += vec1[i] * vec2[i];
  }
  return sum;
}

void vectorMultiplyWithConst(double alfa, double vec1[N], double result[N]) {
  for (int i = 0; i < N; i++) {
    result[i] = alfa * vec1[i];
  }
}

void vectorAdding(double vec1[N], double vec2[N]) {
  for (int i = 0; i < N; i++) {
    vec1[i] += vec2[i];
  }
}

double vectorNorm(double vec1[N], double vec2[N]) {
  return sqrt(scalarMultiply(vec1, vec2));
}


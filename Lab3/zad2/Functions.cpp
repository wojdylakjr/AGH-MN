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
    arr[i] = static_cast<double>(i + 1);
  }
}
void prepareX(double arr[N]) {
  for (int i = 0; i < N; i++) {
    arr[i] = 0;
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

void conjugateGradientMethod(double A[][N], double b[N], double x[N]) {
    double r[N] = {0.0};
    double v[N] = {0.0};
    double temp[N] = {0.0};
    matrixMultiplyWithVector(A, x, temp); //A*x 
    vectorSubtracktion(b, temp, r); //r = b - A*x
    vectorCopy(v,r); // v = r
  FILE *fp = fopen("lab3_x.txt", "w");
  int counter = 0;
 while (vectorNorm(r, r) > 0.00001){

    double temp1[N] = {0.0};
     matrixMultiplyWithVector(A, v, temp1); // A * v = temp1
    double alfa = scalarMultiply(r, r) / scalarMultiply(v, temp1); // r*r / v*Av
    double temp2[N] = {0.0};
    vectorMultiplyWithConst(alfa, v, temp2); // alfa * v = temp2
    double xNorm = vectorNorm(x, x);
    vectorAdding(x, temp2); // x = x + alfa * v
    double temp3[N] = {0.0};
    vectorMultiplyWithConst(alfa, temp1, temp3); //alfa * A*v = alfa * temp1 = temp3
    double r_copy[N] = {0.0};
    vectorCopy(r_copy, r);
    vectorSubtracktion(r, temp3); // r =r -alfa *A*v 
    double beta = scalarMultiply(r, r) / scalarMultiply(r_copy, r_copy);
    double temp4[N] = {0.0};
    vectorMultiplyWithConst(beta, v, temp4); //beta * v = temp4
    vectorAdding(r, temp4, v); //v = r + beta*vec2

    fprintf(fp, "i: %d |r| =%5.10f alfa =%5.10f beta =%5.10f |x| =%5.10f\n", counter,
            vectorNorm(r_copy, r_copy), alfa,beta, xNorm);
    counter++;
   }
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

void vectorSubtracktion(double vec1[N], double vec2[N]) { // odejmowanie
  for (int i = 0; i < N; i++) {
     vec1[i] -= vec2[i];
  }
}
void vectorSubtracktion(double vec1[N], double vec2[N], double result[N]){
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

void vectorAdding(double vec1[N], double vec2[N], double result[N]){
    for (int i = 0; i < N; i++) {
    result[i] = vec1[i] + vec2[i];
  }
}

double vectorNorm(double vec1[N], double vec2[N]) {
  return sqrt(scalarMultiply(vec1, vec2));
}


void vectorCopy(double vec[N] ,double copiedVec [N]){
    for(int i = 0; i < N; i++){
        vec[i] = copiedVec[i];
    }
}
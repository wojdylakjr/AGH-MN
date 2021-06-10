#include <cmath>
#include <iostream>
#define N 21

double function1(double x);
void printMatrix(double matrix[][N]);
void printArray(double arr[N]);
void setX(double vecX[N], double delta, double min);
void setY(double vecY[N], double vecX[N]);
void wyzM(double m[N], double vecX[N], double vecY[N], double alfa, double beta,
          double min, double max);
void solveEquation(double matrix[][N], double b[N], double[N]);
void matrixMultiply(double matrix1[][N], double x[N], double result[N]);
void solve(double matrix[][N], double b[N], double x[N]);
double wyzS(double m[N], double vecX[N], double vecY[N], double x, double min,
            double max);
double deriverate(double x, double dx);

int main() {
  double alfa = 0.0, beta = 0.0;
  double min = -5.0, max = 5.0;
  double vecX[N] = {0.0};
  double vecY[N] = {0.0};
  double vecM[N] = {0.0};

  double delta = (max - min) / ((double)N - 1.0);
  setX(vecX, delta, min);
  setY(vecY, vecX);
  printArray(vecX);
  printArray(vecY);
  wyzM(vecM, vecX, vecY, alfa, beta, min, max);
  printArray(vecM);

  FILE *f1;
  f1 = fopen("s(x).txt", "w");

  FILE *f2;
  f2 = fopen("wezly.txt", "w");
    FILE *f3;
  f3 = fopen("funkcja.txt", "w");
      FILE *f4;
  f4 = fopen("pochodne.txt", "w");
  for(int i = 0; i < N; i++){
      fprintf(f2, "%f %f\n", vecX[i], vecM[i]);
  }
  for (double i = min; i <= max; i += 0.05) {
    double temp = wyzS(vecM, vecX, vecY, i, min, max);
    fprintf(f1, "%f %f\n", i, temp);
  }
for (double i = min; i <= max; i += 0.05) {
    fprintf(f3, "%f %f\n", i, function1(i));
  }

  for(int i = 0; i < N;i++){
      fprintf(f4, "%f %f\n", vecX[i], deriverate(vecX[i], 0.01));
  }
  
  fclose(f1);
  fclose(f2);
  fclose(f3);

  return 0;
}

void printMatrix(double matrix[][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%10.6f    ", matrix[i][j]);
      //   printf("%e    ", matrix[i][j]);
    }
    std::cout << std::endl;
  }
}

void printArray(double arr[N]) {
  std::cout << "[ ";
  for (int i = 0; i < N; i++) {
    printf("%2.4f    ", arr[i]);
    // printf("%e    ", arr[i]);
  }
  std::cout << "]" << std::endl;
}

double function1(double x) {
//   return 1.0 / (1.0 + pow(x, 2));
     return cos(2*x);
}

void setX(double vecX[N], double delta, double min) {
  for (int i = 0; i < N; i++) {
    vecX[i] = min + delta * i;
  }
}

void setY(double vecY[N], double vecX[N]) {
  for (int i = 0; i < N; i++) {
    vecY[i] = function1(vecX[i]);
  }
}

void wyzM(double m[N], double vecX[N], double vecY[N], double alfa, double beta,
          double min, double max) {
  double matrixA[N][N] = {0.0};
  double vecD[N] = {0.0};

  double d = (max - min) / ((double)N - 1.0);
  double lambda = 0.5; //=1/2
  double mi = 1.0 - lambda;
  vecD[0] = alfa;
  vecD[N - 1] = beta;
  // wartosci wektora
  for (int i = 1; i < N - 1; ++i) {
    vecD[i] = (6.0 / (d + d)) *
              ((vecY[i + 1] - vecY[i]) / d - (vecY[i] - vecY[i - 1]) / d);
  }

  // diagonala macierzy A
  for (int i = 0; i < N; ++i) {
    matrixA[i][i] = 2.0;
  }
  matrixA[0][0] = 1.0;
  matrixA[N-1][N-1] = 1.0;
  // pozostale wartosci
  for (int i = 1; i < N - 1; i++) {
    for (int j = 1; j < N - 1; j++) {
      if (i == j) {
        matrixA[i][j - 1] = mi;
        matrixA[i][j + 1] = lambda;
      }
    }
  }

  // rozwiazanie macierzy
  solve(matrixA, vecD, m);
}
void GaussMethod(double matrix[][N], double b[N]) {
  for (int columns = 0; columns < N; columns++) {
    for (int rows = columns + 1; rows < N; rows++) {
      double l = matrix[rows][columns] / matrix[columns][columns];
      for (int i = columns; i < N; i++) { // i iteruje po kolumnach
        matrix[rows][i] = matrix[rows][i] - (matrix[columns][i] * l);
      }
      b[rows] = b[rows] - (b[columns] * l);
    }
  }
}

void solveEquation(double matrix[][N], double b[N], double x[N]) {
  x[N - 1] = b[N - 1] / matrix[N - 1][N - 1];
  for (int rows = N - 2; rows >= 0; rows--) {
    double sum = 0;
    for (int columns = rows + 1; columns < N; columns++) {
      sum += matrix[rows][columns] * x[columns];
    }
    x[rows] = (b[rows] - sum) / matrix[rows][rows];
  }
}

void matrixMultiply(double matrix1[][N], double x[N], double result[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i] += matrix1[i][j] * x[j];
    }
  }
}

void solve(double matrix[][N], double b[N], double x[N]) {
  GaussMethod(matrix, b);
  solveEquation(matrix, b, x);
}

double wyzS(double m[N], double vecX[N], double vecY[N], double x, double min,
            double max) {
  int przedzial = 0;
  double Sx = 0;
  double hi = (max - min) / ((double)N - 1.0);

  for (int i = 1; i < N; i++) {
    if (vecX[i - 1] <= x && x <= vecX[i]) {
      przedzial = i - 1;
      break;
    }
  }

  double Ai;
  double Bi;

  double tmp1 = m[przedzial + 1];
  double tmp2 = m[przedzial];
  Ai = ((vecY[przedzial + 1] - vecY[przedzial]) / hi) -
       (hi / 6.0) * (tmp1 - tmp2);

  tmp1 = m[przedzial];
  Bi = vecY[przedzial] - tmp1 * ((pow(hi, 2)) / 6.0);

  // Wyznaczanie wzoru (8), stopniowo, zeby sie nie pogubic
  Sx = m[przedzial];
  Sx *= (pow((vecX[przedzial + 1] - x), 3) / (6.0 * hi));
  Sx += m[przedzial + 1] * (pow((x - vecX[przedzial]), 3) / (6.0 * hi));
  //    Sx *= pow((x - vecX[przedzial]), 3) / 6.0*hi;
  Sx += Ai * (x - vecX[przedzial]);
  //    Sx *= (x - vecX[przedzial]);
  Sx += Bi;

  return Sx;
}

double deriverate(double x, double dx) {
  return (function1(x - dx) - 2 * function1(x) + function1(x + dx)) /
         pow(dx, 2);
}
#include <cmath>
#include <iostream>
const int N = 5;
void printMatrix(double matrix[][N], double b[N]);
void printArray(double arr[N]);
void GaussMethod(double matrix[][N], double b[N]);
void solveEquation(double matrix[][N], double b[N], double[N]);
void matrixMultiply(double matrix1[][N], double x[N], double result[N]);
double standardDaviation(double array1[N], double array2[N]);
void solve(double matrix[][N], double b[N], double x[N]);

int main() {
    double backupB[N] = {10, 2, 9, 8, 3};


  FILE *fp = fopen("lab1.2_x.txt", "w");

  for (double q = 0.2; q <= 5; q += 0.2001) {
    double b[N] = {10, 2, 9, 8, 3};
    double c[N] = {0, 0, 0, 0, 0};
    double x[N] = {1, 1, 1, 1, 1};

    double originalMatrix[][N] = {{q * 0.0002, 1, 6, 9, 10},
                                  {0.0002, 1, 6, 9, 10},
                                  {1, 6, 6, 8, 6},
                                  {5, 9, 10, 7, 10},
                                  {3, 4, 9, 7, 9}};
    double matrix[][N] = {{q * 0.0002, 1, 6, 9, 10},
                          {0.0002, 1, 6, 9, 10},
                          {1, 6, 6, 8, 6},
                          {5, 9, 10, 7, 10},
                          {3, 4, 9, 7, 9}};

    solve(matrix, b, x);
    std::cout << "Dla q = " << q << ": ";
    matrixMultiply(originalMatrix, x, c);

    fprintf(fp, "%2.12f %2.12f\n", q, standardDaviation(c, backupB));
    std::cout << "Szacowane b: ";
    printArray(c);
    std::cout << "\n \n";
  }

}

void printMatrix(double matrix[][N], double b[N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << "  " << matrix[i][j] << "  ";
    }
    std::cout << "  | " << b[i] << std::endl;
  }
}

void printArray(double arr[N]) {
  std::cout << "[ ";
  for (int i = 0; i < N; i++) {
    std::cout << arr[i] << ",   ";
    // printf("%2.12f    ",arr[i] );
  }
  std::cout << "]" << std::endl;
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

double standardDaviation(double array1[N], double array2[N]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    sum = sum + pow(array1[i] - array2[i], 2);
  }
  return (sqrt(sum)/5.0);
}

void solve(double matrix[][N], double b[N], double x[N]){
       GaussMethod(matrix, b);
    solveEquation(matrix, b, x);
}

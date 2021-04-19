#include <cmath>
#include <iostream>
#define N 2

void prepareMatrix(double matrix[][N], double vecXY[N]);
void prepareVector(double vector[N], double vecXY[N]);
void vectorSubstraction(double vector1[N], double vector2[N]);
void printMatrix(double matrix[][N]);
void printArray(double arr[N]);
double vectorNorm(double vec1[N], double vec2[N]);
void inverse(double A[N][N], double inverse[N][N]);
void adjoint(double A[N][N], double adj[N][N]);
double determinant(double A[N][N], double n);
void getCofactor(double A[N][N], double temp[N][N], double p, double q, double n);
void matrixMultiplyWithVector(double matrix1[][N], double x[N], double result[N]);




int main() {
  double matrix[N][N] = {0.0};
  double vecV[N] = {0.0};
  double vecR[N] = {10, 10};

  prepareMatrix(matrix, vecR);
  prepareVector(vecV, vecR);
//   std::cout << "Macierz wejściowa: "<<std::endl;
//   printMatrix(matrix);
//   std::cout << "Wektor wejściowy: "<<std::endl;
//   printArray(vecV);
//   printArray(vecR);
  double inversedMatrix[N][N] = {0.0};
  double vectorDelta[N]{1.0};
  inverse(matrix, inversedMatrix);
  matrixMultiplyWithVector(inversedMatrix, vecV, vectorDelta);
  double tempVector[N] = {1.0};// posluzy nam do przechowania wektora równego 𝐫𝒌−𝒓𝒌−𝟏


int i = 0;

  FILE *fp = fopen("lab6.txt", "w");
  
  double x,y = 0;
  while(vectorNorm(tempVector, tempVector)  > pow(10,-6)){

      std::cout<<"i: "<<i<<" x: " << vecR[0] <<" y: " << vecR[1]<<" delta_r: "<<vectorNorm(tempVector, tempVector)  <<std::endl;
          fprintf(fp, " %d   %5.6f  %5.6f  %5.6f\n", i,
            vecR[0] , vecR[1],vectorNorm(tempVector, tempVector));
      x = vecR[0]; y = vecR[1];
      prepareMatrix(matrix, vecR);
      prepareVector(vecV, vecR);
      inverse(matrix, inversedMatrix);
      matrixMultiplyWithVector(inversedMatrix, vecV, vectorDelta);
      vecR[0] -=  vectorDelta[0];
      vecR[1] -= vectorDelta[1];
      tempVector[0] = vecR[0] - x;
      tempVector[1] = vecR[1] - y;
      ++i;

  }

}



void prepareMatrix(double matrix[][N], double vecXY[N]) {
    double x = vecXY[0];
    double y = vecXY[1];

    matrix[0][0] = 2 * pow(y, 2) - 6*x*y;
    matrix[1][0] = 2 *x* pow(y, 3) + 2*y;
    matrix[0][1] = 4*x*y - 3*pow(x, 2);
    matrix[1][1] = 3 * pow(y, 2) * pow(x,2) + 2*x;
}

void prepareVector(double vector[N], double vecXY[N]) {
     double x = vecXY[0];
    double y = vecXY[1];
    vector[0] = 2*x*pow(y,2) - 3*pow(x, 2)*y -2;
    vector[1] = pow(x,2)*pow(y,3) + 2*x*y - 12;
}

void matrixMultiplyWithVector(double matrix1[][N], double x[N], double result[N]) {
    result[0] = 0, result[1] =0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i] += matrix1[i][j] * x[j];
    }
  }
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



void vectorMultiplyWithConst(double alfa, double vec1[N], double result[N]) {
  for (int i = 0; i < N; i++) {
    result[i] = alfa * vec1[i];
  }
}

void vectorSubtracktion(double vec1[N], double vec2[N], double result[N]) {
  for (int i = 0; i < N; i++) {
    result[i] = vec1[i] - vec2[i];
  }
}

void vectorMultiply(double vec1[N], double vec2[N], double result[N][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i][j] = vec1[i] * vec2[j];
    }
  }
}

void matrixMultiplyWithConst(double alfa, double matrix[][N]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matrix[i][j] *= alfa;
    }
  }
}



void matrixMultiply(double mat1[][N], double mat2[][N], double res[][N]) {
  int i, j, k;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      // res[i][j] = 0;
      for (k = 0; k < N; k++)
        res[i][j] += mat1[i][k] * mat2[k][j];
    }
  }
}




void getCofactor(double A[N][N], double temp[N][N], double p, double q,
                 double n) {
  int i = 0, j = 0;

  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {

      if (row != p && col != q) {
        temp[i][j++] = A[row][col];

        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

double determinant(double A[N][N], double n) {
  double D = 0; 

  if (n == 1)
    return A[0][0];

  double temp[N][N];
  int sign = 1;

  for (int f = 0; f < n; f++) {
    getCofactor(A, temp, 0, f, n);
    D += sign * A[0][f] * determinant(temp, n - 1);
    sign = -sign;
  }

  return D;
}

void adjoint(double A[N][N], double adj[N][N]) {
  if (N == 1) {
    adj[0][0] = 1;
    return;
  }

  int sign = 1;
  double temp[N][N];

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      getCofactor(A, temp, i, j, N);
      sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj[j][i] = (sign) * (determinant(temp, N - 1));
    }
  }
}

void inverse(double A[N][N], double inverse[N][N]) {

  double det = determinant(A, N);
  if (det == 0) {
    std::cout << "Singular matrix, can't find its inverse";
  }

  double adj[N][N];
  adjoint(A, adj);

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      inverse[i][j] = adj[i][j] / double(det);
}


void matrixCopy(double matrix[N][N], double copy[N][N]){
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      copy[i][j] = matrix[i][j];
    }
  }
}

#include <iostream>
#include<cmath>
#define N 5

double function(double x);
void setMatrix(double matrix[][N + 1], double a, double dx);
double wielomian(double wn[N + 1],double matrix[][N + 1],  double x, double min, double max);
void printMatrix(double matrix[][N + 1]);
void printArray(double arr[N + 1]);
double czebyszew(double min, double max, int m, int n);

int main() {
  double min = -5.0, max = 5.0;
  double matrixF[N + 1][N + 1]{0.0};
  double dx = (max - min)/ N;
  setMatrix(matrixF, min, dx);

  double Wn[N + 1]{0.0};
  for(double x = min; x <= max; x+= 0.02){
      std::cout<<x<< " "<<wielomian(Wn, matrixF, x,  min, max)<<std::endl;
  }

  

}

double function(double x){
    return 1/(1 + x*x);
}

void setMatrix(double matrixF[][N + 1], double a, double dx){
    for(int i = 0; i <= N; i++){
        matrixF[0][i] = function(i);
    }

    for(int i = 1; i <= N; i++){
        for(int j = i; j <= N; j++){
            matrixF[i][j] = (matrixF[i][j -1] - matrixF[i - 1][j - 1]) / ((a + i * dx) - (a + (i - j)*dx));
        }
    }
}

void printMatrix(double matrix[][N + 1]) {
  for (int i = 0; i <= N; i++) {
    for (int j = 0; j <= N; j++) {
      printf("%10.6f    ", matrix[i][j]);
      //   printf("%e    ", matrix[i][j]);
    }
    std::cout << std::endl;
  }
}

void printArray(double arr[N + 1]) {
  std::cout << "[ ";
  for (int i = 0; i <=N; i++) {
    printf("%2.4f    ", arr[i]);
    // printf("%e    ", arr[i]);
  }
  std::cout << "]" << std::endl;
}

double czebyszew(double min, double max, int m, int n)
{
    return 0.5 * ((max - min) * cos(M_PI * ((2.0 * (double)(m) + 1.0) / (2.0 * (double)(n) + 2.0))) + (min + max));
}

double wielomian(double wn[N + 1],double matrix[][N + 1],  double x, double min, double max){
    double sum = 0;
    for(int j = 0; j <= N; j++){
        double iloczyn = 1;
        for(int i = 0; i <= j -1; i++){
            iloczyn *= (x-(min+i*i));
        }
        sum += matrix[j][j] + iloczyn;
    }
    return sum;
}
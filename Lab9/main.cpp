#include <cmath>
#include <iostream>
#define N 10

double function(double x);
double F1(double x1, double x2);
double F2(double x1, double x2, double x3);

int main() {
  double min = -1.5;
  double max = 1.0;
  double h = 0.01;
  double x1 = -0.5;
  double x2 = x1 + h;
  double x3 = x2 + h;
  double x_m = 0;
  double eps = 1e-8;
  double minimum = 0;
  FILE *dane;
  dane = fopen("dane.txt", "w");
  fprintf(dane, "Nr.      x1         x2           x3          x_m          F1   "
                "     F2\n");

  for (int i = 0; i < N; i++) {
    x_m = ((x1 + x2) / 2.0) - (F1(x1, x2) / (2 * F2(x1, x2, x3)));
    fprintf(dane, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i + 1, x1, x2, x3, x_m,
            F1(x1, x2), F2(x1, x2, x3));
    if (fabs(x_m - x1) < eps || fabs(x_m - x2) < eps || fabs(x_m - x3) < eps) {
      minimum = x_m;
      break;
    }
    if (fabs(x_m - x1) < fabs(x_m - x2)) {
      if (fabs(x_m - x2) < fabs(x_m - x3)) {
        x3 = x_m;
      } else {
        x2 = x_m;
      }

    } else {
      if (fabs(x_m - x1) < fabs(x_m - x3)) {
        x3 = x_m;
      } else {
        x1 = x_m;
      }
    }
  }

  printf("Minimum: %.10f \n", minimum);
  fclose(dane);
  return 0;
}

double function(double x) { 
  return log(pow(x, 5) + 3 * pow(x, 2) + x + 9); 
//   return pow(x,6);
}

double F1(double x1, double x2) {
  return (function(x2) - function(x1)) / (x2 - x1);
}

double F2(double x1, double x2, double x3) {
  return (F1(x2, x3) - F1(x1, x2)) / (x3 - x1);
}
#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#define N 6 // wielomian jest stopnina N-1

void vectorMultiply(double x[N], double c[N], double result[N]); 
void matrixMultiply(double m1[][N], double m2[][N], double result[][N]);
void printArray(double arr[N]);
void prepareMatrix(double x[N], double x_matrix[][N]); //przygotowuje macierz A na podstawie wektora x-ów
void printMatrix(double matrix[][N]);
void GaussMethodLU(double matrix[][N], double L[][N]); //tworzy macierz L oraz zamienia macierz x-ów na macierz U
void solve_U(double U_matrix[][N], double result[N], double[N]); //obliczamy wektor pomocniczy z w rownaniu L*z = y
void solve_L(double L_matrix[][N], double result[N], double y[N]); //obliczamy wektor c współczynikow z równania U*c = z
double horner(double wsp[], int stopien, double x); //zwraca wartosc wielomianu dla konkretego x 
double det(double U_matrix[][N]); //zwraca wyznacznik macierzy U
void drawPlotTXT(double c[N], double minRange, double maxRange); //funkcja zajmujaca sie wygenerowaniem danych do wykresu
void findInverseMatrix_L(double L_matrix[][N], double result[][N]);
void findInverseMatrix_U(double U_matrix[][N], double result[][N]);
void findInverseMatrix(double x_matrix [][N], double result[][N]);
double findMaxRow(double matrix [][N]); //zwraca najwieksza suma skladnikow w wierszu w  macierzy
void matrixIndicator(double matrix[][N], double inverseMatrix[][N]);
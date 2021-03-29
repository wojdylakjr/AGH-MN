#pragma once
#include <cmath>
#include <iostream>
#define N 1000

void prepareMatrix(double matrix[][N]);
void prepareB(double arr[N]);
void prepareX(double arr[N]);
void printArray(double arr[N]);
void prepareMatrix(double x[N], double x_matrix[][N]);
void printMatrix(double matrix[][N]);
//void methodOfSteepestDescent(double A[][N], double b[N], double x[N]);
void matrixMultiplyWithVector(double matrix1[][N], double x[N], double result[N]);
void vectorSubtracktion(double vec1[N], double vec2[N]); //wynik odejmowania jest zapisany w wektorze vec1
void vectorSubtracktion(double vec1[N], double vec2[N], double result[N]); //wynik odejmowania jest zapisany w wektorze result
double scalarMultiply(double vec1[N], double vec2[N]);
void vectorMultiplyWithConst(double alfa, double vec1[N], double result[N]);
void vectorAdding(double vec1[N], double vec2[N]); // wynik dodawania jest zapisany w wektorze vec1
void vectorAdding(double vec1[N], double vec2[N], double result[N]); // wynik dodawania jest zapisany w wektorze vec1
double vectorNorm(double vec1[N], double vec2[N]);

void conjugateGradientMethod(double A[][N], double b[N], double x[N]);
void vectorCopy(double vec[N] ,double copiedVec [N]);


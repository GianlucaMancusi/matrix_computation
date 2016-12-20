// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi - Unimore
#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
extern double* compMinor(double *matr, int matr_rows, int matr_cols, int row, int col);
extern double det3x3(double *matr);
extern double laplace(double *matr, int dim);
#endif //!MATRICE_H
// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi - Unimore
#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
extern double det(double *matr, int dim);
extern void fillCompMinor(const double *src_matr, double *dst_matr, int src_dim, int row, int col);
extern double* compMinor(double *matr, int matr_rows, int matr_cols, int row, int col);
extern double det3x3(double *matr);
extern double laplace(double *matr, int dim, double **compm_matrs, int start_dim);
#endif //!MATRICE_H
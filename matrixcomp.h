// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi - Unimore
#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
#include <stdarg.h>

struct matrix {
	size_t rows;
	size_t cols;
	double* data;
};

//====== Useful functions =======

extern double det(struct matrix *matr);//Calculates the determinant of matrices
extern struct matrix *mulmatr(struct matrix *lhs, struct matrix *rhs);//Multiply 2 matrices
extern struct matrix *creatematrix(size_t rows, size_t cols);
extern struct matrix *createemptymatrix(size_t rows, size_t cols); //empty matrix;
extern int matrixrow(struct matrix* matr, size_t rowid, double *row);
extern void destroymatrix(struct matrix* m);

//====== Internal functions =======

extern void fillCompMinor(const double *src_matr, double *dst_matr, int src_dim, int row, int col);
extern double* compMinor(double *matr, int matr_rows, int matr_cols, int row, int col);
extern double det3x3(double *matr);
extern double laplace(double *matr, int dim, double **compm_matrs, int start_dim);
#endif //!MATRIX_H
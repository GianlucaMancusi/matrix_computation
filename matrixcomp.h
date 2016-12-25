// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi -  from the"Enzo Ferrari" Department of Engineering
// http://www.ingmo.unimore.it/site/en/home.html
#ifndef MATRIXCOMP_H
#define MATRIXCOMP_H
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>

struct matrix {
	size_t rows;
	size_t cols;
	double* data;
};

struct matrix_selection
{
	size_t rows;
	size_t cols;
	uint8_t selection;
};

enum
{
	row,
	col
};

//====== Useful functions =======

extern double det(struct matrix *matr);//Calculates the determinant of matrices
extern struct matrix *mulmatr(struct matrix *lhs, struct matrix *rhs);//Multiply 2 matrices
extern struct matrix *creatematr(size_t rows, size_t cols);
extern struct matrix *createemptymatr(size_t rows, size_t cols); //empty matrix;
extern struct matrix *clonematr(struct matrix* matr);
extern struct matrix* matrcompminor(struct matrix *matr, int row, int col); //Find a complementary minor from row and col
extern int matrrow(struct matrix* matr, size_t rowid, double *row); //Change the entire row
extern void destroymatr(struct matrix* matr);

//====== Internal functions =======

extern void fillCompMinor(const double *src_matr, double *dst_matr, size_t src_dim, size_t row, size_t col);
extern double det3x3(double *matr);
extern double laplace(double *matr, size_t dim, double **compm_matrs, size_t start_dim);
extern struct matrix_selection findlinewithmorezeros(struct matrix *matr);
#endif //!MATRIXCOMP_H
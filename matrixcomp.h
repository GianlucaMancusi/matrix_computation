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

enum linetype
{
	row,
	col
};

//====== Useful functions =======

extern double det(const struct matrix *matr);//Calculates the determinant of matrices
extern struct matrix *mulmatr(const struct matrix *lhs, const struct matrix *rhs);//Multiply 2 matrices
extern struct matrix *creatematr(size_t rows, size_t cols);
extern struct matrix *createemptymatr(size_t rows, size_t cols); //empty matrix;
extern struct matrix *creatematrfrom(const double* matr, size_t row, size_t col); //create matrix from 3 parameters
extern struct matrix *clonematr(const struct matrix* matr);
extern struct matrix* matrcompminor(const struct matrix *matr, int row, int col); //Find a complementary minor from row and col
extern int matrrow(const struct matrix* matr, size_t rowid, const double *row); //Change the entire row
extern void destroymatr(struct matrix* matr);

//====== Internal functions =======

extern void fillCompMinor(const double *src_matr, double *dst_matr, size_t src_dim, size_t row, size_t col);
extern double det3x3(const double *matr);
extern double laplace(const double *matr, size_t dim, double **compm_matrs, size_t start_dim, const struct matrix_selection *selection);
extern struct matrix_selection findlinewithmorezeros(const struct matrix *matr);
#endif //!MATRIXCOMP_H
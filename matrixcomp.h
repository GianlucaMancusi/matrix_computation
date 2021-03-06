// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi -  from the"Enzo Ferrari" Department of Engineering
// http://www.ingmo.unimore.it/site/en/home.html
#ifndef MATRIXCOMP_H
#define MATRIXCOMP_H

#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>

struct matrix {
	size_t rows;
	size_t cols;
	double* data;
};

enum linetype
{
	row,
	col
};

/*
	rows: chosen row
	cols: chosen column
	selection: row or column??
*/
struct matrix_selection
{
	size_t rows;
	size_t cols;
	enum linetype selection;
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
extern double *elementAt(struct matrix *matr, size_t rowIndex, size_t colIndex);//Fast way to access an element of a matrix without thinking too much about it.
extern void rowEchelonForm(struct matrix *matr);	//Reduce a mtrix in row echelon form. For more information see https://en.wikipedia.org/wiki/Row_echelon_form
extern void printMatrix(struct matrix *matr, FILE *f);
extern void transposeMatrix(struct matrix *matr);	//Transpose the input matrix.

extern struct matrix *createRandomMatrix(size_t rows, size_t cols, uint16_t randMin, uint16_t randMax);
extern struct matrix *createSequentialMatrix(size_t rows, size_t cols, double start, double step);
//====== Internal functions =======
extern void fillCompMinor(const double *src_matr, double *dst_matr, size_t src_dim, size_t row, size_t col);
extern double det3x3(const double *matr);
extern double laplace(const double *matr, size_t dim, double **compm_matrs, size_t start_dim, const struct matrix_selection *selection);
extern struct matrix_selection findlinewithmorezeros(const struct matrix *matr);
extern void sumTwoRows(struct matrix *matr, size_t rowDest, size_t rowSource, double lambda);		//Sums to rowDest rowSource * lambda
#endif //!MATRIXCOMP_H
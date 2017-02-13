// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi -  from the"Enzo Ferrari" Department of Engineering
// http://www.ingmo.unimore.it/site/en/home.html
#include "matrixcomp.h"
#include <string.h>	

double det(const struct matrix *matr)
{
	if (matr == NULL || matr->rows <= 0 || matr->cols <= 0 || matr->cols != matr->rows) return 0;
	size_t dim = matr->cols;
	if (dim > 3)
	{
		double **compm_matrs = NULL; //allocates once all the space for the complementary minors up to 3x3 (...dim-3)
		compm_matrs = malloc((dim - 3) * sizeof(double*));
		for (int i = 0; i < dim - 3; i++)
		{
			compm_matrs[i] = malloc((i + 3)*(i + 3) * sizeof(double)); //allocate complementary matrix to use: [0] = 3x3... [1] = 4x4... ... [n-1] = (n-1 +3)x(n-1 +3)
		}
		struct matrix_selection selection = findlinewithmorezeros(matr);
		double det = laplace(matr->data, dim, compm_matrs, dim, &selection);
		if (compm_matrs != NULL)
		{
			for (int i = 0; i < dim - 3; i++)
			{
				free(compm_matrs[i]);
			}
			free(compm_matrs);
		}
		return det;
	}
	return laplace(matr->data, dim, NULL, 0, NULL);
}

struct matrix *creatematr(size_t rows, size_t cols)
{
	struct matrix *ris = malloc(sizeof(struct matrix));
	ris->rows = rows;
	ris->cols = cols;
	if(rows*cols > 0)
		ris->data = malloc(rows * cols * sizeof(double));
	return ris;
}

struct matrix *createemptymatr(size_t rows, size_t cols)
{
	struct matrix *ris = malloc(sizeof(struct matrix));
	ris->rows = rows;
	ris->cols = cols;
	if (rows*cols > 0)
		ris->data = calloc(rows * cols * sizeof(double), sizeof(double));
	return ris;
}

void destroymatr(struct matrix* matr)
{
	free(matr->data);
	free(matr);
}

double * elementAt(struct matrix * matr, size_t rowIndex, size_t colIndex)
{
	if (colIndex >= matr->cols || rowIndex >= matr->rows)
		return NULL;

	return matr->data + (matr->cols * rowIndex + colIndex);
}

int matrrow(const struct matrix* matr, size_t rowid, const double *row)
{
	if (matr == NULL && rowid >= matr->rows) return 0;
	for (int c = 0; c < matr->cols; c++)
	{
		matr->data[matr->cols*rowid + c] = row[c];
	}
	return 1;
}

struct matrix *clonematr(const struct matrix* matr)
{
	struct matrix* p = malloc(sizeof(struct matrix));
	p->data = malloc(matr->rows * matr->cols * sizeof(double));
	p->cols = matr->cols;
	p->rows = matr->rows;
	
	memcpy(p->data, matr->data, matr->rows * matr->cols);		//Faster than previous solution.

	return p;
}

struct matrix *creatematrfrom(const double* matr, size_t row, size_t col)
{
	struct matrix* p = creatematr(row, col);
	for (int i = 0; i < row*col; i++)
	{
		p->data[i] = matr[i];
	}
	return p;
}

struct matrix *mulmatr(const struct matrix *lhs, const struct matrix *rhs)
{
	if (lhs->cols != rhs->rows || lhs == NULL || rhs == NULL) return NULL;
	size_t lcols = lhs->cols, rcols = rhs->cols;
	struct matrix* ris = creatematr(lhs->rows, rhs->cols);

	for (size_t r = 0; r < ris->rows; r++)
	{
		for (size_t c = 0; c < ris->cols; c++)
		{
			double somma = 0;
			for (size_t z = 0; z < lcols; z++)
			{
				somma += lhs->data[lcols * r + z] * rhs->data[rcols * z + c];
			}
			ris->data[ris->cols * r + c] = somma;
		}
	}

	return ris;
}

double laplace(const double *matr, size_t dim, double **compm_matrs, size_t start_dim, const struct matrix_selection *selection)
{
	double det = 0;
	if (matr == NULL || dim <= 0) return 0;
	else if (dim == 1) return matr[0];
	else if (dim == 2) return matr[0] * matr[3] - matr[1] * matr[2];
	else if (dim == 3) return det3x3(matr);
	for (int i = 0; i < dim; i++)
	{
		size_t curr_row = 0, curr_col = i;
		if (selection != NULL)
		{
			curr_row = selection->rows + (selection->selection == row ? 0 : i);
			curr_col = selection->cols + (selection->selection == col ? 0 : i);
		}
		if (matr[dim*curr_row+curr_col] != 0)
		{
			const double* current_matr = dim == start_dim ? matr : compm_matrs[dim - 3];
			fillCompMinor(current_matr, compm_matrs[dim-1 - 3], dim, curr_row, curr_col);
			struct matrix_selection r_selection = findlinewithmorezeros(creatematrfrom(compm_matrs[dim-1 -3], curr_row, curr_col));
			det += (i % 2 == 0 ? 1 : -1) * matr[dim*curr_row+curr_col] * laplace(compm_matrs[dim-1 - 3], dim - 1, compm_matrs, start_dim, &r_selection);
		}
	}
	return det;
}

void fillCompMinor(const double *src_matr, double *dst_matr, size_t src_dim, size_t row, size_t col)
{
	if (src_matr == NULL || dst_matr == NULL || src_dim <= 0 || row >= src_dim || col >= src_dim) return;
	for (int i = 0, j = 0; i < src_dim*src_dim; i++)
	{
		size_t r = i / src_dim;
		size_t c = i % src_dim;
		if (r != row && c != col)
		{
			dst_matr[j] = src_matr[i];
			j++;
		}
	}
}

struct matrix* matrcompminor(const struct matrix *matr, int row, int col)
{
	if (matr == NULL || matr->rows == 0 || matr->cols == 0 || row >= matr->rows || col >= matr->cols) return NULL;
	struct matrix *p = creatematr((matr->rows - 1),(matr->cols - 1));
	if (p)
	{
		for (int i = 0, j = 0; i < matr->rows*matr->cols; i++)
		{
			size_t r = i / matr->cols;
			size_t c = i % matr->cols;
			if (r != row && c != col)
			{
				p->data[j] = matr->data[i];
				j++;
			}
		}
	}
	return p;
}

struct matrix_selection findlinewithmorezeros(const struct matrix *matr)
{
	struct matrix_selection selection = { 0, 0, row };
	size_t zeros = 0;
	for (int r = 0; r < matr->rows; r++)
	{
		int count = 0;
		for (int c = 0; c < matr->cols; c++)
		{
			if (matr->data[matr->cols*r + c] == 0)
			{
				count++;
			}
		}
		if (count > zeros)
		{
			zeros = count;
			selection.cols = 0;
			selection.rows = r;
			selection.selection = row;
		}
	}
	for (int c = 0; c < matr->cols; c++)
	{
		int count = 0;
		for (int r = 0; r < matr->rows; r++)
		{
			if (matr->data[matr->cols*r + c] == 0)
			{
				count++;
			}
		}
		if (count > zeros)
		{
			zeros = count;
			selection.cols = c;
			selection.rows = 0;
			selection.selection = col;
		}
	}
	return selection;
}

double det3x3(const double *matr)
{
	if (matr == NULL) return 0;
	return (matr[0] * matr[4] * matr[8]) + (matr[1] * matr[5] * matr[6]) + (matr[2] * matr[3] * matr[7]) - (matr[2] * matr[4] * matr[6]) - (matr[0] * matr[5] * matr[7]) - (matr[1] * matr[3] * matr[8]);
}
// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi - Unimore
#include "matrixcomp.h"

double det(double *matr, int dim)
{
	if (matr == NULL || dim <= 0) return 0;
	if (dim > 3)
	{
		double **compm_matrs = NULL; //allocates once all the space for the complementary minors up to 3x3 (...dim-3)
		compm_matrs = malloc((dim - 3) * sizeof(double*));
		for (int i = 0; i < dim - 3; i++)
		{
			compm_matrs[i] = malloc((i + 3)*(i + 3) * sizeof(double)); //allocate complementary matrix to use: [0] = 3x3... [1] = 4x4... ... [n-1] = (n-1 +3)x(n-1 +3)
		}
		double det = laplace(matr, dim, compm_matrs, dim);
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
	return laplace(matr, dim, NULL, 0);
}

double laplace(double *matr, int dim, double **compm_matrs, int start_dim)
{
	double det = 0;
	if (matr == NULL || dim <= 0) return 0;
	else if (dim == 1) return matr[0];
	else if (dim == 2) return matr[0] * matr[3] - matr[1] * matr[2];
	else if (dim == 3) return det3x3(matr);
	for (int i = 0; i < dim; i++)
	{
		if (matr[i] != 0)
		{
			double* current_matr = dim == start_dim ? matr : compm_matrs[dim - 3];
			fillCompMinor(current_matr, compm_matrs[dim-1 - 3], dim, 0, i);
			det += (i % 2 == 0 ? 1 : -1) * matr[i] * laplace(compm_matrs[dim-1 - 3], dim - 1, compm_matrs, start_dim);
		}
	}
	return det;
}

void fillCompMinor(const double *src_matr, double *dst_matr, int src_dim, int row, int col)
{
	if (src_matr == NULL || dst_matr == NULL || src_dim <= 0 || row < 0 || col < 0 || row >= src_dim || col >= src_dim) return;
	for (int i = 0, j = 0; i < src_dim*src_dim; i++)
	{
		int r = i / src_dim;
		int c = i % src_dim;
		if (r != row && c != col)
		{
			dst_matr[j] = src_matr[i];
			j++;
		}
	}
}

double* compMinor(double *matr, int matr_rows, int matr_cols, int row, int col)
{
	if (matr == NULL || matr_rows == 0 || matr_cols == 0 || row >= matr_rows || col >= matr_cols) return NULL;
	double *n_matr = malloc((matr_cols - 1)*(matr_rows - 1) * sizeof(double));
	if (n_matr)
	{
		for (int i = 0, j = 0; i < matr_rows*matr_cols; i++)
		{
			int r = i / matr_cols;
			int c = i % matr_cols;
			if (r != row && c != col)
			{
				n_matr[j] = matr[i];
				j++;
			}
		}
	}
	return n_matr;
}

double det3x3(double *matr)
{
	if (matr == NULL) return 0;
	return (matr[0] * matr[4] * matr[8]) + (matr[1] * matr[5] * matr[6]) + (matr[2] * matr[3] * matr[7]) - (matr[2] * matr[4] * matr[6]) - (matr[0] * matr[5] * matr[7]) - (matr[1] * matr[3] * matr[8]);
}
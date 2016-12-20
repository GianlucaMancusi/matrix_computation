// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi - Unimore
#include "matrixcomp.h"

double laplace(double *matr, int dim)
{
	double det = 0;
	if (matr == NULL || dim == 0) return 0;
	else if (dim == 1) return matr[0];
	else if (dim == 2) return matr[0] * matr[3] - matr[1] * matr[2];
	else if (dim == 3) return det3x3(matr);
	for (int i = 0; i < dim; i++)
	{
		double* cMinor = compMinor(matr, dim, dim, 0, i);
		det += (i % 2 == 0 ? 1 : -1) * matr[i] * laplace(cMinor, dim - 1);
		free(cMinor);
	}
	return det;
}

double* compMinor(double *matr, int matr_rows, int matr_cols, int row, int col)
{
	if (matr == NULL || matr_rows == 0 || matr_cols == 0 || row >= matr_rows || col >= matr_cols) return 0;
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
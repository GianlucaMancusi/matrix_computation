// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi -  from the"Enzo Ferrari" Department of Engineering
// http://www.ingmo.unimore.it/site/en/home.html
#include "matrixcomp.h"
#include <string.h>	
#include <assert.h>
#include <time.h>
#include <stdbool.h>

//To make this foldable in VS...
#pragma region boolField

struct boolField {
	uint8_t b0 : 1;
	uint8_t b1 : 1;
	uint8_t b2 : 1;
	uint8_t b3 : 1;
	uint8_t b4 : 1;
	uint8_t b5 : 1;
	uint8_t b6 : 1;
	uint8_t b7 : 1;
};

//Returns the value (true if 1 or false if 0) of the bIndex-th bit of the fIndex-th field in the bField array of boolFields
//Produces an undefined behavior if bIndex is not between 0 and 7 (inclusive)!
bool getField(struct boolField *bFields, uint8_t bIndex, size_t fIndex) {
	switch (bIndex)
	{
	case 0:
		return bFields[fIndex].b0 == 1 ? true : false;
	case 1:
		return bFields[fIndex].b1 == 1 ? true : false;
	case 2:
		return bFields[fIndex].b2 == 1 ? true : false;
	case 3:
		return bFields[fIndex].b3 == 1 ? true : false;
	case 4:
		return bFields[fIndex].b4 == 1 ? true : false;
	case 5:
		return bFields[fIndex].b5 == 1 ? true : false;
	case 6:
		return bFields[fIndex].b6 == 1 ? true : false;
	case 7:
		return bFields[fIndex].b7 == 1 ? true : false;
	default:
		break;
	}
}

//Sets the value (true if 1 or false if 0) of the bIndex-th bit of the fIndex-th field in the bField array of boolFields
//Does not do anything if bIndex is not between 0 and 7 (inclusive)!
void setField(struct boolField *bFields, uint8_t bIndex, size_t fIndex, bool value) {
	switch (bIndex)
	{
	case 0:
		bFields[fIndex].b0 = value == true ? 1 : 0; break;
	case 1:
		bFields[fIndex].b1 = value == true ? 1 : 0; break;
	case 2:
		bFields[fIndex].b2 = value == true ? 1 : 0; break;
	case 3:
		bFields[fIndex].b3 = value == true ? 1 : 0; break;
	case 4:
		bFields[fIndex].b4 = value == true ? 1 : 0; break;
	case 5:
		bFields[fIndex].b5 = value == true ? 1 : 0; break;
	case 6:
		bFields[fIndex].b6 = value == true ? 1 : 0; break;
	case 7:
		bFields[fIndex].b7 = value == true ? 1 : 0; break;
	default:
		break;
	}
}

#pragma endregion


//Calculates the determinant with Laplace's rule
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

//Allocates the space for a new matrix 
struct matrix *creatematr(size_t rows, size_t cols)
{
	struct matrix *ris = malloc(sizeof(struct matrix));
	ris->rows = rows;
	ris->cols = cols;
	if (rows*cols > 0)
		ris->data = malloc(rows * cols * sizeof(double));
	return ris;
}
//Allocates the space for a new matrix and initialises its members to 0
struct matrix *createemptymatr(size_t rows, size_t cols)
{
	struct matrix *ris = malloc(sizeof(struct matrix));
	ris->rows = rows;
	ris->cols = cols;
	if (rows*cols > 0)
		ris->data = calloc(rows * cols, sizeof(double));
	return ris;
}

//Frees the memory occupied by an existing matrix
void destroymatr(struct matrix* matr)
{
	free(matr->data);
	free(matr);
}

//Returns the pointer of an element of the matrix given its row and column indicies
double * elementAt(struct matrix * matr, size_t rowIndex, size_t colIndex)
{
	if (colIndex >= matr->cols || rowIndex >= matr->rows)
		return NULL;

	return matr->data + (matr->cols * rowIndex + colIndex);
}

//Converts the index (r,c) in a single index (k) of the element (r,c) of the matrix in matr->data 
size_t rowIndexof(struct matrix *matr, size_t rowIndex, size_t colIndex) {
	return rowIndex * matr->cols + colIndex;
}

//Reduces a matrix in Row Echelon Form. For more information see 
void rowEchelonForm(struct matrix * matr)
{
	size_t cRow = 1;		//Starts from 1 for rows
	//Iterate on columns
	for (size_t c = 0; c < matr->cols && cRow < matr->rows; c++)
	{
		//Find the first element that is non-zero in the column (starting by the c-th row, the others are already reduced)
		size_t firstNonZeroRowIndex = cRow - 1;
		for (; firstNonZeroRowIndex < matr->rows; firstNonZeroRowIndex++)
		{
			if (matr->data[matr->cols * firstNonZeroRowIndex + c] != 0)
				break;
		}

		if (firstNonZeroRowIndex == matr->rows)			//There is no non-zero element in this column.
			continue;

		if (firstNonZeroRowIndex != cRow - 1)
		{
			//Sum to the c-th row the one found (I'm now currently working on the c-th row and the c-th column)
			sumTwoRows(matr, cRow - 1, firstNonZeroRowIndex, 1);
		}

		double a = matr->data[(cRow - 1)*matr->cols + c];
		for (size_t r = cRow; r < matr->rows; r++)
		{
			sumTwoRows(matr, r, cRow - 1, -matr->data[r* matr->cols + c] / a);
		}

		cRow++;
	}
}

//Returns the linear index of the element of the transposed matrix that will be at (elRow,elCol) in the transposed matrix.
size_t inverseTransposedRowIndexOf(struct matrix *matr, size_t elRow, size_t elCol) {
	return elCol * matr->cols + elRow;
}

size_t inverseTransposedLinearIndexOf(struct matrix *matr, size_t linearIndex) {
	return inverseTransposedRowIndexOf(matr, linearIndex / matr->rows, linearIndex % matr->rows);
}

void transposeMatrix(struct matrix * matr)
{
	//Notation tip: (r,c) is the element at row r and column c.
	size_t newRows = matr->cols;
	size_t newCols = matr->rows;

	if (matr->cols == matr->rows) {
		for (size_t r = 0; r < matr->rows; r++)
		{
			for (size_t c = r+1; c < matr->cols; c++)
			{
				double temp = matr->data[r*matr->cols + c];
				matr->data[r*matr->cols + c] = matr->data[c*matr->cols + r];
				matr->data[c*matr->cols + r] = temp;
			}
		}
	}
	else {
		//1. Creating table for element tagging
		size_t elements = matr->cols * matr->rows;
		size_t tableSize = elements / 8 + (elements % 8 != 0);
		struct boolField *auxTable = calloc(tableSize, sizeof(struct boolField));

		//2. Begin to transpose the matrix from the first element
		for (size_t i = 0; i < matr->cols * matr->rows; i++)
		{
			if (getField(auxTable, i % 8, i / 8) == true)
				continue;		//Already truansposed

			size_t linear = i;		//Transpose the matrix starting from i

			size_t transposed = inverseTransposedLinearIndexOf(matr, i);		//What element goes into index i in the transposed matrix?

			if (i == transposed) {
				setField(auxTable, i % 8, i / 8, true);
				continue;		//Fixed element
			}

			double temp = matr->data[i];
			while (i != transposed) {
				if (getField(auxTable, linear % 8, linear / 8))
					continue;		//This element has already been transposed

				//Transpose element linear.
				matr->data[linear] = matr->data[transposed];		//I'm considering the linear element and I put the element that goes in the linear index of the transposed matrix in its place.
				setField(auxTable, linear % 8, linear / 8, true);

				//Update linear and transposed
				linear = transposed;		//Considering now the transposed index (what element goes into the transposed index in the transposed matrix?)
				transposed = inverseTransposedLinearIndexOf(matr, linear);		//This is the index of the element that will go in the traposed index in the tranposed matrix.
			}
			matr->data[linear] = temp;
			setField(auxTable, linear % 8, linear / 8, true);
		}

		matr->cols = newCols;
		matr->rows = newRows;
		free(auxTable);
	}

	
}

//Prints a matrix on a file
void printMatrix(struct matrix * matr, FILE * f)
{
	for (size_t r = 0; r < matr->rows; r++)
	{
		for (size_t c = 0; c < matr->cols; c++)
		{
			fprintf(f, "%5.5lg", matr->data[r*matr->cols + c]);
		}
		fprintf(f, "\n");
	}
	return;
}

//Writes a column on the matrix
int matrrow(const struct matrix* matr, size_t rowid, const double *row)
{
	if (matr == NULL && rowid >= matr->rows) return 0;
	for (int c = 0; c < matr->cols; c++)
	{
		matr->data[matr->cols*rowid + c] = row[c];
	}
	return 1;
}

//Clones a matrix
struct matrix *clonematr(const struct matrix* matr)
{
	struct matrix* p = malloc(sizeof(struct matrix));
	p->data = malloc(matr->rows * matr->cols * sizeof(double));
	p->cols = matr->cols;
	p->rows = matr->rows;

	memcpy(p->data, matr->data, matr->rows * matr->cols);

	return p;
}

//
struct matrix *creatematrfrom(const double* matr, size_t row, size_t col)
{
	struct matrix* p = creatematr(row, col);
	for (int i = 0; i < row*col; i++)
	{
		p->data[i] = matr[i];
	}
	return p;
}

//Matrix product
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
		if (matr[dim*curr_row + curr_col] != 0)
		{
			const double* current_matr = dim == start_dim ? matr : compm_matrs[dim - 3];
			fillCompMinor(current_matr, compm_matrs[dim - 1 - 3], dim, curr_row, curr_col);
			struct matrix_selection r_selection = findlinewithmorezeros(creatematrfrom(compm_matrs[dim - 1 - 3], curr_row, curr_col));
			det += (i % 2 == 0 ? 1 : -1) * matr[dim*curr_row + curr_col] * laplace(compm_matrs[dim - 1 - 3], dim - 1, compm_matrs, start_dim, &r_selection);
		}
	}
	return det;
}

struct matrix * createRandomMatrix(size_t rows, size_t cols, uint16_t randMin, uint16_t randMax)
{
	srand(time(NULL));
	struct matrix *matr = creatematr(rows, cols);
	for (size_t i = 0; i < cols*rows; i++)
	{
		matr->data[i] = (rand() % (randMax - randMin)) + randMin;
	}
	return matr;
}

//Creates a matrix which adjacent elements have all the same absolute difference
struct matrix * createSequentialMatrix(size_t rows, size_t cols, double start, double step)
{
	struct matrix *matr = creatematr(rows, cols);
	double n = start;
	for (size_t i = 0; i < rows*cols; i++)
	{
		matr->data[i] = n;
		n += step;
	}
	return matr;
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
	struct matrix *p = creatematr((matr->rows - 1), (matr->cols - 1));
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

//Sums to rowDest rowSource * lambda
void sumTwoRows(struct matrix * matr, size_t rowDest, size_t rowSource, double lambda)
{
	assert(!(rowDest >= matr->rows || rowSource >= matr->rows));

	if (rowDest >= matr->rows || rowSource >= matr->rows)
		return;

	for (size_t c = 0; c < matr->cols; c++)
	{
		matr->data[rowDest * matr->cols + c] += matr->data[rowSource * matr->cols + c] * lambda;
	}

	return;
}
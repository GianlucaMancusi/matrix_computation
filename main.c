// Gianluca Mancusi, Daniele Manicardi, Gianmarco Lusvardi -  from the"Enzo Ferrari" Department of Engineering
// http://www.ingmo.unimore.it/site/en/home.html
#include "matrixcomp.h"
int main(void)
{
	//EXAMPLES
	//LAPLACE 3x3
	struct matrix *matr = creatematr(3, 3);
	if (matr)
	{
		matr->data[0] = 1;
		matr->data[1] = 2;
		matr->data[2] = 3;
		matr->data[3] = 1;
		matr->data[4] = 1;
		matr->data[5] = 1;
		matr->data[6] = 1;
		matr->data[7] = 2;
		matr->data[8] = 1;
		double ris = det(matr);
		destroymatr(matr);
	}

	//LAPLACE 2x2
	struct matrix *matr2 = creatematr(2, 2);
	if (matr2)
	{
		matr2->data[0] = 1;
		matr2->data[1] = 2;
		matr2->data[2] = 3;
		matr2->data[3] = 1;
		double ris2 = det(matr2);
		destroymatr(matr2);
	}

	//LAPLACE 1x1
	struct matrix *matr3 = creatematr(1, 1);
	if (matr3)
	{
		matr3->data[0] = 4;
		double ris3 = det(matr3);
		destroymatr(matr3);
	}

	//LAPLACE 3x3
	struct matrix *matr4 = creatematr(3, 3);
	if (matr4)
	{
		matr4->data[0] = 1;
		matr4->data[1] = 2;
		matr4->data[2] = 3;
		matr4->data[3] = 1;
		matr4->data[4] = 1;
		matr4->data[5] = 1;
		matr4->data[6] = 1;
		matr4->data[7] = 2;
		matr4->data[8] = 1;
		double ris4 = det(matr4);
		destroymatr(matr4);
	}

	//LAPLACE 4x4
	struct matrix *matr5 = creatematr(4, 4);
	if (matr5)
	{
		matr5->data[0] = 1;
		matr5->data[1] = 2;
		matr5->data[2] = 3;
		matr5->data[3] = 1;
		matr5->data[4] = 0;
		matr5->data[5] = 6;
		matr5->data[6] = 1;
		matr5->data[7] = -2;
		matr5->data[8] = 1;
		matr5->data[9] = 0;
		matr5->data[10] = 2;
		matr5->data[11] = 3;
		matr5->data[12] = 0;
		matr5->data[13] = 1;
		matr5->data[14] = 1;
		matr5->data[15] = 1;
		double ris5 = det(matr5); //DET = -12
		destroymatr(matr5);
	}

	//LAPLACE 5x5 (it's big! i'll write it using mattrow!)
	struct matrix *matr6 = createemptymatr(5, 5);
	if (matr6)
	{
		double row1[] = { 1,3,1,2,1 };
		double row2[] = { 0,1,0,1,1 };
		double row3[] = { 1,-1,1,0,0 };
		double row4[] = { 1,0,2,2,4 };
		double row5[] = { 1,1,0,1,0 };
		matrrow(matr6, 0, row1);
		matrrow(matr6, 1, row2);
		matrrow(matr6, 2, row3);
		matrrow(matr6, 3, row4);
		matrrow(matr6, 4, row5);
		/*
			1  3  1  2  1
			0  1  0  1  1
			1 -1  1  0  0   DET = -2
			1  0  2  2  4
			1  1  0  1  0
		*/
		double ris6 = det(matr6);
		destroymatr(matr6);
	}

	//MULTIPLY TWO MATRICES
	struct matrix *matr7 = createemptymatr(3, 3);
	struct matrix *matr8 = createemptymatr(3, 2);
	if (matr7 && matr8)
	{
		double row1[3] = { 1,0,1 };
		double row2[3] = { 1,5,-1 };
		double row3[3] = { 3,2,0 };
		matrrow(matr7, 0, row1);
		matrrow(matr7, 1, row2);
		matrrow(matr7, 2, row3);

		double row21[2] = { 7,1 };
		double row22[2] = { 1,0 };
		double row23[2] = { 0,4 };
		matrrow(matr8, 0, row21);
		matrrow(matr8, 1, row22);
		matrrow(matr8, 2, row23);
		struct matrix *matr7_8 = mulmatr(matr7, matr8);
		//result matrix:
		/*
			7   5
			12 -3
			23  3
		*/
		destroymatr(matr7);
		destroymatr(matr8);
		destroymatr(matr7_8);
	}
	return 0;
}
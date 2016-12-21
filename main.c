// main.c Example of usage
#include "matrixcomp.h"
int main(void)
{
	//EXAMPLE
	double *p = (double*)malloc(9 * sizeof(double));
	if (p)
	{
		p[0] = 1;
		p[1] = 2;
		p[2] = 3;
		p[3] = 1;
		p[4] = 1;
		p[5] = 1;
		p[6] = 1;
		p[7] = 2;
		p[8] = 1;
		double ris = det(p,3);
		free(p);
	}

	//LAPLACE 2x2
	double *p2 = (double*)malloc(4 * sizeof(double));
	if (p2)
	{
		p2[0] = 1;
		p2[1] = 2;
		p2[2] = 3;
		p2[3] = 1;
		double ris2 = det(p2, 2);
		free(p2);
	}

	//LAPLACE 1x1
	double *p3 = (double*)malloc(sizeof(double));
	if (p3)
	{
		p3[0] = 4;
		double ris3 = det(p3, 1);
		free(p3);
	}

	//LAPLACE 3x3
	double *p4 = (double*)malloc(9 * sizeof(double));
	if (p4)
	{
		p4[0] = 1;
		p4[1] = 2;
		p4[2] = 3;
		p4[3] = 1;
		p4[4] = 1;
		p4[5] = 1;
		p4[6] = 1;
		p4[7] = 2;
		p4[8] = 1;
		double ris4 = det(p4, 3);
		free(p4);
	}

	//LAPLACE 4x4
	double *p5 = (double*)malloc(16 * sizeof(double));
	if (p5)
	{
		p5[0] = 1;
		p5[1] = 2;
		p5[2] = 3;
		p5[3] = 1;
		p5[4] = 0;
		p5[5] = 6;
		p5[6] = 1;
		p5[7] = -2;
		p5[8] = 1;
		p5[9] = 0;
		p5[10] = 2;
		p5[11] = 3;
		p5[12] = 0;
		p5[13] = 1;
		p5[14] = 1;
		p5[15] = 1;
		double ris5 = det(p5, 4);
		free(p5);
	}

	//LAPLACE 5x5
	double *p6 = (double*)malloc(25 * sizeof(double));
	if (p6)
	{
		p6[0] = 1;
		p6[1] = 3;
		p6[2] = 1;
		p6[3] = 2;
		p6[4] = 1;
		p6[5] = 0;
		p6[6] = 1;
		p6[7] = 0;
		p6[8] = 1;
		p6[9] = 1;
		p6[10] = 1;
		p6[11] = -1;
		p6[12] = 1;
		p6[13] = 0;
		p6[14] = 0;
		p6[15] = 1;
		p6[16] = 0;
		p6[17] = 2;
		p6[18] = 2;
		p6[19] = 4;
		p6[20] = 1;
		p6[21] = 1;
		p6[22] = 0;
		p6[23] = 1;
		p6[24] = 0;
		double ris6 = det(p6, 5);
		free(p6);
	}
	return 0;
}
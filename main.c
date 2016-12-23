// main.c ExamPLe of usage
#include "matrixcomp.h"
int main(void)
{
	//EXAMPLES
	//LAPLACE 3x3
	struct matrix *m = creatematrix(3, 3);
	if (m)
	{
		m->data[0] = 1;
		m->data[1] = 2;
		m->data[2] = 3;
		m->data[3] = 1;
		m->data[4] = 1;
		m->data[5] = 1;
		m->data[6] = 1;
		m->data[7] = 2;
		m->data[8] = 1;
		double ris = det(m);
		destroymatrix(m);
	}

	//LAPLACE 2x2
	struct matrix *m2 = creatematrix(2, 2);
	if (m2)
	{
		m2->data[0] = 1;
		m2->data[1] = 2;
		m2->data[2] = 3;
		m2->data[3] = 1;
		double ris2 = det(m2);
		destroymatrix(m2);
	}

	//LAPLACE 1x1
	struct matrix *m3 = creatematrix(1, 1);
	if (m3)
	{
		m3->data[0] = 4;
		double ris3 = det(m3);
		destroymatrix(m3);
	}

	//LAPLACE 3x3
	struct matrix *m4 = creatematrix(3, 3);
	if (m4)
	{
		m4->data[0] = 1;
		m4->data[1] = 2;
		m4->data[2] = 3;
		m4->data[3] = 1;
		m4->data[4] = 1;
		m4->data[5] = 1;
		m4->data[6] = 1;
		m4->data[7] = 2;
		m4->data[8] = 1;
		double ris4 = det(m4);
		destroymatrix(m4);
	}

	//LAPLACE 4x4
	struct matrix *m5 = creatematrix(4, 4);
	if (m5)
	{
		m5->data[0] = 1;
		m5->data[1] = 2;
		m5->data[2] = 3;
		m5->data[3] = 1;
		m5->data[4] = 0;
		m5->data[5] = 6;
		m5->data[6] = 1;
		m5->data[7] = -2;
		m5->data[8] = 1;
		m5->data[9] = 0;
		m5->data[10] = 2;
		m5->data[11] = 3;
		m5->data[12] = 0;
		m5->data[13] = 1;
		m5->data[14] = 1;
		m5->data[15] = 1;
		double ris5 = det(m5);
		destroymatrix(m5);
	}

	//LAPLACE 5x5
	struct matrix *m6 = creatematrix(5, 5);
	if (m6)
	{
		m6->data[0] = 1;
		m6->data[1] = 3;
		m6->data[2] = 1;
		m6->data[3] = 2;
		m6->data[4] = 1;
		m6->data[5] = 0;
		m6->data[6] = 1;
		m6->data[7] = 0;
		m6->data[8] = 1;
		m6->data[9] = 1;
		m6->data[10] = 1;
		m6->data[11] = -1;
		m6->data[12] = 1;
		m6->data[13] = 0;
		m6->data[14] = 0;
		m6->data[15] = 1;
		m6->data[16] = 0;
		m6->data[17] = 2;
		m6->data[18] = 2;
		m6->data[19] = 4;
		m6->data[20] = 1;
		m6->data[21] = 1;
		m6->data[22] = 0;
		m6->data[23] = 1;
		m6->data[24] = 0;
		double ris6 = det(m6, 5);
		destroymatrix(m6);
	}

	//MULTIPLY TWO MATRICES
	struct matrix *m7 = creatematrix(3, 3);
	struct matrix *m8 = creatematrix(3, 2);
	if (m7 && m8)
	{
		m7->data[0] = 1; //3x3
		m7->data[1] = 0;
		m7->data[2] = 1;
		m7->data[3] = 1; 
		m7->data[4] = 5; 
		m7->data[5] = -1;
		m7->data[6] = 3;
		m7->data[7] = 2;
		m7->data[8] = 0;

		m8->data[0] = 7; //3x2
		m8->data[1] = 1;
		m8->data[2] = 1;
		m8->data[3] = 0;
		m8->data[4] = 0;
		m8->data[5] = 4;
		struct matrix *m7_8 = mulmatr(m7, m8);
		//result matrix:
		/*
			7   5
			12 -3
			23  3
		*/
		destroymatrix(m7);
		destroymatrix(m8);
		destroymatrix(m7_8);
	}
	return 0;
}
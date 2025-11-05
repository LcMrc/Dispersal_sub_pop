#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Init(int X1[], int X2[], int X10, int X20, double *t)
{
	X1[0] = X10;
	X1[1] = X10;
	X1[2] = X10;
	X1[3] = X10;
	X1[4] = 0;
	X2[0] = 0;
	X2[1] = 0;
	X2[2] = 0;
	X2[3] = 0;
	X2[4] = X20;	
	*t = 0;
}

int Sum(int X[], int D, int *S)
{

	int i, q = 0;

	for (i = 0; i < D; i++)
	{ 
		q = q+X[i];
	}

	*S = q;

}

double SumDivisionRate(double F, int X1[], int X2[], int K, int D)
{
	double S = 0;
	int i;

	for (i = 0; i < D; i++)
	{
		S = S+F*(1-(double)(X1[i]+X2[i])/K)*X1[i];
	}

	return S;
}

double SumDeathRate(double G, int X[], int D)
{
	double S = 0;
	int i;

	for (i = 0; i < D; i++)
	{
		S = S+G*X[i];
	}

	return S;
}

double SumMigrationRate(double M1, int X1[], double M2, int X2[], int D)
{
	double S = 0;

	S = (D-1)*(M1*X1[0]+M2*X2[0]+M1*X1[1]+M2*X2[1]+M1*X1[2]+M2*X2[2]+M1*X1[3]+M2*X2[3]+M1*X1[4]+M2*X2[4]);

	return S;
}

void TotalTransitionRate(double FA, double GA, double MA, int XA[], double FB, double GB, double MB, int XB[], int K, int D, double *array)
{
	array[0] = SumDivisionRate(FA, XA, XB, K, D);
	array[1] = SumDeathRate(GA, XA, D);
	array[2] = SumDivisionRate(FB, XB, XA, K, D);
	array[3] = SumDeathRate(GB, XB, D);
	array[4] = SumMigrationRate(MA, XA, MB, XB, D);
}

void CumSum(double *array, int index) 
{
    if(index <= 0) return;
    CumSum(array, index-1);
    array[index] += array[index-1];
}

double rnd()
{
    return (double)rand()/(double)RAND_MAX;
}

int SamplingTowerLinear(double *array, int index) 
{
	int i = 0;
	double R1, R2;
	
	R1 = rnd();
	R2 = R1*array[index-1];

	while (array[i] < R2)
	{
		i++;
	}

	return i;
}

double UpdateTime(double *t, double T)
{
	double r, x;
	r = rnd();
	x = 1/r;
	*t += 1/T*log(x);
}

void PartialTransitionRateBis(double *array, int *X1, int *X2, int K, int D, int ir)
{
	int i;

	if (ir == 0)
	{
		for (i = 0; i < D; i++)
		{
			array[i] = (1-(double)(X1[i]+X2[i])/K)*X1[i];
		}
	}
	else if (ir == 1)
	{
		for (i = 0; i < D; i++)
		{
			array[i] = X1[i];
		}
	}
	else if (ir == 2)
	{
		for (i = 0; i < D; i++)
		{
			array[i] = (1-(double)(X1[i]+X2[i])/K)*X2[i];
		}
	}
	else
	{
		for (i = 0; i < D; i++)
		{
			array[i] = X2[i];
		}
	}
}

void PartialTransitionRateTer(double *array, double M1, int *X1, double M2, int *X2)
{
	array[0] = M1*X1[0];
	array[1] = M1*X1[0];
	array[2] = M1*X1[0];
	array[3] = M1*X1[0];
	array[4] = M1*X1[1];
	array[5] = M1*X1[1];
	array[6] = M1*X1[1];
	array[7] = M1*X1[1];
	array[8] = M1*X1[2];
	array[9] = M1*X1[2];
	array[10] = M1*X1[2];
	array[11] = M1*X1[2];
	array[12] = M1*X1[3];
	array[13] = M1*X1[3];
	array[14] = M1*X1[3];
	array[15] = M1*X1[3];
	array[16] = M1*X1[4];
	array[17] = M1*X1[4];
	array[18] = M1*X1[4];
	array[19] = M1*X1[4];
	array[20] = M2*X2[0];
	array[21] = M2*X2[0];
	array[22] = M2*X2[0];
	array[23] = M2*X2[0];
	array[24] = M2*X2[1];
	array[25] = M2*X2[1];
	array[26] = M2*X2[1];
	array[27] = M2*X2[1];
	array[28] = M2*X2[2];
	array[29] = M2*X2[2];
	array[30] = M2*X2[2];
	array[31] = M2*X2[2];
	array[32] = M2*X2[3];
	array[33] = M2*X2[3];
	array[34] = M2*X2[3];
	array[35] = M2*X2[3];	
	array[36] = M2*X2[4];
	array[37] = M2*X2[4];
	array[38] = M2*X2[4];
	array[39] = M2*X2[4];
}

void ReactionBis(int *X1, int *X2, int ir1, int ir2)
{
	if (ir1 == 0)
	{
		X1[ir2] += 1;
	}
	else if (ir1 == 1)
	{
		X1[ir2] -= 1;
	}
	else if (ir1 == 2)
	{
		X2[ir2] += 1;
	}
	else
	{
		X2[ir2] -= 1;
	}
}

void ReactionTer(int *X1, int *X2, int ir2)
{

	if (ir2 == 0)
	{
		X1[0] -= 1;
		X1[1] += 1;
	}
	else if (ir2 == 1)
	{
		X1[0] -= 1;
		X1[2] += 1;
	}
	else if (ir2 == 2)
	{
		X1[0] -= 1;
		X1[3] += 1;
	}
	else if (ir2 == 3)
	{
		X1[0] -= 1;
		X1[4] += 1;
	}
	else if (ir2 == 4)
	{
		X1[1] -= 1;
		X1[0] += 1;
	}
	else if (ir2 == 5)
	{
		X1[1] -= 1;
		X1[2] += 1;
	}
	else if (ir2 == 6)
	{
		X1[1] -= 1;
		X1[3] += 1;
	}
	else if (ir2 == 7)
	{
		X1[1] -= 1;
		X1[4] += 1;
	}
	else if (ir2 == 8)
	{
		X1[2] -= 1;
		X1[0] += 1;
	}
	else if (ir2 == 9)
	{
		X1[2] -= 1;
		X1[1] += 1;
	}
	else if (ir2 == 10)
	{
		X1[2] -= 1;
		X1[3] += 1;
	}
	else if (ir2 == 11)
	{
		X1[2] -= 1;
		X1[4] += 1;
	}
	else if (ir2 == 12)
	{
		X1[3] -= 1;
		X1[0] += 1;
	}
	else if (ir2 == 13)
	{
		X1[3] -= 1;
		X1[1] += 1;
	}
	else if (ir2 == 14)
	{
		X1[3] -= 1;
		X1[2] += 1;
	}
	else if (ir2 == 15)
	{
		X1[3] -= 1;
		X1[4] += 1;
	}
	else if (ir2 == 16)
	{
		X1[4] -= 1;
		X1[0] += 1;
	}
	else if (ir2 == 17)
	{
		X1[4] -= 1;
		X1[1] += 1;
	}
	else if (ir2 == 18)
	{
		X1[4] -= 1;
		X1[2] += 1;
	}
	else if (ir2 == 19)
	{
		X1[4] -= 1;
		X1[3] += 1;
	}
	else if (ir2 == 20)
	{
		X2[0] -= 1;
		X2[1] += 1;
	}
	else if (ir2 == 21)
	{
		X2[0] -= 1;
		X2[2] += 1;
	}
	else if (ir2 == 22)
	{
		X2[0] -= 1;
		X2[3] += 1;
	}
	else if (ir2 == 23)
	{
		X2[0] -= 1;
		X2[4] += 1;
	}
	else if (ir2 == 24)
	{
		X2[1] -= 1;
		X2[0] += 1;
	}
	else if (ir2 == 25)
	{
		X2[1] -= 1;
		X2[2] += 1;
	}
	else if (ir2 == 26)
	{
		X2[1] -= 1;
		X2[3] += 1;
	}
	else if (ir2 == 27)
	{
		X2[1] -= 1;
		X2[4] += 1;
	}
	else if (ir2 == 28)
	{
		X2[2] -= 1;
		X2[0] += 1;
	}
	else if (ir2 == 29)
	{
		X2[2] -= 1;
		X2[1] += 1;
	}
	else if (ir2 == 30)
	{
		X2[2] -= 1;
		X2[3] += 1;
	}
	else if (ir2 == 31)
	{
		X2[2] -= 1;
		X2[4] += 1;
	}
	else if (ir2 == 32)
	{
		X2[3] -= 1;
		X2[0] += 1;
	}
	else if (ir2 == 33)
	{
		X2[3] -= 1;
		X2[1] += 1;
	}
	else if (ir2 == 34)
	{
		X2[3] -= 1;
		X2[2] += 1;
	}
	else if (ir2 == 35)
	{
		X2[3] -= 1;
		X2[4] += 1;
	}
	else if (ir2 == 36)
	{
		X2[4] -= 1;
		X2[0] += 1;
	}
	else if (ir2 == 37)
	{
		X2[4] -= 1;
		X2[1] += 1;
	}
	else if (ir2 == 38)
	{
		X2[4] -= 1;
		X2[2] += 1;
	}
	else if (ir2 == 39)
	{
		X2[4] -= 1;
		X2[3] += 1;
	}

}

int main()
{

	const double FA = 1.0, GA = 0.1, MA = 0.000001, FB = 1.0, GB = 0.1, MB = 0.000001;
	const int K = 100, XA0 = 90, XB0 = 90, D = 5, Nit = 1000;
	double t = 0, Tvect[5], Tvectbis[5], Tvectter[40], tlist[1000];
	int XA[5] = {XA0, XA0, XA0, XA0, 0}, XB[5] = {0, 0, 0, 0, XB0}, k, i, j, ir1, ir2, XASum = 0, XBSum = 0, XBlist[1000];

	for (k = 0; k < Nit; k++)
	{
		Init(XA, XB, XA0, XB0, &t);
		Sum(XA, D, &XASum);
		Sum(XB, D, &XBSum);

		while (XASum != 0 && XBSum != 0)
		{
			TotalTransitionRate(FA, GA, MA, XA, FB, GB, MB, XB, K, D, Tvect);
			CumSum(Tvect, D-1);
			ir1 = SamplingTowerLinear(Tvect, D); 

			if (ir1 < 4)
			{
				PartialTransitionRateBis(Tvectbis, XA, XB, K, D, ir1);
				CumSum(Tvectbis, D-1);
				ir2 = SamplingTowerLinear(Tvectbis, D);
				ReactionBis(XA, XB, ir1, ir2);
			}
			else
			{
				PartialTransitionRateTer(Tvectter, MA, XA, MB, XB);
				CumSum(Tvectter, 40-1);
				ir2 = SamplingTowerLinear(Tvectter, 40);
				ReactionTer(XA, XB, ir2);
			}

			UpdateTime(&t, Tvect[D-1]);

			Sum(XA, D, &XASum);
			Sum(XB, D, &XBSum);	 
		}
		
		XBlist[k] = XB[4];
		tlist[k] = t;
	}

	FILE *fp1;
	fp1 = fopen("XBlist.txt", "w");
	for (i = 0; i < Nit; i++)
	{
		fprintf(fp1, "%d\n", XBlist[i]);
	}
	fclose(fp1); 

	FILE *fp2;
	fp2 = fopen("tlist.txt", "w");
	for (i = 0; i < Nit; i++)
	{
		fprintf(fp2, "%f\n", tlist[i]);
	}
	fclose(fp2);

	return 0;
} 

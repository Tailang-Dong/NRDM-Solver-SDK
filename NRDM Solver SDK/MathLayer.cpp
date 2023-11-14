#include "MathLayer.h"
#include  <iostream>
#include <cmath>
using namespace std;

//Guass-Seidel�����������Դ���������,ϵ������Ϊ����
void Guass_Seidel_Slover(double* x,double** A, double* b, int n, double tolerance)
{
	/******Step1 ���ý�ĵ�����ֵ-�û��Ѹ�******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it_sm = 0;//�����ƴ���

	//�����ʱ�䲻������breakһ��
	//cout << "��������Gusas-Seidel�����������"  << endl;

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				x[i] =x[i] - A[i][j] * x[j] / A[i][i];
				
			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||max���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it_sm = it_sm + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		//cout << "it=" <<it_sm << "  " << "tol=" << tol <<endl;

	} while (tol > tolerance);

	//cout << "it=" << it_sm << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete [] x0;//�ͷŶ�̬����

}

//Guass-Seidel�����������Դ���������-ȱʡԤ�辫�����غ������������������ɾ��
void Guass_Seidel_Slover(double* x, double** A, double* b, int n)
{
	/******Step1 ���ý�ĵ�����ֵ-�û��Ѹ�******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double tolerance=1e-7;//
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it = 0;//�����ƴ���
	cout << "������ʱû��Ԥ�辫��Ҫ��Ĭ��Ԥ�辫��tolerance=" << tolerance << endl;

	//�����ʱ�䲻������breakһ�£������Ȳ��Ӵ˹����ˡ���ֻ��һ������

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x[j] / A[i][i];

			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||1���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete[] x0;//�ͷŶ�̬����

}

//Jacobi�����������Դ���������
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance)
{
	/******Step1 ���ý�ĵ�����ֵ******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it = 0;//�����ƴ���

	//�����ʱ�䲻������breakһ��
	cout << "��������Jacobi�����������" << endl;

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x0[j] / A[i][i];

			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||max���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol > temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete[] x0;//�ͷŶ�̬����

}

//��m>��n,��С���˽�
/*void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//��ʱ����VTV��VTphi
	double** ATA = new double* [n];
	double* ATb = new double[n];
	for (int i = 0; i < n; i++)
		ATA[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < m; k++)
				ATA[i][j] += A[k][i] * A[k][j];
		}

		ATb[i] = 0;
		for (int k = 0; k < m; k++)
			ATb[i] += A[k][i] * b[k];
	}
	//����GuassSeidel�����������theta
	Guass_Seidel_Slover(x, ATA, ATb, n);
	//�ͷ���ʱ����VTV��VTphi
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}*/

//��row=m>��n=column,��С���˽�
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance)
{
	//��ʱ����VTV��VTb
	double** ATA = new double* [n];
	double* ATb = new double[n];
	for (int i = 0; i < n; i++)
		ATA[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < m; k++)
				ATA[i][j] += A[k][i] * A[k][j];
		}

		ATb[i] = 0;
		for (int k = 0; k < m; k++)
			ATb[i] += A[k][i] * b[k];
	}
	//����GuassSeidel�����������x
	Guass_Seidel_Slover(x, ATA, ATb, n, tolerance);
	//�ͷ���ʱ����VTV��VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//��������New Functions
//���������ʽ��ֵDeterminant of Matrix A
double Determinant_of_Matrix(double** A, int n)//n=2or3
{
	if (n == 2)
	{
		return  A[0][0] * A[1][1] - A[0][1] * A[1][0];
	}
	else if (n == 3)
	{
		return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][2] - A[1][1] * A[2][0]);
	}
	else
		return 1.0;
}

//�������n=2or3
void Inverse_of_Matrix(double** InverseA, double** A, int n)
{
	if (n == 2)
	{
		InverseA[0][0] = A[1][1] / Determinant_of_Matrix(A, n);
		InverseA[1][1] = A[0][0] / Determinant_of_Matrix(A, n);
		InverseA[0][1] = -A[0][1] / Determinant_of_Matrix(A, n);
		InverseA[1][0] = -A[1][0] / Determinant_of_Matrix(A, n);
	}
	else if (n == 3)
	{
		InverseA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / Determinant_of_Matrix(A, n);//��������ʽalgebraic complement
		InverseA[1][0] = - (A[1][0] * A[2][2] - A[1][2] * A[2][0]) / Determinant_of_Matrix(A, n);//���ͳ��������,�����κ�������,������ʱ�������봫ֵ�ͷ�
		InverseA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / Determinant_of_Matrix(A, n);//
		InverseA[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / Determinant_of_Matrix(A, n);//
		InverseA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / Determinant_of_Matrix(A, n);
		InverseA[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / Determinant_of_Matrix(A, n);
		InverseA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / Determinant_of_Matrix(A, n);
		InverseA[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / Determinant_of_Matrix(A, n);
		InverseA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / Determinant_of_Matrix(A, n);
	}
	else;
}

//�����������y[m]=A[m][n]x[n] Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n)
{		
	for (int i = 0; i < m; i++)
	{
		y[i] = 0.0;
		for (int j = 0; j < n; j++)
			y[i] = y[i] + A[i][j] * x[j];
	}		
}

//����������Է�����x=A-1* b
void Matrix_Equation_Solver(double* x, double** A, double* b, int n)
{
	double** InverseA = new double* [n];
	for (int i = 0; i < n; i++)
		InverseA[i] = new double[n];
	Inverse_of_Matrix(InverseA, A, n);
	Matrix_Vector_Multiplier(x, InverseA, b, n, n);

	for (int i = 0; i < n; i++)
		delete[] InverseA[i];
	delete[] InverseA;

}

//��m>��n,��С���˽�,���淨
void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//��ʱ����VTV��VTb
	double** ATA = new double* [n];
	double* ATb = new double[n];
	for (int i = 0; i < n; i++)
		ATA[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < m; k++)
				ATA[i][j] += A[k][i] * A[k][j];
		}

		ATb[i] = 0;
		for (int k = 0; k < m; k++)
			ATb[i] += A[k][i] * b[k];
	}
	//����Matrix_Equation_Solver�����������x
	Matrix_Equation_Solver(x, ATA, ATb, n);
	//�ͷ���ʱ����VTV��VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//������Ծ���Cik=AijBjk Matrix Matrix Multiplication
void Matrix_Matrix_Multiplier(double** C, double** A, double** B, int I, int J, int K)
{
	for(int i=0; i<I; i++)
		for (int k = 0; k < K; k++)
		{
			C[i][k] = 0.0;
			for (int j = 0; j < J; j++)
				C[i][k] = C[i][k] + A[i][j] * B[j][k];
		}
}

//����ת�ó��Ծ���VTV=VT[][]*V[m][n]m>n
void VT_V_Multiplier(double** VTV, double** V, int m, int n)
{
	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
		{
			VTV[i][k] = 0.0;
			for (int j = 0; j < m; j++)
				VTV[i][k] = VTV[i][k] + V[j][i] * V[j][k];
		}
}
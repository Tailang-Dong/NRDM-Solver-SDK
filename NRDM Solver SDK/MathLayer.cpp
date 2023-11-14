#include "MathLayer.h"
#include  <iostream>
#include <cmath>
using namespace std;//functions definitions, functions in header "MathLayer.h"

//4-2Guass-Seidel�����������������. Guass-Seidel iterative method solver, function definition.
void Guass_Seidel_Slover(double* x,double** A, double* b, int n, double tolerance)
{
	/******Step1 设置解的迭代初�?用户已给******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变�?*****/
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it_sm = 0;//迭代计次�?

	//如果长时间不收敛，break一�?
	//cout << "您调用了Gusas-Seidel迭代法求解器"  << endl;

	/******Step3 开始迭�?*****/
	do
	{
		//存一下第it步的�?
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				x[i] =x[i] - A[i][j] * x[j] / A[i][i];
				
			}

		}
		//第it次迭代完�?

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||max范数，所有分量的绝对值必须小于容忍公�?
		tol = x[0] - x0[0];//�?分量误差
		tol = abs(tol);//取绝对�?
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it_sm = it_sm + 1;//迭代计次�?1

		/******Step5 输出计算log日志it与tol******/
		//cout << "it=" <<it_sm << "  " << "tol=" << tol <<endl;

	} while (tol > tolerance);

	//cout << "it=" << it_sm << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内�?*****/
	delete [] x0;//释放动态数�?

}

//4-2-2�ع�һ��toleranceȱʡֵΪ1e-7.���԰���������. overload of 4-2,  default tolerance=1e-7, Test version, include output
void Guass_Seidel_Slover(double* x, double** A, double* b, int n)
{
	/******Step1 设置解的迭代初�?用户已给******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变�?*****/
	double tolerance=1e-7;//
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it = 0;//迭代计次�?
	cout << "您调用时没有预设精度要求，默认预设精度tolerance=" << tolerance << endl;

	//如果长时间不收敛，break一下，算了先不加此功能了。这只是一个函�?

	/******Step3 开始迭�?*****/
	do
	{
		//存一下第it步的�?
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x[j] / A[i][i];

			}

		}
		//第it次迭代完�?

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||1范数，所有分量的绝对值必须小于容忍公�?
		tol = x[0] - x0[0];//�?分量误差
		tol = abs(tol);//取绝对�?
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//迭代计次�?1

		/******Step5 输出计算log日志it与tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内�?*****/
	delete[] x0;//释放动态数�?

}

//4-1Jacobi���������������. Jacobi iterative method solver, function definition.
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance)
{
	/******Step1 设置解的迭代初�?*****/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变�?*****/
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it = 0;//迭代计次�?

	//如果长时间不收敛，break一�?
	cout << "您调用了Jacobi迭代法求解器" << endl;

	/******Step3 开始迭�?*****/
	do
	{
		//存一下第it步的�?
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x0[j] / A[i][i];

			}

		}
		//第it次迭代完�?

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||max范数，所有分量的绝对值必须小于容忍公�?
		tol = x[0] - x0[0];//�?分量误差
		tol = abs(tol);//取绝对�?
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol > temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//迭代计次�?1

		/******Step5 输出计算log日志it与tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内�?*****/
	delete[] x0;//释放动态数�?

}

//行m>列n,最小二乘解
/*void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//临时变量VTV和VTphi
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
	//调用GuassSeidel矩阵求解器求theta
	Guass_Seidel_Slover(x, ATA, ATb, n);
	//释放临时变量VTV和VTphi
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}*/

//4-3��m>��n,��С���˽�. row=m>n=column, Matrix equation solver, Least squares solution
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance)
{
	//临时变量VTV和VTb
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
	//调用GuassSeidel矩阵求解器求x
	Guass_Seidel_Slover(x, ATA, ATb, n, tolerance);
	//释放临时变量VTV和VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//新增函数New Functions
//4-4���������ʽ��ֵDeterminant of Matrix A
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

//4-5�������n=2or3. Inverse of Matrix A.
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
		InverseA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / Determinant_of_Matrix(A, n);//代数余子式algebraic complement
		InverseA[1][0] = - (A[1][0] * A[2][2] - A[1][2] * A[2][0]) / Determinant_of_Matrix(A, n);//降低程序计算�?避免多次函数调用,避免临时变量申请传值释�?
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

//4-6�����������y[m]=A[m][n]x[n] Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n)
{		
	for (int i = 0; i < m; i++)
	{
		y[i] = 0.0;
		for (int j = 0; j < n; j++)
			y[i] = y[i] + A[i][j] * x[j];
	}		
}

//4-9����������Է�����x=A-1* b. Matrix Equation Solver, Inverse Matrix Method.
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

//4-10��m>��n,��С���˽�,���淨. row=m>n=column, Matrix equation solver, Least squares solution, Inverse Matrix Method.
void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//临时变量VTV和VTb
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
	//调用Matrix_Equation_Solver矩阵求解器求x
	Matrix_Equation_Solver(x, ATA, ATb, n);
	//释放临时变量VTV和VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//4-7������Ծ���Cik=AijBjk Matrix Matrix Multiplication
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

//4-8����ת�ó��Ծ���VTV=VT[][]*V[m][n]m>n. Matrix Transpose times Matrix.
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
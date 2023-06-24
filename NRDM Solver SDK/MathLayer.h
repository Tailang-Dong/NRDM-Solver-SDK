#pragma once
//����ϵ������ľ��������Ax=b,��������. Matrix equation solver, row = column, function declaration.
//xΪ��������ָ��AΪϵ�������ָ��bΪ����������ָ��nΪ����ά��tolranceΪԤ�辫��.x is the pointer to the solution vector, A is the pointer to the coefficient matrix, b is a pointer to a constant vector, n is the vector dimension, tolrance is preset precision.

//4-1Jacobi����������. Jacobi iterative method solver function declaration.
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel�����������������. Guass-Seidel iterative method solver function declaration.
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2�ع�һ��toleranceȱʡֵΪ1e-7. overload of 4-2,  default tolerance=1e-7.
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3��m>��n,��С���˽�,ATA��������ݲ�. row=m>n=column, Matrix equation solver, Least squares solution
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2�ع���m>��n,��С���˽�ATA��������ݲ�Ĭ��1e-7
//void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//��������New Functions
//4-4���������ʽ��ֵ. Determinant of Matrix A.
double Determinant_of_Matrix(double** A, int n);
//4-5�������,������. Inverse of Matrix A.
void Inverse_of_Matrix(double** InverseA, double** A, int n);
//4-6�����������y=Ax. Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n);
//4-9����������Է�����. Matrix Equation Solver, Inverse Matrix Method.
void Matrix_Equation_Solver(double* x, double** A, double* b, int n);
//4-10��m>��n,��С���˽�,���淨. row=m>n=column, Matrix equation solver, Least squares solution, Inverse Matrix Method.
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//4-7������Ծ���Cik=AijBjk. Matrix Matrix Multiplication
void Matrix_Matrix_Multiplier(double** C, double** A, double** B, int i, int j, int k);

//4-8����ת�ó��Ծ���VTV=VT[][]*V[m][n]m>n. Matrix Transpose times Matrix.
void VT_V_Multiplier(double** VTV, double** V, int m, int n);
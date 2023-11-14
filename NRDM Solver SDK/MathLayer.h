#pragma once
//����ϵ������ľ��������Ax=b
//xΪ��������ָ��AΪϵ�������ָ��bΪ����������ָ��nΪ����ά��tolranceΪԤ�辫��

//4-1Jacobi����������
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel�����������������
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2�ع�һ��toleranceȱʡֵΪ1e-7
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3��m>��n,��С���˽�,ATA��������ݲ�
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2�ع���m>��n,��С���˽�ATA��������ݲ�Ĭ��1e-7
//void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//��������New Functions
//���������ʽ��ֵDeterminant of Matrix A
double Determinant_of_Matrix(double** A, int n);
//�������,������
void Inverse_of_Matrix(double** InverseA, double** A, int n);
//�����������y=Ax Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n);
//����������Է�����
void Matrix_Equation_Solver(double* x, double** A, double* b, int n);
//��m>��n,��С���˽�,���淨
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//������Ծ���Cik=AijBjk Matrix Matrix Multiplication
void Matrix_Matrix_Multiplier(double** C, double** A, double** B, int i, int j, int k);

//����ת�ó��Ծ���VTV=VT[][]*V[m][n]m>n
void VT_V_Multiplier(double** VTV, double** V, int m, int n);
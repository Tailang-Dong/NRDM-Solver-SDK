#pragma once
//方阵系数矩阵的矩阵求解器Ax=b
//x为解向量的指针A为系数矩阵的指针b为常数向量的指针n为向量维数tolrance为预设精度

//4-1Jacobi迭代法声明
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel迭代求解器函数声明
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2重构一下tolerance缺省值为1e-7
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3行m>列n,最小二乘解,ATA方阵迭代容差
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2重构行m>列n,最小二乘解ATA方阵迭代容差默认1e-7
//void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//新增函数New Functions
//矩阵的行列式的值Determinant of Matrix A
double Determinant_of_Matrix(double** A, int n);
//矩阵的逆,待补充
void Inverse_of_Matrix(double** InverseA, double** A, int n);
//矩阵乘以向量y=Ax Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n);
//逆矩阵法求线性方程组
void Matrix_Equation_Solver(double* x, double** A, double* b, int n);
//行m>列n,最小二乘解,求逆法
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//矩阵乘以矩阵Cik=AijBjk Matrix Matrix Multiplication
void Matrix_Matrix_Multiplier(double** C, double** A, double** B, int i, int j, int k);

//矩阵转置乘以矩阵VTV=VT[][]*V[m][n]m>n
void VT_V_Multiplier(double** VTV, double** V, int m, int n);
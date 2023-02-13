#pragma once
//方阵系数矩阵的矩阵求解器Ax=b,函数声明. Matrix equation solver, row = column, function declaration.
//x为解向量的指针A为系数矩阵的指针b为常数向量的指针n为向量维数tolrance为预设精度.x is the pointer to the solution vector, A is the pointer to the coefficient matrix, b is a pointer to a constant vector, n is the vector dimension, tolrance is preset precision.

//4-1Jacobi迭代法声明. Jacobi iterative method solver function declaration.
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel迭代求解器函数声明. Guass-Seidel iterative method solver function declaration.
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2重构一下tolerance缺省值为1e-7. overload of 4-2,  default tolerance=1e-7.
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3行m>列n,最小二乘解,ATA方阵迭代容差. row=m>n=column, Matrix equation solver, Least squares solution
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2重构行m>列n,最小二乘解ATA方阵迭代容差默认1e-7
//void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//新增函数New Functions
//4-4矩阵的行列式的值. Determinant of Matrix A.
double Determinant_of_Matrix(double** A, int n);
//4-5矩阵的逆,待补充. Inverse of Matrix A.
void Inverse_of_Matrix(double** InverseA, double** A, int n);
//4-6矩阵乘以向量y=Ax. Matrix Vector Multiplication
void Matrix_Vector_Multiplier(double* y, double** A, double* x, int m, int n);
//4-9逆矩阵法求线性方程组. Matrix Equation Solver, Inverse Matrix Method.
void Matrix_Equation_Solver(double* x, double** A, double* b, int n);
//4-10行m>列n,最小二乘解,求逆法. row=m>n=column, Matrix equation solver, Least squares solution, Inverse Matrix Method.
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);

//4-7矩阵乘以矩阵Cik=AijBjk. Matrix Matrix Multiplication
void Matrix_Matrix_Multiplier(double** C, double** A, double** B, int i, int j, int k);

//4-8矩阵转置乘以矩阵VTV=VT[][]*V[m][n]m>n. Matrix Transpose times Matrix.
void VT_V_Multiplier(double** VTV, double** V, int m, int n);
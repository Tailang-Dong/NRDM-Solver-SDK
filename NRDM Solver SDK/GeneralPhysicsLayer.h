#pragma once
/*****通用物理层头文件, 函数声明. Header of GeneralPhysicsLayer, Founctions Declarations.******/

/***最常用，被力学层真正调用的.  most commonly called, Called by the mechanics layer***/
/*解线性方程调用求逆法, 是否采用加权最小二乘选项, Dim=2 or 3. call inverse method to solve matrix equation, with weighted least squares option.*/
//3-1标量场的梯度算子. 梯度Gradphi[Dim], 求导矩阵Vec[n_nei][Dim], 临近点属性phi_nei[n_nei], 中心点属性phi_m, 临近点数n_nei, 维度Dim, WeightedLeastSquares为加权最小二乘选项
//3-1gradient of scalar field. gradient=Gradphi[Dim], derivative Matix=Vec[n_nei][Dim], scalar on neighbors=phi_nei[n_nei], scalar on central point=phi_m, number of neighbors=n_nei, 维度Dim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2矢量场的梯度算子.梯度GradE[Dim][Dim], 求导矩阵Vec[n_nei][Dim], 临近点属性E_nei[n_nei][Dim], 中心点属性E_m[Dim], 临近点数n_nei, 维度Dim, WeightedLeastSquares为加权最小二乘选项
//3-2gradient of vector field. gradient=GradE[Dim][Dim], derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, 维度Dim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3矢量场的散度算子.散度DivE, 求导矩阵Vec[n_nei][Dim], 临近点属性E_nei[n_nei][Dim], 中心点属性E_m[Dim], 临近点数n_nei, 维度Dim, WeightedLeastSquares为加权最小二乘选项
//3-3divergence of vector field. gradient=DivE, derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, 维度Dim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4张量场的散度算子.散度GradT[Dim], 求导矩阵Vec[n_nei][Dim], 临近点属性T_nei[n_nei][Dim][Dim], 中心点属性T_m[Dim][Dim], 临近点数n_nei, 维度Dim, WeightedLeastSquares为加权最小二乘选项
//3-4divergence of tensor field. gradient=GradT[Dim], derivative Matix=Vec[n_nei][Dim], tensor on neighbors=T_nei[n_nei][Dim][Dim], tensor on central point=T_m[Dim][Dim], number of neighbors=n_nei, 维度Dim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);

/***重载III, 不被调用早期版本. Overload III, early version, not be called***/
/*解线性方程调用迭代法, 需要预设精度. call iterative method to solve matrix equation, need preset tolerance*/
//3-1-III
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);
//3-2-III
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-3-III
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-4-III
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);

/***重载I,不被调用早期版本. Overload I, early version, not be called***/
/*解线性方程调用求逆法, 无加权最小二乘选项, Dim=2 or 3. call inverse method to solve matrix equation, No weighted least squares option*/
//3-1-I
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-I
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-I
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-I
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim);

/***重载II,已不被调用的早期版本. Overload II, early version, not be called***/
/*解线性方程调用求逆法, 无加权最小二乘选项,Dim=2 or 3.  需要Inverse_VecTVec[Dim][Dim]
call inverse method to solve matrix equation, No weighted least squares option, need Inverse_VecTVec[Dim][Dim]*/
//3-1-II
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-II
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-II
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-II
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim);

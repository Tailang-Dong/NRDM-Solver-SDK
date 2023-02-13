#pragma once
//3-1想求Node m的属性标量phi的梯度, 返回梯度Gradphi[]，指向矩阵为Vec[][]，临近点属性phi_nei[]，所求节点属性phi_m，临近点数n_nei，维度Dim，
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);

//3-2想求Node m的属性矢量E的梯度, 返回梯度GradE[][]，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-3向量的散度想求Node m的属性矢量E的散度, 返回散度DivE，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-4张量的散度想求Node m的属性张量T的散度, 返回散度GradE，指向矩阵为Vec[][]，临近点属性T_nei[][][]，所求节点属性T_m[][]，临近点数n_nei，维度Dim，
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);

/*一次重构-方阵求逆进行线性方程组求解，不需要迭代容差*/
//3-1想求Node m的属性标量phi的梯度, 返回梯度Gradphi[]，指向矩阵为Vec[][]，临近点属性phi_nei[]，所求节点属性phi_m，临近点数n_nei，维度Dim，
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim);

//3-2想求Node m的属性矢量E的梯度, 返回梯度GradE[][]，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-3向量的散度想求Node m的属性矢量E的散度, 返回散度DivE，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-4张量的散度想求Node m的属性张量T的散度, 返回散度GradE，指向矩阵为Vec[][]，临近点属性T_nei[][][]，所求节点属性T_m[][]，临近点数n_nei，维度Dim，
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*二次重构-方阵求逆进行线性方程组求解，不需要迭代容差，预存VTV的逆*/
//3-1想求Node m的属性标量phi的梯度, 返回梯度Gradphi[]，指向矩阵为Vec[][]，临近点属性phi_nei[]，所求节点属性phi_m，临近点数n_nei，维度Dim，
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim);

//3-2想求Node m的属性矢量E的梯度, 返回梯度GradE[][]，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-3向量的散度想求Node m的属性矢量E的散度, 返回散度DivE，指向矩阵为Vec[][]，临近点属性phi_nei[][]，所求节点属性phi_m[]，临近点数n_nei，维度Dim，
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-4张量的散度想求Node m的属性张量T的散度, 返回散度GradE，指向矩阵为Vec[][]，临近点属性T_nei[][][]，所求节点属性T_m[][]，临近点数n_nei，维度Dim，
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*一-2次重构-方阵求逆进行线性方程组求解，不需要迭代容差-增加加权最小二乘选项*/
//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3矢量E的散度
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);
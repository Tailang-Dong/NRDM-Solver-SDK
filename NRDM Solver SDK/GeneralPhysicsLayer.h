#pragma once
/*****ͨ�������ͷ�ļ�, ��������. Header of GeneralPhysicsLayer, Founctions Declarations.******/

/***��ã�����ѧ���������õ�.  most commonly called, Called by the mechanics layer***/
/*�����Է��̵������淨, �Ƿ���ü�Ȩ��С����ѡ��, Dim=2 or 3. call inverse method to solve matrix equation, with weighted least squares option.*/
//3-1���������ݶ�����. �ݶ�Gradphi[Dim], �󵼾���Vec[n_nei][Dim], �ٽ�������phi_nei[n_nei], ���ĵ�����phi_m, �ٽ�����n_nei, ά��Dim, WeightedLeastSquaresΪ��Ȩ��С����ѡ��
//3-1gradient of scalar field. gradient=Gradphi[Dim], derivative Matix=Vec[n_nei][Dim], scalar on neighbors=phi_nei[n_nei], scalar on central point=phi_m, number of neighbors=n_nei, ά��Dim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2ʸ�������ݶ�����.�ݶ�GradE[Dim][Dim], �󵼾���Vec[n_nei][Dim], �ٽ�������E_nei[n_nei][Dim], ���ĵ�����E_m[Dim], �ٽ�����n_nei, ά��Dim, WeightedLeastSquaresΪ��Ȩ��С����ѡ��
//3-2gradient of vector field. gradient=GradE[Dim][Dim], derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, ά��Dim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3ʸ������ɢ������.ɢ��DivE, �󵼾���Vec[n_nei][Dim], �ٽ�������E_nei[n_nei][Dim], ���ĵ�����E_m[Dim], �ٽ�����n_nei, ά��Dim, WeightedLeastSquaresΪ��Ȩ��С����ѡ��
//3-3divergence of vector field. gradient=DivE, derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, ά��Dim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4��������ɢ������.ɢ��GradT[Dim], �󵼾���Vec[n_nei][Dim], �ٽ�������T_nei[n_nei][Dim][Dim], ���ĵ�����T_m[Dim][Dim], �ٽ�����n_nei, ά��Dim, WeightedLeastSquaresΪ��Ȩ��С����ѡ��
//3-4divergence of tensor field. gradient=GradT[Dim], derivative Matix=Vec[n_nei][Dim], tensor on neighbors=T_nei[n_nei][Dim][Dim], tensor on central point=T_m[Dim][Dim], number of neighbors=n_nei, ά��Dim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);

/***����III, �����������ڰ汾. Overload III, early version, not be called***/
/*�����Է��̵��õ�����, ��ҪԤ�辫��. call iterative method to solve matrix equation, need preset tolerance*/
//3-1-III
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);
//3-2-III
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-3-III
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-4-III
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);

/***����I,�����������ڰ汾. Overload I, early version, not be called***/
/*�����Է��̵������淨, �޼�Ȩ��С����ѡ��, Dim=2 or 3. call inverse method to solve matrix equation, No weighted least squares option*/
//3-1-I
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-I
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-I
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-I
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim);

/***����II,�Ѳ������õ����ڰ汾. Overload II, early version, not be called***/
/*�����Է��̵������淨, �޼�Ȩ��С����ѡ��,Dim=2 or 3.  ��ҪInverse_VecTVec[Dim][Dim]
call inverse method to solve matrix equation, No weighted least squares option, need Inverse_VecTVec[Dim][Dim]*/
//3-1-II
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-II
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-II
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-II
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*一-2次重�?方阵求逆进行线性方程组求解，不需要迭代容�?增加加权最小二乘选项*/
//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2矢量E的梯�?E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3矢量E的散�?
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4张量的散�?
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);
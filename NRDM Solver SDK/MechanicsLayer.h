#pragma once
//节点Node2D类
//类内成员数据
class Node2D
{
public:
	double X[2];//1-1节点坐标. coordinate
	bool IfBoundary;//1-2-1内部点/边界点
	double n[2];//1-2-2方向 out-of-boundary normal vector
	double r;//2-1临近域半径 influence radius
	int Nei;//2-2临近点数 number of neighbors
	bool Fixed[2];//3-1是否固定 if fixed
	double Force[2];//3-2载荷 load
	double R[2];//4-1余量 residual
	double displacement[2];//节点位移5-1. displacement
	double epsilon[2][2];//节点应变5-2. strain
	double sigma[2][2];//节点应力5-3. stress
};

//节点Node3D类
//类内成员数据
class Node3D
{
public:
	double X[3];//1-1节点坐标coordinate
	bool IfBoundary;//1-2-1内部点/边界点
	double n[3];//1-2-2方向out-of-boundary normal vector
	double r;//2-1临近域半径influence radius
	int Nei;//2-2临近点数number of neighbors
	bool Fixed[3];//3-1是否固定if fixed
	double Force[3];//3-2载荷load
	double R[3];//4-1余量residual
	double displacement[3];//节点位移5-1displacement
	double epsilon[3][3];//节点应变5-2strain
	double sigma[3][3];//节点应力5-3stress
};

//类外全局变量
//2-3number of Nodes
extern int Num_m;
//2-4Number of Neighbors
extern int* N_Nei;
//2-5GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Node
extern int** GlobalIndex_Nei;
//2-6Vect[m][Nei][2] is the connectivity matrix of the m-th Node
extern double*** Vect;
extern double E;//0-1弹性模量Young's modulus
extern double nu;//0-2泊松比Poisson's ratio
extern double G;//0-3剪切模量shear modulus

//类外函数Functions out of class
//2-1几何方程geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, double tolerance);

//2-2本构方程constitutive equation solver
void Stress_from_Strain(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//2D平面应力2D plane stress
void Stress_from_Strain2D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//3D
void Stress_from_Strain3D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);

//2-3应力算内力evaluate internal force from stress
//2-3-1应力算内力-内部点 on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, double tolerance);
//2-3-2应力算内力-边界点 on boundary nodes
void InForce_from_Stress_bdy(double* inforce, double** sigma_m, double* n, int Dim);

/*overload I. call inverse matrix method to solve matrix equation*/
/*重构 I. 方阵求逆进行线性方程组求解*/
//2-1-I几何方程geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim);

//2-3-1-I应力算内力-内部点. evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim);

/*overload II. call inverse matrix method to solve matrix equation, need Inverse_VecTVec[][]*/
/*重构 II. 方阵求逆进行线性方程组求解,需要Inverse_VecTVec[][]*/
//2-1-II几何方程geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** Inverse_VecTVec, double** disp_nei, double* disp_m, int n_nei, int Dim);

//2-3-1-II应力算内力-内部点evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double** Inverse_VecTVec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim);

/*overload III. call inverse matrix method to solve matrix equation, weighted least square option*/
/*重构III-方阵求逆进行线性方程组求解,具有加权最小二乘选项*/
//2-1-III几何方程geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, bool WeightedLeastSquares);

//2-3-1-III应力算内力-内部点evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, bool WeightedLeastSquares);
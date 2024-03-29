#include <iostream>
#include "MechanicsLayer.h"
#include "GeneralPhysicsLayer.h"
//类外全局变量定义
//2-3number of Nodes
int Num_m = 1;
//2-4Number of Neighbors
int* N_Nei=NULL;
//2-4GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Node
int** GlobalIndex_Nei = NULL;
//2-5Vect[m][Nei][2] is the connectivity matrix of the m-th Node
double*** Vect = NULL;
double E = 2.0e11;//0-1����ģ��Young's modulus
double nu = 0.3;//0-2���ɱ�Poisson's ratio
double G = E / 2 / (1 + nu);//0-3����ģ��shear modulus

//2-1���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, double tolerance)
{
	double** grad_disp = new double* [Dim];//2-1-1位移梯度grad_disp[][]
	for (int i = 0; i < Dim; i++)
		grad_disp[i] = new double[Dim];
	Gradient_of_Vector(grad_disp, Vec, disp_nei, disp_m, n_nei, Dim, tolerance);//2-1-2调用3-2求位移梯�?
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			epsilon[i][j] = (grad_disp[i][j] + grad_disp[j][i]) / 2;//2-1-3根据位移梯度求应�?
	for (int i = 0; i < Dim; i++)//2-1-4释放临时变量grad_disp[][]
		delete[] grad_disp[i];
	delete[] grad_disp;
}

//2-2��������constitutive equation solver
void Stress_from_Strain(double** sigma, double E, double nu, double G, double** epsilon, int Dim)
{
	if (Dim == 2)
	{
		sigma[0][0] = E * (epsilon[0][0] + nu * epsilon[1][1]) / (1 - nu * nu);
		sigma[1][1] = E * (epsilon[1][1] + nu * epsilon[0][0]) / (1 - nu * nu);
		sigma[0][1] = 2 * G * epsilon[0][1];
		sigma[1][0] = 2 * G * epsilon[1][0];
	}
	else
	{
		double traceepsilon = 0;
		for (int i = 0; i < Dim; i++)
			traceepsilon = traceepsilon + epsilon[i][i];
		for (int i = 0; i < Dim; i++)
		{
			for (int j = 0; j < Dim; j++)
			{
				sigma[i][j] = E * epsilon[i][j] / (1 + nu);
				if (i == j)
					sigma[i][j] = sigma[i][j] + nu * E * traceepsilon / (1 + nu) / (1 - 2.0 * nu);
			}
		}
	}
}
//2Dƽ��Ӧ��plane stress
void Stress_from_Strain2D(double** sigma, double E, double nu, double G, double** epsilon, int Dim)
{
	sigma[0][0] = E * (epsilon[0][0] + nu * epsilon[1][1]) / (1 - nu * nu);
	sigma[1][1] = E * (epsilon[1][1] + nu * epsilon[0][0]) / (1 - nu * nu);
	sigma[0][1] = 2 * G * epsilon[0][1];
	sigma[1][0] = 2 * G * epsilon[1][0];
}
//3D
void Stress_from_Strain3D(double** sigma, double E, double nu, double G, double** epsilon, int Dim)
{
	double traceepsilon = 0;
	for (int i = 0; i < Dim; i++)
		traceepsilon = traceepsilon + epsilon[i][i];
	for (int i = 0; i < Dim; i++)
	{
		for (int j = 0; j < Dim; j++)
		{
			sigma[i][j] = E * epsilon[i][j] / (1 + nu);
			if (i == j)
				sigma[i][j] = sigma[i][j] + nu * E * traceepsilon / (1 + nu) / (1 - 2.0 * nu);
		}
	}
}

//2-3-1Ӧ��������-�ڲ���evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, double tolerance)
{
	Divergence_of_Tensor(inforce, Vec, sigma_nei, sigma_m, n_nei, Dim, tolerance);
}
//2-3-2Ӧ��������-�߽��evaluate internal force from stress, on boundary nodes
void InForce_from_Stress_bdy(double* inforce, double** sigma_m, double* n_m, int Dim)
{
	for (int j = 0; j < Dim; j++)
	{
		inforce[j] = 0;
		for (int i = 0; i < Dim; i++)
			inforce[j] = inforce[j] - n_m[i] * sigma_m[i][j];//适用外法向算内力
	}
}

/*overload I. call inverse matrix method to solve matrix equation*/
/*�ع� I. ��������������Է��������*/
//2-1-I���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim)
{
	double** grad_disp = new double* [Dim];//2-1-1位移梯度grad_disp[][]
	for (int i = 0; i < Dim; i++)
		grad_disp[i] = new double[Dim];
	Gradient_of_Vector(grad_disp, Vec, disp_nei, disp_m, n_nei, Dim);//2-1-2调用3-2求位移梯�?
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			epsilon[i][j] = (grad_disp[i][j] + grad_disp[j][i]) / 2;//2-1-3根据位移梯度求应�?
	for (int i = 0; i < Dim; i++)//2-1-4释放临时变量grad_disp[][]
		delete[] grad_disp[i];
	delete[] grad_disp;
}

//2-3-1-IӦ��������-�ڲ���. evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim)
{
	Divergence_of_Tensor(inforce, Vec, sigma_nei, sigma_m, n_nei, Dim);
}

/*overload II. call inverse matrix method to solve matrix equation, need Inverse_VecTVec[][]*/
/*�ع� II. ��������������Է��������,��ҪInverse_VecTVec[][]*/
//2-1-II���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** Inverse_VecTVec, double** disp_nei, double* disp_m, int n_nei, int Dim)
{
	double** grad_disp = new double* [Dim];//2-1-1位移梯度grad_disp[][]
	for (int i = 0; i < Dim; i++)
		grad_disp[i] = new double[Dim];
	Gradient_of_Vector(grad_disp, Vec, Inverse_VecTVec, disp_nei, disp_m, n_nei, Dim);//2-1-2调用3-2求位移梯�?
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			epsilon[i][j] = (grad_disp[i][j] + grad_disp[j][i]) / 2;//2-1-3根据位移梯度求应�?
	for (int i = 0; i < Dim; i++)//2-1-4释放临时变量grad_disp[][]
		delete[] grad_disp[i];
	delete[] grad_disp;
}

//2-3-1-IIӦ��������-�ڲ���evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double** Inverse_VecTVec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim)
{
	Divergence_of_Tensor(inforce, Vec, Inverse_VecTVec, sigma_nei, sigma_m, n_nei, Dim);
}

/*overload III. call inverse matrix method to solve matrix equation, weighted least square option*/
/*�ع�III-��������������Է��������,���м�Ȩ��С����ѡ��*/
//2-1-III���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	double** grad_disp = new double* [Dim];//2-1-1位移梯度grad_disp[][]
	for (int i = 0; i < Dim; i++)
		grad_disp[i] = new double[Dim];
	Gradient_of_Vector(grad_disp, Vec, disp_nei, disp_m, n_nei, Dim, WeightedLeastSquares);//2-1-2调用3-2求位移梯�?
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			epsilon[i][j] = (grad_disp[i][j] + grad_disp[j][i]) / 2;//2-1-3根据位移梯度求应�?
	for (int i = 0; i < Dim; i++)//2-1-4释放临时变量grad_disp[][]
		delete[] grad_disp[i];
	delete[] grad_disp;
}

//2-3-1-IIIӦ��������-�ڲ���evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	Divergence_of_Tensor(inforce, Vec, sigma_nei, sigma_m, n_nei, Dim, WeightedLeastSquares);
}
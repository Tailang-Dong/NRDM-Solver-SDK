#pragma once
//�ڵ�Node2D��
//���ڳ�Ա����
class Node2D
{
public:
	double X[2];//1-1�ڵ�����. coordinate
	bool IfBoundary;//1-2-1�ڲ���/�߽��
	double n[2];//1-2-2���� out-of-boundary normal vector
	double r;//2-1�ٽ���뾶 influence radius
	int Nei;//2-2�ٽ����� number of neighbors
	bool Fixed[2];//3-1�Ƿ�̶� if fixed
	double Force[2];//3-2�غ� load
	double R[2];//4-1���� residual
	double displacement[2];//�ڵ�λ��5-1. displacement
	double epsilon[2][2];//�ڵ�Ӧ��5-2. strain
	double sigma[2][2];//�ڵ�Ӧ��5-3. stress
};

//�ڵ�Node3D��
//���ڳ�Ա����
class Node3D
{
public:
	double X[3];//1-1�ڵ�����coordinate
	bool IfBoundary;//1-2-1�ڲ���/�߽��
	double n[3];//1-2-2����out-of-boundary normal vector
	double r;//2-1�ٽ���뾶influence radius
	int Nei;//2-2�ٽ�����number of neighbors
	bool Fixed[3];//3-1�Ƿ�̶�if fixed
	double Force[3];//3-2�غ�load
	double R[3];//4-1����residual
	double displacement[3];//�ڵ�λ��5-1displacement
	double epsilon[3][3];//�ڵ�Ӧ��5-2strain
	double sigma[3][3];//�ڵ�Ӧ��5-3stress
};

//����ȫ�ֱ���
//2-3number of Nodes
extern int Num_m;
//2-4Number of Neighbors
extern int* N_Nei;
//2-5GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Node
extern int** GlobalIndex_Nei;
//2-6Vect[m][Nei][2] is the connectivity matrix of the m-th Node
extern double*** Vect;
extern double E;//0-1����ģ��Young's modulus
extern double nu;//0-2���ɱ�Poisson's ratio
extern double G;//0-3����ģ��shear modulus

//���⺯��Functions out of class
//2-1���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, double tolerance);

//2-2��������constitutive equation solver
void Stress_from_Strain(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//2Dƽ��Ӧ��2D plane stress
void Stress_from_Strain2D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//3D
void Stress_from_Strain3D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);

//2-3Ӧ��������evaluate internal force from stress
//2-3-1Ӧ��������-�ڲ��� on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, double tolerance);
//2-3-2Ӧ��������-�߽�� on boundary nodes
void InForce_from_Stress_bdy(double* inforce, double** sigma_m, double* n, int Dim);

/*overload I. call inverse matrix method to solve matrix equation*/
/*�ع� I. ��������������Է��������*/
//2-1-I���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim);

//2-3-1-IӦ��������-�ڲ���. evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim);

/*overload II. call inverse matrix method to solve matrix equation, need Inverse_VecTVec[][]*/
/*�ع� II. ��������������Է��������,��ҪInverse_VecTVec[][]*/
//2-1-II���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** Inverse_VecTVec, double** disp_nei, double* disp_m, int n_nei, int Dim);

//2-3-1-IIӦ��������-�ڲ���evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double** Inverse_VecTVec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim);

/*overload III. call inverse matrix method to solve matrix equation, weighted least square option*/
/*�ع�III-��������������Է��������,���м�Ȩ��С����ѡ��*/
//2-1-III���η���geometric equation solver
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, bool WeightedLeastSquares);

//2-3-1-IIIӦ��������-�ڲ���evaluate internal force from stress, on internal nodes
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, bool WeightedLeastSquares);
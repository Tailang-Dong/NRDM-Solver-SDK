#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include <iostream>

//3-1�����ݶ�
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance)
{
	double* dphi = new double[n_nei];//3-1-1����ʱ����dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim, tolerance);//3-1-2����4-3����ݶ�
	delete[] dphi;//3-1-3�ͷ���ʱ����dphi[]

}

//3-2ʸ��E���ݶ� E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-2-1����ʱ����Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//ʸ����ÿ��������Ϊ�������ֱ�����ݶ�
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2����3-1����ÿ���������ݶ�
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute��d_attriV������,�����������
	{
		for (int n = 0; n < n_nei; n++)//Vector�����ĵ�d_attriV����������һ������phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, tolerance);//Vector�����ĵ�d_attriV���������ݶ�
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //��ֵ������GradE[][]�ĵ�d_attriV��
	}
	//3-2-3�ͷžֲ������ڴ�
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3ʸ��E��ɢ��
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-3-1����ʱ����GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2����3-2����ʸ���ݶ�GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, tolerance);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4�ͷ���ʱ�ڴ�GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4������ɢ��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance)
{
	//ÿһ�е�ɢ��T=[E0, E1,]
	//3-4-1��ʱ����Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj=0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//��j��
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2����3-3���j�е�ɢ��
		Divergence_of_Vector(DivEj,Vec, Ej_nei, Ej_m, n_nei, Dim, tolerance);

		//3-4-3������ֵ
		GradT[j] = DivEj;
	}
	//3-4-4�ͷ���ʱ����
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*һ���ع�-��������������Է�������⣬����Ҫ�����ݲ�*/
//3-1�����ݶ�
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1����ʱ����dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2����4-3����ݶ�
	delete[] dphi;//3-1-3�ͷ���ʱ����dphi[]

}

//3-2ʸ��E���ݶ� E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1����ʱ����Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//ʸ����ÿ��������Ϊ�������ֱ�����ݶ�
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2����3-1����ÿ���������ݶ�
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute��d_attriV������,�����������
	{
		for (int n = 0; n < n_nei; n++)//Vector�����ĵ�d_attriV����������һ������phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim);//Vector�����ĵ�d_attriV���������ݶ�
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //��ֵ������GradE[][]�ĵ�d_attriV��
	}
	//3-2-3�ͷžֲ������ڴ�
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3ʸ��E��ɢ��
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1����ʱ����GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2����3-2����ʸ���ݶ�GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4�ͷ���ʱ�ڴ�GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4������ɢ��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//ÿһ�е�ɢ��T=[E0, E1,]
	//3-4-1��ʱ����Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//��j��
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2����3-3���j�е�ɢ��
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3������ֵ
		GradT[j] = DivEj;
	}
	//3-4-4�ͷ���ʱ����
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*�����ع�-��������������Է�������⣬����Ҫ�����ݲԤ��VTV����*/
//3-1����marker m�����Ա���phi���ݶ�, �����ݶ�Gradphi[]����������ΪVec[][]���ٽ�������phi_nei[]������ڵ�����phi_m���ٽ�����n_nei��ά��Dim��
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1����ʱ����dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	double* VecTdphi = new double[Dim];//3-1-1����ʱ����VecTdphi[]
	//VecTdphi[Dim]=Vec[n_nei][Dim]*dphi[n_nei]
	for (int j = 0; j < Dim; j++)
	{
		VecTdphi[j] = 0.0;
		for (int i = 0; i < n_nei; i++)
			VecTdphi[j] = VecTdphi[j] + dphi[i] * Vec[i][j];			
	}
	Matrix_Vector_Multiplier(Gradphi, Inverse_VecTVec, VecTdphi, Dim, Dim);//3-1-2���ú���

	delete[] dphi;//3-1-3�ͷ���ʱ����dphi[]
	delete[] VecTdphi;//3-1-3�ͷ���ʱ����VecTdphi[]

}

//3-2ʸ��E���ݶ� E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1����ʱ����Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//ʸ����ÿ��������Ϊ�������ֱ�����ݶ�
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2����3-1����ÿ���������ݶ�
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute��d_attriV������,�����������
	{
		for (int n = 0; n < n_nei; n++)//Vector�����ĵ�d_attriV����������һ������phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, Inverse_VecTVec, phi_nei, phi_m, n_nei, Dim);//Vector�����ĵ�d_attriV���������ݶ�
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //��ֵ������GradE[][]�ĵ�d_attriV��
	}
	//3-2-3�ͷžֲ������ڴ�
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3ʸ��E��ɢ��
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1����ʱ����GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2����3-2����ʸ���ݶ�GradE
	Gradient_of_Vector(GradE, Vec, Inverse_VecTVec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4�ͷ���ʱ�ڴ�GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4������ɢ��
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//ÿһ�е�ɢ��T=[E0, E1,]
	//3-4-1��ʱ����Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//��j��
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2����3-3���j�е�ɢ��
		Divergence_of_Vector(DivEj, Vec, Inverse_VecTVec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3������ֵ
		GradT[j] = DivEj;
	}
	//3-4-4�ͷ���ʱ����
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*һ-2���ع�-��������������Է�������⣬����Ҫ�����ݲ�-���Ӽ�Ȩ��С����ѡ��*/
//3-1�����ݶ�
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	double* dphi = new double[n_nei];//3-1-1����ʱ����dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	if (WeightedLeastSquares == 1)//������ü�Ȩ��С���˹���
	{
		//׼����Ȩ�ĵ�������ͼ�Ȩ�ĳ�������
		double** WeightVec = new double* [n_nei];
		for (int n = 0; n < n_nei; n++)
			WeightVec[n] = new double[Dim];//���ü�Ȩ�ĵ�������-�ֲ�����
		double DistanceSquare_CenterNeighbor;//��ʱ�ľֲ��������ĵ㵼�ٽ���
		for (int n = 0; n < n_nei; n++)
		{
			DistanceSquare_CenterNeighbor=0;
			for (int i = 0; i < Dim; i++)
				DistanceSquare_CenterNeighbor += Vec[n][i] * Vec[n][i];//�����n��Neighbor�ľ���ƽ��
			dphi[n] = dphi[n] / DistanceSquare_CenterNeighbor;//��Ȩ�ĵ�������
			for (int i = 0; i < Dim; i++)
				WeightVec[n][i] = Vec[n][i] /DistanceSquare_CenterNeighbor;//��Ȩ�ĳ�������
		}
		LeastSquareSolution(Gradphi,WeightVec,dphi,n_nei,Dim);
		for (int n = 0; n < n_nei; n++)//�ͷż�Ȩ�ĵ�������
			delete[] WeightVec[n];
		delete[] WeightVec;
	}
	else
		LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2����4-3����ݶ�//�����ü�Ȩ��С���˹���
	delete[] dphi;//3-1-3�ͷ���ʱ����dphi[]

}

//3-2ʸ��E���ݶ� E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-2-1����ʱ����Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//ʸ����ÿ��������Ϊ�������ֱ�����ݶ�
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2����3-1����ÿ���������ݶ�
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute��d_attriV������,�����������
	{
		for (int n = 0; n < n_nei; n++)//Vector�����ĵ�d_attriV����������һ������phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, WeightedLeastSquares);//Vector�����ĵ�d_attriV���������ݶ�,���þ��м�Ȩ��С����ѡ��ĺ���
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //��ֵ������GradE[][]�ĵ�d_attriV��
	}
	//3-2-3�ͷžֲ������ڴ�
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3ʸ��E��ɢ��
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-3-1����ʱ����GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2����3-2����ʸ���ݶ�GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, WeightedLeastSquares);//���þ��м�Ȩ��С����ѡ��ĺ���
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4�ͷ���ʱ�ڴ�GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4������ɢ��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//ÿһ�е�ɢ��T=[E0, E1,]
	//3-4-1��ʱ����Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//��j��
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2����3-3���j�е�ɢ��
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim, WeightedLeastSquares);//���þ��м�Ȩ��С����ѡ��ĺ���

		//3-4-3������ֵ
		GradT[j] = DivEj;
	}
	//3-4-4�ͷ���ʱ����
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}
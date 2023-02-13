#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include <iostream>

//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance)
{
	double* dphi = new double[n_nei];//3-1-1设临时变量dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim, tolerance);//3-1-2调用4-3求解梯度
	delete[] dphi;//3-1-3释放临时变量dphi[]

}

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-2-1设临时变量Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//矢量的每个分量作为标量来分别计算梯度
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2调用3-1计算每个分量的梯度
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute第d_attriV个分量,逐个分量计算
	{
		for (int n = 0; n < n_nei; n++)//Vector变量的第d_attriV个分量看作一个标量phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, tolerance);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //赋值给张量GradE[][]的第d_attriV行
	}
	//3-2-3释放局部变量内存
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3矢量E的散度
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-3-1设临时变量GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2调用3-2计算矢量梯度GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, tolerance);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4释放临时内存GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance)
{
	//每一列的散度T=[E0, E1,]
	//3-4-1临时变量Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj=0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//第j列
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2调用3-3求第j列的散度
		Divergence_of_Vector(DivEj,Vec, Ej_nei, Ej_m, n_nei, Dim, tolerance);

		//3-4-3分量赋值
		GradT[j] = DivEj;
	}
	//3-4-4释放临时变量
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*一次重构-方阵求逆进行线性方程组求解，不需要迭代容差*/
//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1设临时变量dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2调用4-3求解梯度
	delete[] dphi;//3-1-3释放临时变量dphi[]

}

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1设临时变量Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//矢量的每个分量作为标量来分别计算梯度
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2调用3-1计算每个分量的梯度
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute第d_attriV个分量,逐个分量计算
	{
		for (int n = 0; n < n_nei; n++)//Vector变量的第d_attriV个分量看作一个标量phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //赋值给张量GradE[][]的第d_attriV行
	}
	//3-2-3释放局部变量内存
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3矢量E的散度
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1设临时变量GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2调用3-2计算矢量梯度GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4释放临时内存GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//每一列的散度T=[E0, E1,]
	//3-4-1临时变量Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//第j列
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2调用3-3求第j列的散度
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3分量赋值
		GradT[j] = DivEj;
	}
	//3-4-4释放临时变量
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*二次重构-方阵求逆进行线性方程组求解，不需要迭代容差，预存VTV的逆*/
//3-1想求marker m的属性标量phi的梯度, 返回梯度Gradphi[]，导数矩阵为Vec[][]，临近点属性phi_nei[]，所求节点属性phi_m，临近点数n_nei，维度Dim，
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1设临时变量dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	double* VecTdphi = new double[Dim];//3-1-1设临时变量VecTdphi[]
	//VecTdphi[Dim]=Vec[n_nei][Dim]*dphi[n_nei]
	for (int j = 0; j < Dim; j++)
	{
		VecTdphi[j] = 0.0;
		for (int i = 0; i < n_nei; i++)
			VecTdphi[j] = VecTdphi[j] + dphi[i] * Vec[i][j];			
	}
	Matrix_Vector_Multiplier(Gradphi, Inverse_VecTVec, VecTdphi, Dim, Dim);//3-1-2调用函数

	delete[] dphi;//3-1-3释放临时变量dphi[]
	delete[] VecTdphi;//3-1-3释放临时变量VecTdphi[]

}

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1设临时变量Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//矢量的每个分量作为标量来分别计算梯度
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2调用3-1计算每个分量的梯度
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute第d_attriV个分量,逐个分量计算
	{
		for (int n = 0; n < n_nei; n++)//Vector变量的第d_attriV个分量看作一个标量phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, Inverse_VecTVec, phi_nei, phi_m, n_nei, Dim);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //赋值给张量GradE[][]的第d_attriV行
	}
	//3-2-3释放局部变量内存
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3矢量E的散度
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1设临时变量GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2调用3-2计算矢量梯度GradE
	Gradient_of_Vector(GradE, Vec, Inverse_VecTVec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4释放临时内存GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//每一列的散度T=[E0, E1,]
	//3-4-1临时变量Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//第j列
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2调用3-3求第j列的散度
		Divergence_of_Vector(DivEj, Vec, Inverse_VecTVec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3分量赋值
		GradT[j] = DivEj;
	}
	//3-4-4释放临时变量
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/*一-2次重构-方阵求逆进行线性方程组求解，不需要迭代容差-增加加权最小二乘选项*/
//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	double* dphi = new double[n_nei];//3-1-1设临时变量dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	if (WeightedLeastSquares == 1)//如果采用加权最小二乘构造
	{
		//准备加权的导数矩阵和加权的常数向量
		double** WeightVec = new double* [n_nei];
		for (int n = 0; n < n_nei; n++)
			WeightVec[n] = new double[Dim];//采用加权的导数矩阵-局部变量
		double DistanceSquare_CenterNeighbor;//临时的局部变量中心点导临近点
		for (int n = 0; n < n_nei; n++)
		{
			DistanceSquare_CenterNeighbor=0;
			for (int i = 0; i < Dim; i++)
				DistanceSquare_CenterNeighbor += Vec[n][i] * Vec[n][i];//计算第n个Neighbor的距离平方
			dphi[n] = dphi[n] / DistanceSquare_CenterNeighbor;//加权的导数矩阵
			for (int i = 0; i < Dim; i++)
				WeightVec[n][i] = Vec[n][i] /DistanceSquare_CenterNeighbor;//加权的常数向量
		}
		LeastSquareSolution(Gradphi,WeightVec,dphi,n_nei,Dim);
		for (int n = 0; n < n_nei; n++)//释放加权的导数矩阵
			delete[] WeightVec[n];
		delete[] WeightVec;
	}
	else
		LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2调用4-3求解梯度//不采用加权最小二乘构造
	delete[] dphi;//3-1-3释放临时变量dphi[]

}

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-2-1设临时变量Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//矢量的每个分量作为标量来分别计算梯度
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2调用3-1计算每个分量的梯度
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute第d_attriV个分量,逐个分量计算
	{
		for (int n = 0; n < n_nei; n++)//Vector变量的第d_attriV个分量看作一个标量phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, WeightedLeastSquares);//Vector变量的第d_attriV个分量的梯度,调用具有加权最小二乘选项的函数
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //赋值给张量GradE[][]的第d_attriV行
	}
	//3-2-3释放局部变量内存
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3矢量E的散度
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-3-1设临时变量GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2调用3-2计算矢量梯度GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, WeightedLeastSquares);//调用具有加权最小二乘选项的函数
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4释放临时内存GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//每一列的散度T=[E0, E1,]
	//3-4-1临时变量Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//第j列
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2调用3-3求第j列的散度
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim, WeightedLeastSquares);//调用具有加权最小二乘选项的函数

		//3-4-3分量赋值
		GradT[j] = DivEj;
	}
	//3-4-4释放临时变量
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}
#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include <iostream>

/*****通用物理层源文件, 函数声明. Source of GeneralPhysicsLayer, Founctions Definitions.******/

/***最常用，被力学层真正调用的.  most commonly called, Called by the mechanics layer***/
/*解线性方程调用求逆法, 是否采用加权最小二乘选项, Dim=2 or 3. call inverse method to solve matrix equation, with weighted least squares option.*/
//3-1
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	double* dphi = new double[n_nei];//3-1-1璁句复鏃跺彉閲廳phi[]
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
			DistanceSquare_CenterNeighbor = 0;
			for (int i = 0; i < Dim; i++)
				DistanceSquare_CenterNeighbor += Vec[n][i] * Vec[n][i];//计算第n个Neighbor的距离平方
			dphi[n] = dphi[n] / DistanceSquare_CenterNeighbor;//加权的导数矩阵
			for (int i = 0; i < Dim; i++)
				WeightVec[n][i] = Vec[n][i] / DistanceSquare_CenterNeighbor;//加权的常数向量
		}
		LeastSquareSolution(Gradphi, WeightVec, dphi, n_nei, Dim);
		for (int n = 0; n < n_nei; n++)//释放加权的导数矩阵
			delete[] WeightVec[n];
		delete[] WeightVec;
	}
	else
		LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2调用4-3求解梯度//不采用加权最小二乘构造
	delete[] dphi;//3-1-3閲婃斁涓存椂鍙橀噺dphi[]

}


//3-2
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-2-1璁句复鏃跺彉閲廏radphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//鐭㈤噺鐨勬瘡涓垎閲忎綔涓烘爣閲忔潵鍒嗗埆璁＄畻姊害
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2璋冪敤3-1璁＄畻姣忎釜鍒嗛噺鐨勬搴?
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute绗琩_attriV涓垎閲?閫愪釜鍒嗛噺璁＄畻
	{
		for (int n = 0; n < n_nei; n++)//Vector鍙橀噺鐨勭d_attriV涓垎閲忕湅浣滀竴涓爣閲弍hi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, WeightedLeastSquares);//Vector变量的第d_attriV个分量的梯度,调用具有加权最小二乘选项的函数
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //璧嬪�肩粰寮犻噺GradE[][]鐨勭d_attriV琛?
	}
	//3-2-3閲婃斁灞�閮ㄥ彉閲忓唴瀛?
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//3-3-1璁句复鏃跺彉閲廏radE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2璋冪敤3-2璁＄畻鐭㈤噺姊害GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, WeightedLeastSquares);//调用具有加权最小二乘选项的函数
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4閲婃斁涓存椂鍐呭瓨GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares)
{
	//姣忎竴鍒楃殑鏁ｅ害T=[E0, E1,]
	//3-4-1涓存椂鍙橀噺Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//绗琷鍒?
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2璋冪敤3-3姹傜j鍒楃殑鏁ｅ害
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim, WeightedLeastSquares);//调用具有加权最小二乘选项的函数

		//3-4-3鍒嗛噺璧嬪�?
		GradT[j] = DivEj;
	}
	//3-4-4閲婃斁涓存椂鍙橀噺
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/***重载III, 不被调用早期版本. Overload III, early version, not be called***/
/*解线性方程调用迭代法, 需要预设精度. call iterative method to solve matrix equation, need preset tolerance*/
//3-1-III
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance)
{
	double* dphi = new double[n_nei];//3-1-1璁句复鏃跺彉閲廳phi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim, tolerance);//3-1-2调用4-3求解梯度
	delete[] dphi;//3-1-3閲婃斁涓存椂鍙橀噺dphi[]

}
//3-2-III
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-2-1璁句复鏃跺彉閲廏radphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//鐭㈤噺鐨勬瘡涓垎閲忎綔涓烘爣閲忔潵鍒嗗埆璁＄畻姊害
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2璋冪敤3-1璁＄畻姣忎釜鍒嗛噺鐨勬搴?
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute绗琩_attriV涓垎閲?閫愪釜鍒嗛噺璁＄畻
	{
		for (int n = 0; n < n_nei; n++)//Vector鍙橀噺鐨勭d_attriV涓垎閲忕湅浣滀竴涓爣閲弍hi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, tolerance);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //璧嬪�肩粰寮犻噺GradE[][]鐨勭d_attriV琛?
	}
	//3-2-3閲婃斁灞�閮ㄥ彉閲忓唴瀛?
	delete[] Gradphi;
	delete[] phi_nei;

}
//3-3-III
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-3-1璁句复鏃跺彉閲廏radE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2璋冪敤3-2璁＄畻鐭㈤噺姊害GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, tolerance);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4閲婃斁涓存椂鍐呭瓨GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}
//3-4-III
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance)
{
	//姣忎竴鍒楃殑鏁ｅ害T=[E0, E1,]
	//3-4-1涓存椂鍙橀噺Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj=0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//绗琷鍒?
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2璋冪敤3-3姹傜j鍒楃殑鏁ｅ害
		Divergence_of_Vector(DivEj,Vec, Ej_nei, Ej_m, n_nei, Dim, tolerance);

		//3-4-3鍒嗛噺璧嬪�?
		GradT[j] = DivEj;
	}
	//3-4-4閲婃斁涓存椂鍙橀噺
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/***重载I,不被调用早期版本. Overload I, early version, not be called***/
/*解线性方程调用求逆法, 无加权最小二乘选项, Dim=2 or 3. call inverse method to solve matrix equation, No weighted least squares option*/
//3-1-I
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1璁句复鏃跺彉閲廳phi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim);//3-1-2调用4-3求解梯度
	delete[] dphi;//3-1-3閲婃斁涓存椂鍙橀噺dphi[]
	delete[] VecTdphi;//3-1-3閲婃斁涓存椂鍙橀噺VecTdphi[]

}
//3-2-I
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1璁句复鏃跺彉閲廏radphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//鐭㈤噺鐨勬瘡涓垎閲忎綔涓烘爣閲忔潵鍒嗗埆璁＄畻姊害
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2璋冪敤3-1璁＄畻姣忎釜鍒嗛噺鐨勬搴?
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute绗琩_attriV涓垎閲?閫愪釜鍒嗛噺璁＄畻
	{
		for (int n = 0; n < n_nei; n++)//Vector鍙橀噺鐨勭d_attriV涓垎閲忕湅浣滀竴涓爣閲弍hi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //璧嬪�肩粰寮犻噺GradE[][]鐨勭d_attriV琛?
	}
	//3-2-3閲婃斁灞�閮ㄥ彉閲忓唴瀛?
	delete[] Gradphi;
	delete[] phi_nei;

}
//3-3-I
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1璁句复鏃跺彉閲廏radE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2璋冪敤3-2璁＄畻鐭㈤噺姊害GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4閲婃斁涓存椂鍐呭瓨GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}
//3-4-I
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//姣忎竴鍒楃殑鏁ｅ害T=[E0, E1,]
	//3-4-1涓存椂鍙橀噺Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//绗琷鍒?
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2璋冪敤3-3姹傜j鍒楃殑鏁ｅ害
		Divergence_of_Vector(DivEj, Vec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3鍒嗛噺璧嬪�?
		GradT[j] = DivEj;
	}
	//3-4-4閲婃斁涓存椂鍙橀噺
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}

/***重载II,已不被调用的早期版本. Overload II, early version, not be called***/
/*解线性方程调用求逆法, 无加权最小二乘选项,Dim=2 or 3.  需要Inverse_VecTVec[Dim][Dim]
call inverse method to solve matrix equation, No weighted least squares option, need Inverse_VecTVec[Dim][Dim]*/
//3-1-II
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim)
{
	double* dphi = new double[n_nei];//3-1-1璁句复鏃跺彉閲廳phi[]
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

	delete[] dphi;//3-1-3閲婃斁涓存椂鍙橀噺dphi[]
	delete[] VecTdphi;//3-1-3閲婃斁涓存椂鍙橀噺VecTdphi[]

}
//3-2-II
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-2-1璁句复鏃跺彉閲廏radphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//鐭㈤噺鐨勬瘡涓垎閲忎綔涓烘爣閲忔潵鍒嗗埆璁＄畻姊害
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2璋冪敤3-1璁＄畻姣忎釜鍒嗛噺鐨勬搴?
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute绗琩_attriV涓垎閲?閫愪釜鍒嗛噺璁＄畻
	{
		for (int n = 0; n < n_nei; n++)//Vector鍙橀噺鐨勭d_attriV涓垎閲忕湅浣滀竴涓爣閲弍hi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, Inverse_VecTVec, phi_nei, phi_m, n_nei, Dim);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //璧嬪�肩粰寮犻噺GradE[][]鐨勭d_attriV琛?
	}
	//3-2-3閲婃斁灞�閮ㄥ彉閲忓唴瀛?
	delete[] Gradphi;
	delete[] phi_nei;

}
//3-3-II
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim)
{
	//3-3-1璁句复鏃跺彉閲廏radE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2璋冪敤3-2璁＄畻鐭㈤噺姊害GradE
	Gradient_of_Vector(GradE, Vec, Inverse_VecTVec, E_nei, E_m, n_nei, Dim);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4閲婃斁涓存椂鍐呭瓨GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}
//3-4-II
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim)
{
	//姣忎竴鍒楃殑鏁ｅ害T=[E0, E1,]
	//3-4-1涓存椂鍙橀噺Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj = 0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//绗琷鍒?
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2璋冪敤3-3姹傜j鍒楃殑鏁ｅ害
		Divergence_of_Vector(DivEj, Vec, Inverse_VecTVec, Ej_nei, Ej_m, n_nei, Dim);

		//3-4-3鍒嗛噺璧嬪�?
		GradT[j] = DivEj;
	}
	//3-4-4閲婃斁涓存椂鍙橀噺
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;


}

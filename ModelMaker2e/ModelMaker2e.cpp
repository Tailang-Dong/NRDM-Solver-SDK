//受弯矩的悬臂梁模型的模型制作程序:1)模型制作 2)读取验证 3)理论解输出
//左端固定，右端弯矩
//固定的给理论解
#include <iostream>
#include <fstream>
#include "Marker2D.h"

using namespace std;
//制作model文件
int main()
{
	/******0-制作前准备******/
	//0-1计算域
	double L = 2.0, h = 1.0;//长高
	double Mesh00X = 0.0, Mesh00Y = -h / 2;//坐下角坐标
	//0-2网格和节点
	int nx = 20, ny = 10;//X21个点,Y11个点
	double dh = L / nx;
	Num_m = (nx + 1) * (ny + 1);//Marker数目
	//0-3载荷
	double LoadM = 2.0e5;//Pa
	/******0-制作前准备完成******/


	/******1-Marker声明与设置m[Num_m]******/
	/*****1-1申请Marker数组*****/
	Marker2D* m = new Marker2D[Num_m];

	/*****1-2Marker数据/成员变量设置*****/

	/****1-2-1几何****/
	/***1-2-1-0几何相关准备***/
	//坐标准备
	double** X = new double* [nx + 1];//计算二维点坐标
	double** Y = new double* [nx + 1];
	for (int i = 0; i <= nx; i++)
	{
		X[i] = new double[ny + 1];
		Y[i] = new double[ny + 1];
	}
	for (int i = 0; i <= nx; i++)
	{
		for (int j = 0; j <= ny; j++)
		{
			X[i][j] = dh * i + Mesh00X;
			Y[i][j] = dh * j + Mesh00Y;
		}
	}
	double* X1 = new double[Num_m];//二维变一维
	double* Y1 = new double[Num_m];	
	for (int i = 0; i <= nx; i++)
	{
		for (int j = 0; j <= ny; j++)
		{
			X1[j * (nx + 1) + i] = X[i][j];
			Y1[j * (nx + 1) + i] = Y[i][j];
		}
	}	
	for (int i = 0; i <= nx; i++)//释放二维描述的内存
	{
		delete[] X[i];
		delete[] Y[i];
	}
	delete[] X;
	delete[] Y;
	/***1-2-1-0几何相关准备完成***/
	/***1-2-1-1几何-坐标***/
	/***1-2-1-2-1几何-边界点?***/
	/***1-2-1-2-2几何-边界外法向***/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].X[0] = X1[k];//1-1坐标
		m[k].X[1] = Y1[k];
		if (X1[k] == L || X1[k] == 0 || (Y1[k] == -h / 2) || (Y1[k] == h / 2))//1-2-1边界
			m[k].IfBoundary = 1;
		else
			m[k].IfBoundary = 0;
		//1-2-2边界外法向
		if (X1[k] == 0)//左
		{
			m[k].n[0] = -1.0;
			m[k].n[1] = 0;
		}
		else if (X1[k] == L)//右
		{
			m[k].n[0] = 1.0;
			m[k].n[1] = 0;
		}
		else if (Y1[k] == h/2.0)//上
		{
			m[k].n[0] = 0;
			m[k].n[1] = 1.0;
		}
		else if (Y1[k] == -h/2.0)//下
		{
			m[k].n[0] = 0;
			m[k].n[1] = -1.0;
		}
		else//内部点
		{
			m[k].n[0] = 0;
			m[k].n[1] = 0;
		}
	}
	/****1-2-1几何完成****/

	/****1-2-2离散连接性****/
	/***1-2-2-0离散连接性相关准备***/
	double* r_m = new double[Num_m];
	int* N_Nei = new int[Num_m];
	for (int k = 0; k < Num_m; k++)
	{
		r_m[k] = L / nx;
		N_Nei[k] = 0;
	}
	for (int i = 0; i < Num_m; i++)
	{
		for (int j = 0; j < Num_m; j++)
		{
			double distance2 = (X1[j] - X1[i]) * (X1[j] - X1[i]) + (Y1[j] - Y1[i]) * (Y1[j] - Y1[i]);
			if (distance2 > 0 && distance2 <= (1.00004*r_m[i] * r_m[i]))
				N_Nei[i]++;
		}
	}
	/***1-2-2-1离散连接性-临近域半径***/
	/***1-2-2-2离散连接性-临近点数***/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].r = r_m[k];//2-1离散连接性-临近域半径
		m[k].Nei = N_Nei[k];//2-2离散连接性-临近点数
	}
	/****1-2-2离散连接性完成****/

	/****1-2-3约束和载荷****/
	/***1-2-3-1约束***/
	/***1-2-3-2载荷***/
	for (int k = 0; k < Num_m; k++)
	{
		if (X1[k] == 0)//3-1约束
		{
			m[k].Fixed[0] = 1;//左边固定
			m[k].Fixed[1] = 1;
		}
		else
		{
			m[k].Fixed[0] = 0;
			m[k].Fixed[1] = 0;
		}
		if (X1[k] == L)//3-2载荷
		{
			m[k].Force[0] = 12*LoadM*Y1[k]/(h*h*h);
			m[k].Force[1] = 0;
		}
		else
		{
			m[k].Force[0] = 0;
			m[k].Force[1] = 0;
		}
	}
	/****1-2-3约束和载荷完成****/

	/****1-2-4余量****/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].R[0] = 0;//约束点余量为0
		m[k].R[1] = 0;
	}
	/****1-2-4余量完成****/

	/****1-2-5结果****/
	for (int k = 0; k < Num_m; k++)
	{
		for (int i = 0; i < 2; i++)
		{
			m[k].u[i] = 0;//5-1位移
			for (int j = 0; j < 2; j++)
			{
				m[k].epsilon[i][j] = 0;//5-2应力
				m[k].sigma[i][j] = 0;//5-3内力
			}
		}
	}
	for (int k = 0; k < Num_m; k++)//左边位移
	{
		if (X1[k] == 0)
			m[k].u[1] = -6 * LoadM * (X1[k] * X1[k] + nu * Y1[k] * Y1[k]) / (E * h * h * h);
	}
	/****1-2-5结果完成****/

	/*****1-2Marker数据/成员变量设置完成*****/
	/******1-Marker声明与设置m[Num_m]完成******/
	

	/******2补充计算离散连接性******/
	GlobalIndex_Nei = new int*[Num_m];
	for (int i = 0; i < Num_m; i++)
		GlobalIndex_Nei[i] = new int[N_Nei[i]];
	for (int i = 0; i < Num_m; i++)
	{
		int nei = 0;
		for (int j = 0; j < Num_m; j++)
		{
			double distance2 = (X1[j] - X1[i]) * (X1[j] - X1[i]) + (Y1[j] - Y1[i]) * (Y1[j] - Y1[i]);
			if (distance2 > 0 && distance2 <= (1.00004 * r_m[i] * r_m[i]))
			{
				GlobalIndex_Nei[i][nei] = j;
				nei++;
			}
		}
	}
	/******2补充计算离散连接性完成******/


	/******3文件输出******/
	/*****3-0文件(流)准备*****/
	ofstream mfout,bmfout;//定义文件输出流
	mfout.open("model2.dat");//连接输出文件
	bmfout.open("model2.b", ios::binary);//连接二进制输出文件

	/*****3-1输出Marker的数据*****/
	mfout << Num_m<<"\n";
	for (int k = 0; k < Num_m; k++)
	{
		mfout << m[k].X[0] << " " << m[k].X[1]<<" ";//1-1几何坐标
		mfout << m[k].IfBoundary << " ";//1-2-1几何是否边界
		mfout << m[k].n[0] << " " << m[k].n[1] << " ";//1-2-2几何边界外法向
		mfout << m[k].r << " ";//2-1离散连接性-临近域半径
		mfout << m[k].Nei << " ";//2-2离散连接性-临近点数
		mfout << m[k].Fixed[0] << " "<< m[k].Fixed[1] << " ";//3-1约束-固定?
		mfout << m[k].Force[0] << " " << m[k].Force[1] << " ";//3-2载荷-载荷
		mfout << m[k].R[0] << " " << m[k].R[1] << " ";//4余量
		mfout << m[k].u[0] << " " << m[k].u[1];//5-1结果-位移
		for(int i=0; i<2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfout  << " "<< m[k].epsilon[i][j];
			}
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfout << " "<< m[k].sigma[i][j] ;
			}
		mfout<< "\n";
	}
	//二进制
	bmfout.write((char*)&Num_m, sizeof Num_m);
	bmfout.write((char*)&(m[0]), (long long)Num_m * sizeof(m[0]));

	/*****3-2connectivity数据输出*****/
	for (int k = 0; k < Num_m; k++)
		mfout << N_Nei[k] << "\n";
	for (int i = 0; i < Num_m; i++)
	{
		for (int j = 0; j < N_Nei[i]; j++)
			mfout << GlobalIndex_Nei[i][j] << " ";
		mfout << "\n";
	}
	//二进制
	bmfout.write((char*)N_Nei, (long long)Num_m * sizeof(N_Nei[0]));
	for (int i = 0; i < Num_m; i++)
	{
		bmfout.write((char*)GlobalIndex_Nei[i], (long long)N_Nei[i] * sizeof(GlobalIndex_Nei[i][0]));
	}

	/*****3-3关闭文件释放内存*****/
	mfout.close();
	bmfout.close();
	delete[] m;//释放对象数组内存
	delete[] r_m;
	delete[] N_Nei;
	for (int i = 0; i < Num_m; i++)
		delete[] GlobalIndex_Nei[i];
	delete[] GlobalIndex_Nei;
	/******3文件输出完成******/


	/******4测试读取******/
		/*4-0打开数据文件*/
	ifstream mfin;
	mfin.open("model2.dat");
	/*4-1读取到内存*/
	int num_m;//4-1-1Marker数量
	mfin >> num_m;
	cout << num_m;
	Marker2D* m2 = new Marker2D[num_m];//4-1-2Marker数据
	for (int k = 0; k < num_m; k++)
	{
		mfin >> m2[k].X[0] >> m2[k].X[1];//1-1几何坐标
		mfin >> m2[k].IfBoundary;//1-2-1几何是否边界
		mfin >> m2[k].n[0] >> m2[k].n[1];//1-2-2几何边界外法向
		mfin >> m2[k].r;//2-1离散连接性-临近域半径
		mfin >> m2[k].Nei;//2-2离散连接性-临近点数
		mfin >> m2[k].Fixed[0] >> m2[k].Fixed[1];//3-1约束-固定?
		mfin >> m2[k].Force[0] >> m2[k].Force[1];//3-2载荷-载荷
		mfin >> m2[k].R[0] >> m2[k].R[1];//4余量
		mfin >> m2[k].u[0] >> m2[k].u[1];//5-1结果-位移
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfin >> m2[k].epsilon[i][j];
			}
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfin >> m2[k].sigma[i][j];
			}
	}
	int* n_nei = new int[num_m];//4-1-3connectivity数据
	for (int k = 0; k < num_m; k++)
		mfin >> n_nei[k];
	int** globalindex_nei = new int* [num_m];
	for (int i = 0; i < num_m; i++)
		globalindex_nei[i] = new int [n_nei[i]];
	for (int i = 0; i < num_m; i++)
	{
		for (int j = 0; j < n_nei[i]; j++)
			mfin >> globalindex_nei[i][j];
	}
	/*4-2关闭数据文件*/
	mfin.close();

	/*4-3测试读取数据输出新文件*/
	ofstream mfout2;
	mfout2.open("test.dat");//4-3-0打开文件
	mfout2 << num_m << "\n";//4-3-1输出marker数据
	for (int k = 0; k < num_m; k++)
	{
		mfout2 << m2[k].X[0] << " " << m2[k].X[1] << " ";//1-1几何坐标
		mfout2 << m2[k].IfBoundary << " ";//1-2-1几何是否边界
		mfout2 << m2[k].n[0] << " " << m2[k].n[1] << " ";//1-2-2几何边界外法向
		mfout2 << m2[k].r << " ";//2-1离散连接性-临近域半径
		mfout2 << m2[k].Nei << " ";//2-2离散连接性-临近点数
		mfout2 << m2[k].Fixed[0] << " " << m2[k].Fixed[1] << " ";//3-1约束-固定?
		mfout2 << m2[k].Force[0] << " " << m2[k].Force[1] << " ";//3-2载荷-载荷
		mfout2 << m2[k].R[0] << " " << m2[k].R[1] << " ";//4余量
		mfout2 << m2[k].u[0] << " " << m2[k].u[1];//5-1结果-位移
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfout2 << " " << m2[k].epsilon[i][j];
			}
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				mfout2 << " " << m2[k].sigma[i][j];
			}
		mfout2 << "\n";
	}
	for (int k = 0; k < num_m; k++)//4-3-2输出connectivity数据
		mfout2 << n_nei[k] << "\n";
	for (int i = 0; i < num_m; i++)
	{
		for (int j = 0; j < n_nei[i]; j++)
			mfout2 << globalindex_nei[i][j] << " ";
		mfout2 << "\n";
	}
	mfout2.close();//4-3-3关闭文件释放内存
	delete[] m2;
	delete[] n_nei;
	for (int i = 0; i < num_m; i++)
		delete[] globalindex_nei[i];
	delete[] globalindex_nei;
	/******4测试读取完成******/

	/******5理论解******/
	double* displacementx = new double[Num_m];
	double* displacementy = new double[Num_m];
	double* sigmaxx = new double[Num_m];
	double* sigmayy = new double[Num_m];
	double* sigmaxy = new double[Num_m];
	
	for (int p = 0; p < Num_m; p++)
	{
		displacementx[p] = 12 * LoadM * X1[p] * Y1[p] / (E * h * h * h);
		displacementy[p] = -6 * LoadM * (X1[p] * X1[p] + nu * Y1[p] * Y1[p]) / (E * h * h * h);
		sigmaxx[p] = 12 * LoadM * Y1[p] / (h * h * h);
		sigmayy[p] = 0;
		sigmaxy[p] = 0;
	}
	ofstream fout;
	fout.open("TheoryStress.dat");
	fout << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\n";
	for (int p = 0; p < Num_m; p++)
	{
		fout << X1[p] << " " << Y1[p] << " " << displacementx[p] << " " << displacementy[p] << " " << sigmaxx[p] << " " << sigmayy[p] << " " << sigmaxy[p] << "\n";
	}
	fout.close();

	delete[] X1;//释放一维描述临时内存
	delete[] Y1;
	delete[] sigmaxx;
	delete[] sigmayy;
	delete[] sigmaxy;
	delete[] displacementx;
	delete[] displacementy;
	/******5理论解完成******/

	return 0;
}
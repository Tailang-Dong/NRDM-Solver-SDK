// NRDSolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束. include "main" function. the program starts here and end here.

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "ApplicationLayer.h"
#include "MechanicsLayer.h"
#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include "FileLayer.h"

using namespace std;
using namespace std::filesystem;

int main()
{
	/******0Get ProjectFile & ConfigureFile******/
	/*0-1-ProjectFile*/
	 //1-字符串输入项目文件路径、转化为path、打开读取ProjectFile
	cout << "请输入项目文件路径!Please Input the ProjectFile Path!\n";
	string ProjectFilePathString;//项目文件路径字符串，接收用户输入ProjectFilePathString from user
	getline(std::cin, ProjectFilePathString);
	cout << "您输入的路径为You have input:" << ProjectFilePathString << "\n";
	path ProjectFilePath(ProjectFilePathString);//ProjectFilePath Assignment项目文件路径，使用字符串变量初始化,eg: E:/NRDProject/Test2/Project.nrdmpro
	cout << "ProjectFilePath=" << ProjectFilePath << "\n";
	//项目文件输入流全局已经定义have been definited globally
	InputProjectFile.open(ProjectFilePath);//打开项目文件并open the ProjectFile.
	if (InputProjectFile.is_open())
		cout << "打开项目文件成功！open the ProjectFile successfully!\n";
	else
	{
		cout << "打开项目文件失败！程序结束！failed to open the ProjectFile! program exit!\n";
		exit(0);
	}
	//逐行读取ProjectFile内容，等号右边给全局变量.read the ProjectFile line by line
	string StringFromFileL, StringFromFileR;//临时变量
	istringstream instr(StringFromFileR); //使用字符串流转化
	//等待读取的全局变量
	while (!InputProjectFile.eof())//判断是否达到文件尾部，循环取变量
	{
		getline(InputProjectFile, StringFromFileL, '=');
		getline(InputProjectFile, StringFromFileR);
		istringstream instr(StringFromFileR);
		if (StringFromFileL == "ProjectName")
			ProjectName = StringFromFileR;
		else if (StringFromFileL == "Dimension")
			instr >> Dimension;
		else if (StringFromFileL == "ProjectPath")
			ProjectPath = StringFromFileR;
		else if (StringFromFileL == "BinaryModelFilePath")
			BinaryModelFilePath = StringFromFileR;
		else if (StringFromFileL == "ASCIIModelFilePath")
			ASCIIModelFilePath = StringFromFileR;
		else if (StringFromFileL == "ConfigureFilePath")
			ConfigureFilePath = StringFromFileR;
		else;
	}

	InputProjectFile.close();
	//输出项目文件内容output the content of ProjectFile
	cout << "ProjectName=" << ProjectName << "\n";
	cout << "Dimension=" << Dimension << "\n";
	cout << "ProjectPath=" << ProjectPath << "\n";
	cout << "BinaryModelFilePath=" << BinaryModelFilePath << "\n";
	cout << "ASCIIModelFilePath=" << ASCIIModelFilePath << "\n";
	cout << "ConfigureFilePath=" << ConfigureFilePath << "\n";

	/*0-2-ConfigureFile*/
	//2-打开计算配置文件Open the ConfigureFile
	InputConfigureFile;//项目文件输入流
	InputConfigureFile.open(ConfigureFilePath);//打开配置文件Open the ConfigureFile
	if (InputConfigureFile.is_open())
		cout << "打开配置文件成功！Open the ConfigureFile successfully!\n";
	else
	{
		cout << "打开配置文件失败！程序结束！failed to open the ConfigureFile! program exit!\n";
		exit(0);
	}
	//待读取的全局变量varibles waiting to be read
	while (!InputConfigureFile.eof())//判断是否达到文件尾部，循环取变量
	{
		getline(InputConfigureFile, StringFromFileL, '=');
		getline(InputConfigureFile, StringFromFileR);
		istringstream instr(StringFromFileR);
		if (StringFromFileL == "UseBinaryModel")
			instr >> UseBinaryModel;
		else if (StringFromFileL == "YoungModulus")
			instr >> E;
		else if (StringFromFileL == "PoissonRatio")
			instr >> nu;
		else if (StringFromFileL == "MoveCoefficient")
			instr >> MoveCoefficient;
		else if (StringFromFileL == "ResiToleranceInternal")
			instr >> ResiToleranceInternal;
		else if (StringFromFileL == "ResiToleranceBounadry")
			instr >> ResiToleranceBounadry;
		else if (StringFromFileL == "IfWeightLeastSquares")
			instr >> IfWeightLeastSquares;
		else if (StringFromFileL == "MaxIteration")
			instr >> MaxIteration;
		else if (StringFromFileL == "IfOutputMoveStiffness")
			instr >> IfOutputMoveStiffness;
		else if (StringFromFileL == "RecordRelaxation")
			instr >> RecordRelaxation;
		else if (StringFromFileL == "OutputLongitudinalSection")
			instr >> OutputLongitudinalSection;
		else if (StringFromFileL == "OutputCrossSection")
			instr >> OutputCrossSection;
		else;
	}
	InputConfigureFile.close();
	cout << "UseBinaryModel=" << UseBinaryModel << "\n";
	cout << "YoungModulus=" << E << "\n";
	cout << "PoissonRatio=" << nu << "\n";
	cout << "MoveCoefficient=" << MoveCoefficient << "\n";
	cout << "ResiToleranceInternal=" << ResiToleranceInternal << "\n";
	cout << "ResiToleranceBounadry=" << ResiToleranceBounadry << "\n";
	cout << "IfWeightLeastSquares=" << IfWeightLeastSquares << "\n";
	cout << "MaxIteration=" << MaxIteration << "\n";
	cout << "IfOutputMoveStiffness=" << IfOutputMoveStiffness << "\n";
	cout << "RecordRelaxation=" << RecordRelaxation << "\n";
	cout << "OutputLongitudinalSection=" << OutputLongitudinalSection << "\n";
	cout << "OutputCrossSection=" << OutputCrossSection << "\n";
	G = E / 2 / (1 + nu);//计算剪切模量

	/*0-3-ModelFile*/
	//3-打开模型文件Open the ModelFile
	if (UseBinaryModel == 1)
		mfin.open(BinaryModelFilePath, ios::binary);
	else
		mfin.open(ASCIIModelFilePath);
	if (mfin.is_open())
		cout << "打开模型文件成功！Open the ModelFile successfully!\n";
	else
	{
		cout << "打开模型文件失败！程序结束！failed to open the ModelFile! program exit!\n";
		exit(0);
	}

	//0-4粗糙移动率相关,一般用不上rough moving rate. useless
	double K_in_Mutiply_InfluenceRadiusSqure = E * (3.0 - nu) / (1 - nu * nu) / 2;
	double K_b_Mutiply_InfluenceRadius = E / (1 - nu * nu);
	if (Dimension == 3)
		double K_in_Mutiply_InfluenceRadiusSqure = E * (2.0 - nu) / (1 - nu * nu);
	double L = 2.0;//仅验证算例用得上. useful only for verification examples 1

	cout << "G=" << G << "\n";
	cout << "E=" << E << "\n";
	cout << "nu=" << nu << "\n";

	/******1初始化模型信息initialize the model information******/
	/****1-1Node数量****/
	if (UseBinaryModel)
		mfin.read((char*)&Num_m, sizeof Num_m);
	else
		mfin >> Num_m;
	cout << "Num_m=" << Num_m << "\n";
	/****1-2Node数据node data****/
	m = new Node2D[Num_m];//在应用层提前申请好
	if (UseBinaryModel)
		mfin.read((char*)&(m[0]), (long long)Num_m * sizeof(m[0]));
	else
	{
		for (int k = 0; k < Num_m; k++)
		{
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].X[d];//1-1几何坐标
			mfin >> m[k].IfBoundary;//1-2-1几何是否边界
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].n[d];//1-2-2几何边界外法向
			mfin >> m[k].r;//2-1离散连接性-临近域半径
			mfin >> m[k].Nei;//2-2离散连接性-临近点数
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].Fixed[d];////3-1约束-固定?
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].Force[d];//3-2载荷-载荷
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].R[d];//4余量
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].displacement[d];//5-1结果-位移
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
				{
					mfin >> m[k].epsilon[i][j];//5-2结果-应变
				}
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
				{
					mfin >> m[k].sigma[i][j];//5-3结果-应力
				}
		}
	}
	
	/****1-3connectivity数据****/
	N_Nei = new int[Num_m];
	if (UseBinaryModel)
		mfin.read((char*)N_Nei, (long long)Num_m * sizeof(N_Nei[0]));
	else
	{
		for (int k = 0; k < Num_m; k++)
			mfin >> N_Nei[k];
	}

	GlobalIndex_Nei = new int* [Num_m];
	for (int i = 0; i < Num_m; i++)
		GlobalIndex_Nei[i] = new int[N_Nei[i]];
	if (UseBinaryModel)
	{
		for (int i = 0; i < Num_m; i++)//**********************
		{
			mfin.read((char*)GlobalIndex_Nei[i], (long long)N_Nei[i] * sizeof(GlobalIndex_Nei[i][0]));
		}
	}
	else
	{
		for (int i = 0; i < Num_m; i++)
		{
			for (int j = 0; j < N_Nei[i]; j++)
				mfin >> GlobalIndex_Nei[i][j];
		}
	}

	mfin.close();//关闭数据文件.close the model file
	/****1-4connectivity数据导数矩阵组array of derivative matrices****/
	Vect = new double** [Num_m];
	for (int p = 0; p < Num_m; p++)
		Vect[p] = new double* [N_Nei[p]];
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < N_Nei[p]; i++)
			Vect[p][i] = new double[Dimension];
	for (int p = 0; p < Num_m; p++)
	{
		for (int i = 0; i < N_Nei[p]; i++)
		{
			int glodex = GlobalIndex_Nei[p][i];//p的第i个neighbor的全局序号globalindex
			for (int j = 0; j < Dimension; j++)//Dim=2
				Vect[p][i][j] = m[glodex].X[j] - m[p].X[j];
		}
	}
	/*是否最小二乘且加速,权值向量Weig[Ntotal][N_Nei[p]]*/
	/*计算Inverse_VectTVect*/
	double*** Inverse_VectTVect = new double** [Num_m];//全局,最后释放
	for (int p = 0; p < Num_m; p++)
		Inverse_VectTVect[p] = new double* [Dimension];
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < Dimension; i++)
			Inverse_VectTVect[p][i] = new double[Dimension];
	double** temp_VTV = new double* [Dimension];//临时VTV,计算完释放
	double** temp_InverseVTV = new double* [Dimension];//临时InverseVTV,计算完释放
	for (int i = 0; i < Dimension; i++)
	{
		temp_VTV[i] = new double[Dimension];
		temp_InverseVTV[i] = new double[Dimension];
	}
	double** temp_V = NULL;//临时,计算完释放
	for (int p = 0; p < Num_m; p++)
	{
		temp_V = new double* [N_Nei[p]];//分配V的内存
		for (int i = 0; i < N_Nei[p]; i++)
			temp_V[i] = new double[Dimension];
		for (int i = 0; i < N_Nei[p]; i++)
			for (int j = 0; j < Dimension; j++)
				temp_V[i][j] = Vect[p][i][j];//取V
		VT_V_Multiplier(temp_VTV, temp_V, N_Nei[p], Dimension);//VTV
		Inverse_of_Matrix(temp_InverseVTV, temp_VTV, Dimension);//VTV^-1
		for (int i = 0; i < Dimension; i++)//保存VTV^-1
			for (int j = 0; j < Dimension; j++)
				Inverse_VectTVect[p][i][j] = temp_InverseVTV[i][j];
		for (int i = 0; i < N_Nei[p]; i++)//释放V的内存,下个节点会重新申请
			delete[] temp_V[i];
		delete[] temp_V;
	}
	for (int i = 0; i < Dimension; i++)
	{
		delete[] temp_VTV[i];
		delete[] temp_InverseVTV[i];
	}
	delete[] temp_VTV;//换节点不需要重新申请,循环完再释放
	delete[] temp_InverseVTV;//换节点不需要重新申请,循环完再释放
	/*计算Inverse_VectTVect完成*/

	/****1-5计算移动刚度矩阵K[p][Dim][Dim]****/
	std::cout << "开始计算MoveStiffness\n";
	//1-5-0设矩阵序列K
	double** K = new double* [Num_m];
	for (int p = 0; p < Num_m; p++)
		K[p] = new double [Dimension];
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < Dimension; i++)
			K[p][i] = 0.0;

	for (int P = 0; P < Num_m; P++)
	{
		//1-5-1全部位移设为0
		for (int Q = 0; Q < Num_m; Q++)
			for (int i = 0; i < Dimension; i++)
				m[Q].displacement[i] = 0.0;
		for (int J = 0; J < Dimension; J++)
		{
			//1-5-2p,j位移1.0
			m[P].displacement[J] = 1.0;

			//1-5-3更新全场信息
			// //1-5-3-1申请临时变量
			//申请随时释放的临时变量
			//临时变量-(nei0)临近点个数n_nei#######随时释放不用释放每次重新赋值
			int n_nei = 0;
			//临时变量-(nei1)指向矩阵V###########随时释放
			double** V = NULL;
			//临时变量-(nei2)临近点位移disp_nei#########随时释放
			double** disp_nei = NULL;
			//(nei3)临时临近点应力
			double*** sigma_nei = NULL;
			//申请最后释放的临时变量
			//临时变量-(m1)标记点位移dis_m!!!!!!!!!!!!!!最后释放
			double* disp_m = new double[Dimension];
			//临时变量-(m2)应变epsi!!!!!!!!最后释放
			double** epsi = new double* [Dimension];
			for (int i = 0; i < Dimension; i++)
				epsi[i] = new double[Dimension];
			//(m3)临时变量应力sigm[Dim][Dim]!!!!!!!!!最后释放
			double** sigm = new double* [Dimension];
			for (int i = 0; i < Dimension; i++)
				sigm[i] = new double[Dimension];
			//(m4)边界点外法向!!!!!!!!!最后释放
			double* n = new double[Dimension];
			//(m5)临时标记点内力
			double* inforce = new double[Dimension];

			//1-5-3-2-1更新全场应变&应力
			for (int p = 0; p < Num_m; p++)
			{
				/****2-1根据位移算应变****/
				/**2-1-1临时变量赋值**/
				//临时变量-(nei0)临近点个数n_nei
				n_nei = N_Nei[p];
				//临时变量-(nei1)指向矩阵V
				V = new double* [n_nei];
				for (int i = 0; i < n_nei; i++)
					V[i] = new double[Dimension];
				for (int i = 0; i < n_nei; i++)
					for (int j = 0; j < Dimension; j++)
						V[i][j] = Vect[p][i][j];
				//临时变量-(nei2)临近点位移disp_nei
				disp_nei = new double* [n_nei];
				for (int i = 0; i < n_nei; i++)
					disp_nei[i] = new double[Dimension];
				for (int i = 0; i < n_nei; i++)
				{
					int glodex = GlobalIndex_Nei[p][i];
					for (int d = 0; d < Dimension; d++)
						disp_nei[i][d] = m[glodex].displacement[d];
				}
				//临时变量-(m1)标记点位移dis_m
				for (int d = 0; d < Dimension; d++)
					disp_m[d] = m[p].displacement[d];
				//(m2)epsi[][]局外

				/**2-1-2计算应变存储在临时变量中**/
				Strain_from_Displacement(epsi, V, disp_nei, disp_m, n_nei, Dimension, IfWeightLeastSquares);

				/**2-1-3把应变传值给标记点**/
				for (int i = 0; i < Dimension; i++)//算完传值给m
					for (int j = 0; j < Dimension; j++)
						m[p].epsilon[i][j] = epsi[i][j];

				/**2-1-4释放需要随时释放的临时变量的内存**/
				for (int i = 0; i < n_nei; i++)
					delete[] V[i];
				delete[] V;//(nei1)释放方向矩阵V
				for (int i = 0; i < n_nei; i++)//(nei2)释放临时临近点位移disp_nei
					delete[] disp_nei[i];
				delete[] disp_nei;

				/****2-1根据位移算应变完成****/

				/****2-2应变算应力****/
				/**2-2-1准备临时变量赋值**/
				//(m2)epsi[][]还在
				//(m3)sigm[][]局外变量

				/**2-2-2本构方程计算应力**/
				Stress_from_Strain2D(sigm, E, nu, G, epsi, Dimension);

				/**2-2-3临时变量应力传给标记点**/
				for (int i = 0; i < Dimension; i++)//算完传值给m
					for (int j = 0; j < Dimension; j++)
						m[p].sigma[i][j] = sigm[i][j];
				/**2-2-4需要随时释放的内存**/

				/****2-2应变算应力完成****/
			}

			//1-5-3-2-2更新全场内力&合力
			for (int p = 0; p < Num_m; p++)
			{
				/****2-3应力算内力****/
				/**2-3-1临时变量赋值**/
				n_nei = N_Nei[p];//(nei0)临近点个数n_nei				
				V = new double* [n_nei];
				for (int i = 0; i < n_nei; i++)//(nei1)指向矩阵V
					V[i] = new double[Dimension];
				for (int i = 0; i < n_nei; i++)
					for (int j = 0; j < Dimension; j++)
						V[i][j] = Vect[p][i][j];
				sigma_nei = new double** [n_nei];//(nei3)临时临近点应力
				for (int n = 0; n < n_nei; n++)
					sigma_nei[n] = new double* [Dimension];
				for (int n = 0; n < n_nei; n++)
					for (int i = 0; i < Dimension; i++)
						sigma_nei[n][i] = new double[Dimension];
				for (int n = 0; n < n_nei; n++)
				{
					int glodex = GlobalIndex_Nei[p][n];
					for (int i = 0; i < Dimension; i++)
						for (int j = 0; j < Dimension; j++)
							sigma_nei[n][i][j] = m[glodex].sigma[i][j];
				}
				for (int i = 0; i < Dimension; i++)//(m3)-sigm[][]
					for (int j = 0; j < Dimension; j++)
						sigm[i][j] = m[p].sigma[i][j];
				for (int i = 0; i < Dimension; i++)//(m4)-n[]
					n[i] = m[p].n[i];

				/**2-3-2计算内力**/
				if ((m[p].IfBoundary) == 0)//内部点
				{
					InForce_from_Stress_in(inforce, V, sigma_nei, sigm, n_nei, Dimension, IfWeightLeastSquares);
				}
				else//外部点
				{
					InForce_from_Stress_bdy(inforce, sigm, n, Dimension);
				}
				/**2-3-3释放临时变量**/
				for (int n = 0; n < n_nei; n++)
					for (int i = 0; i < Dimension; i++)
						delete[] sigma_nei[n][i];
				for (int n = 0; n < n_nei; n++)
					delete[] sigma_nei[n];
				delete[] sigma_nei;//(nei3)释放临近点应力				
				for (int i = 0; i < n_nei; i++)
					delete[] V[i];
				delete[] V;//(nei1)释放方向矩阵V在2-1步2-3步用完了
				/****2-3应力算内力完成****/

				/****2-4令余量=内力****/
				for (int i = 0; i < Dimension; i++)
					m[p].R[i] = inforce[i];//不需要固定余量为零//不需要外力
				/****2-4令余量=内力完成****/
			}

			//1-5-3-3最后释放
			delete[] disp_m; //(m1)释放临时标记点位移dis_m
			for (int i = 0; i < Dimension; i++)//(m2)释放临时应变epsi[][]
				delete[] epsi[i];
			delete[] epsi;
			for (int i = 0; i < Dimension; i++)
				delete[] sigm[i];
			delete[] sigm;//(m3)释放标记点应力
			delete[] n;//(m4)释放边界点外法向
			delete[] inforce;//(m5)释放临时内力inforce[Dim]

			//1-5-4计算刚度
			K[P][J] = -m[P].R[J];
			//1-5-5位移归零
			m[P].displacement[J] = 0.0;
		}
	}
	//1-5-6输出移动刚度
	if (IfOutputMoveStiffness == 1)
	{
		ofstream Kout;
		Kout.open(ProjectPath/"results/K.dat");
		if(Dimension==2)
			Kout << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"MoveStiffnessX\"" << "\"MoveStiffnessY\"" << "\n";
		else
			Kout << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"Z\"" << "\"MoveStiffnessX\"" << "\"MoveStiffnessY\"" << "\"MoveStiffnessZ\"" << "\n";
		for (int P = 0; P < Num_m; P++)
		{
			for (int i = 0; i < Dimension; i++)
			{
				Kout << m[P].X[i] << " ";
			}
			for (int i = 0; i < Dimension; i++)
			{
				Kout << K[P][i] << " ";
			}
			Kout << "\n";
		}
		Kout.close();
	}


	//1-5-7消除影响
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < Dimension; i++)
			for (int j = 0; j < Dimension; j++)
				m[p].sigma[i][j] = 0.0;

	/****1-5计算移动刚度矩阵K[p][Dim][Dim]结束****/

	/****1-6移动率计算****/
	//Rough Stiffness
	double** MoveRate = new double* [Num_m];
	for (int p = 0; p < Num_m; p++)
		MoveRate[p] = new double[Dimension];
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < Dimension; i++)
		{
			MoveRate[p][i] = MoveCoefficient / K[p][i];
			if (IfUseRoughStiffness == 1)
			{
				if (m[p].IfBoundary == 0)
					MoveRate[p][i] = m[p].r * m[p].r * MoveCoefficient_in / K_in_Mutiply_InfluenceRadiusSqure;
				else
					MoveRate[p][i] = m[p].r * MoveCoefficient_b / K_b_Mutiply_InfluenceRadius;
			}
		}

	/******1初始化模型信息完成******/
	/*部分文件设置*/
	//收敛记录文件
	Itfout.open(ProjectPath/"results/Iteration.dat");
	Itfout << "VARIABLES=" << "\"Iteration\"" << "\"N_Residuals\"" << "\"Rate_Residuals\"" << "\"SpecificR_Max\"" << "\n";
	//松弛记录文件
	if (RecordRelaxation == 1)
	{
		RelaxOut.open(ProjectPath / "results/Deformation_Iteration.dat");
		RelaxOut << "TITLE = \"Results of Deformation\"" << "\n";//file header
		RelaxOut << "FILETYPE = FULL" << "\n";
		RelaxOut << "variables= \"x\" \"y\" \"DisplacementX\" \"DisplacementY\" \"sigmaXX\" \"sigmaYY\" \"sigmaXY\" \"ResiX\" \"ResiY\" \"ReX\" \"ReY\"\n";
	}

	/******2~4主循环Main Loop******/
	do
	{
		int Residuals = 0;
		/******2更新节点信息******/
			// //2-1申请临时变量
			//申请随时释放的临时变量
			//临时变量-(nei0)临近点个数n_nei#######随时释放不用释放每次重新赋值
		int n_nei = 0;
		//临时变量-(nei1)导数矩阵V###########随时释放
		double** V = NULL;
		double** Inverse_VTV = NULL;//VTV^-1
		//临时变量-(nei2)临近点位移disp_nei#########随时释放
		double** disp_nei = NULL;
		//(nei3)临时临近点应力
		double*** sigma_nei = NULL;
		//申请最后释放的临时变量
		//临时变量-(m1)节点位移dis_m!!!!!!!!!!!!!!最后释放
		double* disp_m = new double[Dimension];
		//临时变量-(m2)应变epsi!!!!!!!!最后释放
		double** epsi = new double* [Dimension];
		for (int i = 0; i < Dimension; i++)
			epsi[i] = new double[Dimension];
		//(m3)临时变量应力sigm[Dim][Dim]!!!!!!!!!最后释放
		double** sigm = new double* [Dimension];
		for (int i = 0; i < Dimension; i++)
			sigm[i] = new double[Dimension];
		//(m4)边界点外法向!!!!!!!!!最后释放
		double* n = new double[Dimension];
		//(m5)临时节点内力
		double* inforce = new double[Dimension];

		//2-2-1更新全场应变&应力
		for (int p = 0; p < Num_m; p++)
		{
			/****2-1根据位移算应变****/
			/**2-1-1临时变量赋值**/
			//临时变量-(nei0)临近点个数n_nei
			n_nei = N_Nei[p];
			//临时变量-(nei1)导数矩阵V
			V = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)
				V[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
				for (int j = 0; j < Dimension; j++)
					V[i][j] = Vect[p][i][j];
			Inverse_VTV = new double* [Dimension];
			for (int i = 0; i < Dimension; i++)
				Inverse_VTV[i] = new double[Dimension];
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
					Inverse_VTV[i][j] = Inverse_VectTVect[p][i][j];
			//临时变量-(nei2)临近点位移disp_nei
			disp_nei = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)
				disp_nei[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
			{
				int glodex = GlobalIndex_Nei[p][i];
				for (int d = 0; d < Dimension; d++)
					disp_nei[i][d] = m[glodex].displacement[d];
			}
			//IF Weight 
			// 
			//临时变量-(m1)节点位移dis_m
			for (int d = 0; d < Dimension; d++)
				disp_m[d] = m[p].displacement[d];
			//(m2)epsi[][]局外

			/**2-1-2计算应变存储在临时变量中**/
			if(IfPrestoreInverseVTV==0)//不预存VTV^-1
			Strain_from_Displacement(epsi, V, disp_nei, disp_m, n_nei, Dimension, IfWeightLeastSquares);
			else
				Strain_from_Displacement(epsi, V, Inverse_VTV, disp_nei, disp_m, n_nei, Dimension);
			

			/**2-1-3把应变传值给节点**/
			for (int i = 0; i < Dimension; i++)//算完传值给m
				for (int j = 0; j < Dimension; j++)
					m[p].epsilon[i][j] = epsi[i][j];

			/**2-1-4释放需要随时释放的临时变量的内存**/
			for (int i = 0; i < n_nei; i++)
				delete[] V[i];
			delete[] V;//(nei1)释放导数矩阵V
			for (int i = 0; i < Dimension; i++)
				delete[] Inverse_VTV[i];
			delete[] Inverse_VTV;//(nei1)释放VTV^-1
			for (int i = 0; i < n_nei; i++)//(nei2)释放临时临近点位移disp_nei
				delete[] disp_nei[i];
			delete[] disp_nei;

			/****2-1根据位移算应变完成****/

			/****2-2应变算应力****/
			/**2-2-1准备临时变量赋值**/
			//(m2)epsi[][]还在
			//(m3)sigm[][]局外变量

			/**2-2-2本构方程计算应力**/
			Stress_from_Strain(sigm, E, nu, G, epsi, Dimension);

			/**2-2-3临时变量应力传给节点**/
			for (int i = 0; i < Dimension; i++)//算完传值给m
				for (int j = 0; j < Dimension; j++)
					m[p].sigma[i][j] = sigm[i][j];
			/**2-2-4需要随时释放的内存**/

			/****2-2应变算应力完成****/
		}

		//2-2-2更新全场内力&合力
		for (int p = 0; p < Num_m; p++)
		{
			/****2-3应力算内力****/
			/**2-3-1临时变量赋值**/
			n_nei = N_Nei[p];//(nei0)临近点个数n_nei				
			V = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)//(nei1)导数矩阵V
				V[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
				for (int j = 0; j < Dimension; j++)
					V[i][j] = Vect[p][i][j];
			Inverse_VTV = new double* [Dimension];
			for (int i = 0; i < Dimension; i++)
				Inverse_VTV[i] = new double[Dimension];
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
					Inverse_VTV[i][j] = Inverse_VectTVect[p][i][j];
			sigma_nei = new double** [n_nei];//(nei3)临时临近点应力
			for (int n = 0; n < n_nei; n++)
				sigma_nei[n] = new double* [Dimension];
			for (int n = 0; n < n_nei; n++)
				for (int i = 0; i < Dimension; i++)
					sigma_nei[n][i] = new double[Dimension];
			for (int n = 0; n < n_nei; n++)
			{
				int glodex = GlobalIndex_Nei[p][n];
				for (int i = 0; i < Dimension; i++)
					for (int j = 0; j < Dimension; j++)
						sigma_nei[n][i][j] = m[glodex].sigma[i][j];
			}
			for (int i = 0; i < Dimension; i++)//(m3)-sigm[][]
				for (int j = 0; j < Dimension; j++)
					sigm[i][j] = m[p].sigma[i][j];
			for (int i = 0; i < Dimension; i++)//(m4)-n[]
				n[i] = m[p].n[i];

			/**2-3-2计算内力**/
			if ((m[p].IfBoundary) == 0)//内部点
			{
				if (IfPrestoreInverseVTV == 0)//不预存VTV^-1
				InForce_from_Stress_in(inforce, V, sigma_nei, sigm, n_nei, Dimension, IfWeightLeastSquares);
				else
					InForce_from_Stress_in(inforce, V, Inverse_VTV, sigma_nei, sigm, n_nei, Dimension);
			}
			else//外部点
			{
				InForce_from_Stress_bdy(inforce, sigm, n, Dimension);
			}
			/**2-3-3释放临时变量**/
			for (int n = 0; n < n_nei; n++)
				for (int i = 0; i < Dimension; i++)
					delete[] sigma_nei[n][i];
			for (int n = 0; n < n_nei; n++)
				delete[] sigma_nei[n];
			delete[] sigma_nei;//(nei3)释放临近点应力				
			for (int i = 0; i < n_nei; i++)
				delete[] V[i];
			delete[] V;//(nei1)释放导数矩阵V
			for (int i = 0; i < Dimension; i++)
				delete[] Inverse_VTV[i];
			delete[] Inverse_VTV;//(nei1)释放VTV^-1
			/****2-3应力算内力完成****/

			/****2-4令余量=内力+外力****/
			for (int i = 0; i < Dimension; i++)
			{
				if (m[p].Fixed[i] == 1)//固定余量为零
					m[p].R[i] = 0;
				else
					m[p].R[i] = inforce[i] + m[p].Force[i];
			}
			/****2-4令余量=内力完成****/
		}

		//2-3最后释放
		delete[] disp_m; //(m1)释放临时节点位移dis_m
		for (int i = 0; i < Dimension; i++)//(m2)释放临时应变epsi[][]
			delete[] epsi[i];
		delete[] epsi;
		for (int i = 0; i < Dimension; i++)
			delete[] sigm[i];
		delete[] sigm;//(m3)释放节点应力
		delete[] n;//(m4)释放边界点外法向
		delete[] inforce;//(m5)释放临时内力inforce[Dim]
		/******2更新节点信息完成******/

		/******2+记录当前云图******/
		if (RecordRelaxation == 1)
		{
			RelaxOut << "zone T=\"Frame1\"\n";//zone header
			RelaxOut << "STRANDID=1\n";
			RelaxOut << "I=" << Num_m << "\n";
			RelaxOut << "ZONETYPE=Ordered\n";
			RelaxOut << "DATAPACKING=POINT\n";
			RelaxOut << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
			RelaxOut << "SOLUTIONTIME=" << it_orien << "\n";
			for (int p = 0; p < Num_m; p++)//data
			{
				RelaxOut << m[p].X[0] << " " << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << " ";
				if (m[p].IfBoundary == 1)
					RelaxOut << m[p].R[0] / ResiToleranceBounadry << " " << m[p].R[1] / ResiToleranceBounadry << "\n";
				else
					RelaxOut << m[p].R[0] / ResiToleranceInternal << " " << m[p].R[1] / ResiToleranceInternal << "\n";
			}
		}
		/******2+记录当前云图完成******/

		/******3检测余量数&无量纲余量分量最大值******/
		SpecificResidual_Max = 0.0;
		for (int p = 0; p < Num_m; p++)
			for (int i = 0; i < Dimension; i++)
			{
				//3-1检测余量数
				if (m[p].IfBoundary == 1)//区分边界点
					ResiTolerance = ResiToleranceBounadry;
				else
					ResiTolerance = ResiToleranceInternal;
				if (abs(m[p].R[i]) > ResiTolerance)
					Residuals++;
				//3-2SpecificResidual_Max
				if (abs(m[p].R[i]) / ResiTolerance > SpecificResidual_Max)
					SpecificResidual_Max = abs(m[p].R[i]) / ResiTolerance;

			}
		cout << "it_orien=" << it_orien << " Residuals=" << Residuals << " Rate_Resi=" << 100.0 * (double)Residuals / ((double)Num_m * (float)Dimension) <<"%" << "\n";
		Itfout << it_orien << " " << Residuals << " " << 100.0 * (double)Residuals / ((double)Num_m * (double)Dimension) <<" " <<SpecificResidual_Max<< "\n";
		/******4对每个Node循环Move******/
		if (Residuals != 0)
		{
			for (int p = 0; p < Num_m; p++)
			{
				for (int i = 0; i < Dimension; i++)
				{
					if (m[p].Fixed[i] == 0)//非固定
					{
						m[p].displacement[i] += MoveRate[p][i] * m[p].R[i];//内部点deltadisplacement alpha_in[p]=(MoveCoefficient/(E * (3.0 - nu) / (1 - nu * nu) / m[p].r / m[p].r / 2))
					}
					else;//固定点
				}
			}
			it_orien++;
		}
		else
			break;
		/******4对每个Node循环Move完成******/
		
	} while (it_orien<= MaxIteration);
	Itfout.close();
	/******2~4主循环完成Main Loop Over******/

	/******5输出结果Output Results******/
	//松弛记录关闭文件
	if (RecordRelaxation == 1)
		RelaxOut.close();
	/****5-1全场结果field results****/
	fout.open(ProjectPath / "results/SimulatedStress.dat");
	fout << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\"ReX\"" << "\"ReY\"" << "\n";
	for (int p = 0; p < Num_m; p++)
	{
		fout << m[p].X[0] << " " << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << " ";
		if(m[p].IfBoundary==1)
			fout << m[p].R[0]/ResiToleranceBounadry << " " << m[p].R[1] / ResiToleranceBounadry << "\n";
		else
			fout << m[p].R[0] / ResiToleranceInternal << " " << m[p].R[1] / ResiToleranceInternal << "\n";
	}
	fout.close();

	/****5-2纵截面Y=0结果****/
	if (OutputLongitudinalSection == 1)
	{
		LongSecfout.open(ProjectPath / "results/LongitudinalSection.dat");
		LongSecfout << "VARIABLES=" << "\"X\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\n";
		for (int p = 0; p < Num_m; p++)
		{
			if (abs(m[p].X[1]) <= 0.02)
				LongSecfout << m[p].X[0] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << "\n";
		}
		LongSecfout.close();
	}
	
	/****5-3横截面X=L/2结果****/
	if (OutputCrossSection == 1)
	{
		CrossSecfout.open(ProjectPath / "results/CrossSection.dat");
		CrossSecfout << "VARIABLES=" << "\"Y\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\n";
		for (int p = 0; p < Num_m; p++)
		{
			if (abs(m[p].X[0] - L / 2) <= 0.02)
				CrossSecfout << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << "\n";
		}
		CrossSecfout.close();
	}
	

	//释放全局变量内存free the memory
	delete[] m;
	
	for (int i = 0; i < Num_m; i++)
		delete[] GlobalIndex_Nei[i];
	delete[] GlobalIndex_Nei;
	//释放V[][][]
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < N_Nei[p]; i++)
			delete[] Vect[p][i];
	for (int p = 0; p < Num_m; p++)
		delete[] Vect[p];
	delete[] Vect;
	//释放预计算VTV的逆
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < Dimension; i++)
			delete[] Inverse_VectTVect[p][i];
	for (int p = 0; p < Num_m; p++)
		delete[] Inverse_VectTVect[p];
	delete[] Inverse_VectTVect;
	//1-5-7释放移动刚度矩阵内存
	for (int p = 0; p < Num_m; p++)
		delete[] K[p];
	delete[] K;
	for (int p = 0; p < Num_m; p++)//释放MoveRate
		delete[] MoveRate[p];
	delete[] MoveRate;

	delete[] N_Nei;

	LogFout.open(ProjectPath / "results/Log.dat");
	//cout << "The run time is:" << (double)clock() / CLOCKS_PER_SEC << "s" << endl;ResiToleranceInternal
	cout << "计算完成！\n";
	LogFout << "MoveCoefficient=" << MoveCoefficient << "\n";
	LogFout << "ResiToleranceInternal=" << ResiToleranceInternal << "\n";
	LogFout << "ResiToleranceBounadry=" << ResiToleranceBounadry << "\n";
	LogFout << "The run time is:" << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
	LogFout.close();
	return 0;
}

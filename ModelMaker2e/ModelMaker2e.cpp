//����ص�������ģ�͵�ģ����������:1)ģ������ 2)��ȡ��֤ 3)���۽����
//��˹̶����Ҷ����
//�̶��ĸ����۽�
#include <iostream>
#include <fstream>
#include "Marker2D.h"

using namespace std;
//����model�ļ�
int main()
{
	/******0-����ǰ׼��******/
	//0-1������
	double L = 2.0, h = 1.0;//����
	double Mesh00X = 0.0, Mesh00Y = -h / 2;//���½�����
	//0-2����ͽڵ�
	int nx = 20, ny = 10;//X21����,Y11����
	double dh = L / nx;
	Num_m = (nx + 1) * (ny + 1);//Marker��Ŀ
	//0-3�غ�
	double LoadM = 2.0e5;//Pa
	/******0-����ǰ׼�����******/


	/******1-Marker����������m[Num_m]******/
	/*****1-1����Marker����*****/
	Marker2D* m = new Marker2D[Num_m];

	/*****1-2Marker����/��Ա��������*****/

	/****1-2-1����****/
	/***1-2-1-0�������׼��***/
	//����׼��
	double** X = new double* [nx + 1];//�����ά������
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
	double* X1 = new double[Num_m];//��ά��һά
	double* Y1 = new double[Num_m];	
	for (int i = 0; i <= nx; i++)
	{
		for (int j = 0; j <= ny; j++)
		{
			X1[j * (nx + 1) + i] = X[i][j];
			Y1[j * (nx + 1) + i] = Y[i][j];
		}
	}	
	for (int i = 0; i <= nx; i++)//�ͷŶ�ά�������ڴ�
	{
		delete[] X[i];
		delete[] Y[i];
	}
	delete[] X;
	delete[] Y;
	/***1-2-1-0�������׼�����***/
	/***1-2-1-1����-����***/
	/***1-2-1-2-1����-�߽��?***/
	/***1-2-1-2-2����-�߽��ⷨ��***/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].X[0] = X1[k];//1-1����
		m[k].X[1] = Y1[k];
		if (X1[k] == L || X1[k] == 0 || (Y1[k] == -h / 2) || (Y1[k] == h / 2))//1-2-1�߽�
			m[k].IfBoundary = 1;
		else
			m[k].IfBoundary = 0;
		//1-2-2�߽��ⷨ��
		if (X1[k] == 0)//��
		{
			m[k].n[0] = -1.0;
			m[k].n[1] = 0;
		}
		else if (X1[k] == L)//��
		{
			m[k].n[0] = 1.0;
			m[k].n[1] = 0;
		}
		else if (Y1[k] == h/2.0)//��
		{
			m[k].n[0] = 0;
			m[k].n[1] = 1.0;
		}
		else if (Y1[k] == -h/2.0)//��
		{
			m[k].n[0] = 0;
			m[k].n[1] = -1.0;
		}
		else//�ڲ���
		{
			m[k].n[0] = 0;
			m[k].n[1] = 0;
		}
	}
	/****1-2-1�������****/

	/****1-2-2��ɢ������****/
	/***1-2-2-0��ɢ���������׼��***/
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
	/***1-2-2-1��ɢ������-�ٽ���뾶***/
	/***1-2-2-2��ɢ������-�ٽ�����***/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].r = r_m[k];//2-1��ɢ������-�ٽ���뾶
		m[k].Nei = N_Nei[k];//2-2��ɢ������-�ٽ�����
	}
	/****1-2-2��ɢ���������****/

	/****1-2-3Լ�����غ�****/
	/***1-2-3-1Լ��***/
	/***1-2-3-2�غ�***/
	for (int k = 0; k < Num_m; k++)
	{
		if (X1[k] == 0)//3-1Լ��
		{
			m[k].Fixed[0] = 1;//��߹̶�
			m[k].Fixed[1] = 1;
		}
		else
		{
			m[k].Fixed[0] = 0;
			m[k].Fixed[1] = 0;
		}
		if (X1[k] == L)//3-2�غ�
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
	/****1-2-3Լ�����غ����****/

	/****1-2-4����****/
	for (int k = 0; k < Num_m; k++)
	{
		m[k].R[0] = 0;//Լ��������Ϊ0
		m[k].R[1] = 0;
	}
	/****1-2-4�������****/

	/****1-2-5���****/
	for (int k = 0; k < Num_m; k++)
	{
		for (int i = 0; i < 2; i++)
		{
			m[k].u[i] = 0;//5-1λ��
			for (int j = 0; j < 2; j++)
			{
				m[k].epsilon[i][j] = 0;//5-2Ӧ��
				m[k].sigma[i][j] = 0;//5-3����
			}
		}
	}
	for (int k = 0; k < Num_m; k++)//���λ��
	{
		if (X1[k] == 0)
			m[k].u[1] = -6 * LoadM * (X1[k] * X1[k] + nu * Y1[k] * Y1[k]) / (E * h * h * h);
	}
	/****1-2-5������****/

	/*****1-2Marker����/��Ա�����������*****/
	/******1-Marker����������m[Num_m]���******/
	

	/******2���������ɢ������******/
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
	/******2���������ɢ���������******/


	/******3�ļ����******/
	/*****3-0�ļ�(��)׼��*****/
	ofstream mfout,bmfout;//�����ļ������
	mfout.open("model2.dat");//��������ļ�
	bmfout.open("model2.b", ios::binary);//���Ӷ���������ļ�

	/*****3-1���Marker������*****/
	mfout << Num_m<<"\n";
	for (int k = 0; k < Num_m; k++)
	{
		mfout << m[k].X[0] << " " << m[k].X[1]<<" ";//1-1��������
		mfout << m[k].IfBoundary << " ";//1-2-1�����Ƿ�߽�
		mfout << m[k].n[0] << " " << m[k].n[1] << " ";//1-2-2���α߽��ⷨ��
		mfout << m[k].r << " ";//2-1��ɢ������-�ٽ���뾶
		mfout << m[k].Nei << " ";//2-2��ɢ������-�ٽ�����
		mfout << m[k].Fixed[0] << " "<< m[k].Fixed[1] << " ";//3-1Լ��-�̶�?
		mfout << m[k].Force[0] << " " << m[k].Force[1] << " ";//3-2�غ�-�غ�
		mfout << m[k].R[0] << " " << m[k].R[1] << " ";//4����
		mfout << m[k].u[0] << " " << m[k].u[1];//5-1���-λ��
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
	//������
	bmfout.write((char*)&Num_m, sizeof Num_m);
	bmfout.write((char*)&(m[0]), (long long)Num_m * sizeof(m[0]));

	/*****3-2connectivity�������*****/
	for (int k = 0; k < Num_m; k++)
		mfout << N_Nei[k] << "\n";
	for (int i = 0; i < Num_m; i++)
	{
		for (int j = 0; j < N_Nei[i]; j++)
			mfout << GlobalIndex_Nei[i][j] << " ";
		mfout << "\n";
	}
	//������
	bmfout.write((char*)N_Nei, (long long)Num_m * sizeof(N_Nei[0]));
	for (int i = 0; i < Num_m; i++)
	{
		bmfout.write((char*)GlobalIndex_Nei[i], (long long)N_Nei[i] * sizeof(GlobalIndex_Nei[i][0]));
	}

	/*****3-3�ر��ļ��ͷ��ڴ�*****/
	mfout.close();
	bmfout.close();
	delete[] m;//�ͷŶ��������ڴ�
	delete[] r_m;
	delete[] N_Nei;
	for (int i = 0; i < Num_m; i++)
		delete[] GlobalIndex_Nei[i];
	delete[] GlobalIndex_Nei;
	/******3�ļ�������******/


	/******4���Զ�ȡ******/
		/*4-0�������ļ�*/
	ifstream mfin;
	mfin.open("model2.dat");
	/*4-1��ȡ���ڴ�*/
	int num_m;//4-1-1Marker����
	mfin >> num_m;
	cout << num_m;
	Marker2D* m2 = new Marker2D[num_m];//4-1-2Marker����
	for (int k = 0; k < num_m; k++)
	{
		mfin >> m2[k].X[0] >> m2[k].X[1];//1-1��������
		mfin >> m2[k].IfBoundary;//1-2-1�����Ƿ�߽�
		mfin >> m2[k].n[0] >> m2[k].n[1];//1-2-2���α߽��ⷨ��
		mfin >> m2[k].r;//2-1��ɢ������-�ٽ���뾶
		mfin >> m2[k].Nei;//2-2��ɢ������-�ٽ�����
		mfin >> m2[k].Fixed[0] >> m2[k].Fixed[1];//3-1Լ��-�̶�?
		mfin >> m2[k].Force[0] >> m2[k].Force[1];//3-2�غ�-�غ�
		mfin >> m2[k].R[0] >> m2[k].R[1];//4����
		mfin >> m2[k].u[0] >> m2[k].u[1];//5-1���-λ��
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
	int* n_nei = new int[num_m];//4-1-3connectivity����
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
	/*4-2�ر������ļ�*/
	mfin.close();

	/*4-3���Զ�ȡ����������ļ�*/
	ofstream mfout2;
	mfout2.open("test.dat");//4-3-0���ļ�
	mfout2 << num_m << "\n";//4-3-1���marker����
	for (int k = 0; k < num_m; k++)
	{
		mfout2 << m2[k].X[0] << " " << m2[k].X[1] << " ";//1-1��������
		mfout2 << m2[k].IfBoundary << " ";//1-2-1�����Ƿ�߽�
		mfout2 << m2[k].n[0] << " " << m2[k].n[1] << " ";//1-2-2���α߽��ⷨ��
		mfout2 << m2[k].r << " ";//2-1��ɢ������-�ٽ���뾶
		mfout2 << m2[k].Nei << " ";//2-2��ɢ������-�ٽ�����
		mfout2 << m2[k].Fixed[0] << " " << m2[k].Fixed[1] << " ";//3-1Լ��-�̶�?
		mfout2 << m2[k].Force[0] << " " << m2[k].Force[1] << " ";//3-2�غ�-�غ�
		mfout2 << m2[k].R[0] << " " << m2[k].R[1] << " ";//4����
		mfout2 << m2[k].u[0] << " " << m2[k].u[1];//5-1���-λ��
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
	for (int k = 0; k < num_m; k++)//4-3-2���connectivity����
		mfout2 << n_nei[k] << "\n";
	for (int i = 0; i < num_m; i++)
	{
		for (int j = 0; j < n_nei[i]; j++)
			mfout2 << globalindex_nei[i][j] << " ";
		mfout2 << "\n";
	}
	mfout2.close();//4-3-3�ر��ļ��ͷ��ڴ�
	delete[] m2;
	delete[] n_nei;
	for (int i = 0; i < num_m; i++)
		delete[] globalindex_nei[i];
	delete[] globalindex_nei;
	/******4���Զ�ȡ���******/

	/******5���۽�******/
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

	delete[] X1;//�ͷ�һά������ʱ�ڴ�
	delete[] Y1;
	delete[] sigmaxx;
	delete[] sigmayy;
	delete[] sigmaxy;
	delete[] displacementx;
	delete[] displacementy;
	/******5���۽����******/

	return 0;
}
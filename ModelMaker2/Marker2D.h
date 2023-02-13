#pragma once
//��ǵ�Marker��
//���ڳ�Ա����
class Marker2D
{
public:
	double X[2];//1-1��ǵ�����
	bool IfBoundary;//1-2-1�ڲ���||�߽��
	double n[2];//1-2-2�߽��ⷨ��
	double r;//2-1�ٽ���뾶
	int Nei;//2-2�ٽ�����
	bool Fixed[2];//3-1�Ƿ�̶�
	double Force[2];//3-2-1�غ�	
	double R[2];//4-1����
	double u[2];//��ǵ�λ��5-1
	double epsilon[2][2];//��ǵ�Ӧ��5-2
	double sigma[2][2];//��ǵ�Ӧ��5-3
};

//����ȫ�ֱ�������
//2-3number of Markers
extern int Num_m;
//2-4GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Marker
extern int** GlobalIndex_Nei;
//2-5Vect[m][Nei][2] is the connectivity matrix of the m-th Marker
extern double*** Vect;
extern double E;//0-1����ģ��
extern double nu;//0-2���ɱ�
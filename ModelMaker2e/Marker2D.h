#pragma once
//标记点Marker类
//类内成员数据
class Marker2D
{
public:
	double X[2];//1-1标记点坐标
	bool IfBoundary;//1-2-1内部点||边界点
	double n[2];//1-2-2边界外法向
	double r;//2-1临近域半径
	int Nei;//2-2临近点数
	bool Fixed[2];//3-1是否固定
	double Force[2];//3-2-1载荷	
	double R[2];//4-1余量
	double u[2];//标记点位移5-1
	double epsilon[2][2];//标记点应变5-2
	double sigma[2][2];//标记点应力5-3
};

//类外全局变量声明
//2-3number of Markers
extern int Num_m;
//2-4GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Marker
extern int** GlobalIndex_Nei;
//2-5Vect[m][Nei][2] is the connectivity matrix of the m-th Marker
extern double*** Vect;
extern double E;//0-1弹性模量
extern double nu;//0-2泊松比
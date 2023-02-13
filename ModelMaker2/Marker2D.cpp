#include "Marker2D.h"
#include <iostream>



//类外全局变量定义
//2-3number of Markers
int Num_m=1;
//2-4GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Marker
int** GlobalIndex_Nei=NULL;
//2-5Vect[m][Nei][2] is the connectivity matrix of the m-th Marker
double*** Vect=NULL;
double E=2.0e11;//0-1弹性模量
double nu=0.3;//0-2泊松比

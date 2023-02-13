#include "ApplicationLayer.h"
#include <cstddef>

string ProjectName;
path ProjectPath, BinaryModelFilePath, ASCIIModelFilePath, ConfigureFilePath;

int Dimension = 2;
Node2D* m=NULL;
//下降率
double alpha_in=1.0e-13;
double alpha_b=1.0e-12;

bool UseBinaryModel = 1;
double MoveCoefficient = 0.5;//移动系数

double ResiTolerance = 100;//余量判断容差##########全局
double ResiToleranceInternal = 1000;//内部节点余量判断容差##########全局
double ResiToleranceBounadry = 60;//内部节点余量判断容差##########全局

bool IfWeightLeastSquares = 0, IfOutputMoveStiffness = 0, RecordRelaxation = 0, OutputLongitudinalSection = 0, OutputCrossSection = 0;

int it_orien = 0;//记录orientation次数##########全局
int MaxIteration = 10000;

/*不需要设置的变量*/
//暂时不可修改
bool IfPrestoreInverseVTV = 0;//是否预存储InverseVTV!!!!!!线性问题节省计算量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//其他遗留的一些变量，仅验证算例时有用
double MoveCoefficient_in = 0.64;//内部节点移动系数0.64InternalMoveCoefficient
double MoveCoefficient_b = 0.8;//0.8boundaryMoveCoefficient
const bool IfUseRoughStiffness = 0;//1使用粗糙移动刚度!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
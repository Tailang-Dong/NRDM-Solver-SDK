#include "ApplicationLayer.h"
#include <cstddef>

string ProjectName;
path ProjectPath, BinaryModelFilePath, ASCIIModelFilePath, ConfigureFilePath;

int Dimension = 2;
Node2D* m=NULL;//pointer to the Node Array
//移动率moving rate
double alpha_in=1.0e-13;
double alpha_b=1.0e-12;

bool UseBinaryModel = 1;//是否使用二进制模型if Use Binary Model
double MoveCoefficient = 0.5;//绉诲姩绯绘暟

double ResiTolerance = 100;//余量判断容差##########全局residual tolerance. global
double ResiToleranceInternal = 1000;//内部节点余量判断容差##########全局internal residual tolerance. global
double ResiToleranceBounadry = 60;//内部节点余量判断容差##########全局boundary residual tolerance. global

bool IfWeightLeastSquares = 0, IfOutputMoveStiffness = 0, RecordRelaxation = 0, OutputLongitudinalSection = 0, OutputCrossSection = 0;//计算设置configuration

int it_orien = 0;//记录orientation次数##########全局record the iteration. global
int MaxIteration = 10000;//防止死循环Max Iterations to forbid endless loop

/*不需要设置的变量 parameters, need not set*/
//鏆傛椂涓嶅彲淇敼
//bool IfPrestoreInverseVTV = 0;//是否预存储InverseVTV!!!!!!线性问题节省计算量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


//鍏朵粬閬楃暀鐨勪竴浜涘彉閲忥紝浠呴獙璇佺畻渚嬫椂鏈夌敤
double MoveCoefficient_in = 0.64;//鍐呴儴鑺傜偣绉诲姩绯绘暟0.64InternalMoveCoefficient
double MoveCoefficient_b = 0.8;//0.8boundaryMoveCoefficient
const bool IfUseRoughStiffness = 0;//1浣跨敤绮楃硻绉诲姩鍒氬害!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
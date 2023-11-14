#pragma once
#include "MechanicsLayer.h"
#include <string>
#include <filesystem>
using namespace std;
using namespace std::filesystem;

extern string ProjectName;
extern path ProjectPath, BinaryModelFilePath, ASCIIModelFilePath, ConfigureFilePath;

extern int Dimension;
extern Node2D* m;//pointer to the Node Array
//移动率moving rate
extern double alpha_in;
extern double alpha_b;

extern bool UseBinaryModel;//是否使用二进制模型if Use Binary Model
extern double MoveCoefficient;//移动系数

extern double ResiTolerance;//余量判断容差##########全局 residual tolerance. global
extern double ResiToleranceInternal;//内部节点余量判断容差##########全局internal residual tolerance. global
extern double ResiToleranceBounadry;//内部节点余量判断容差##########全局boundary residual tolerance. global

extern bool IfWeightLeastSquares, IfOutputMoveStiffness, RecordRelaxation, OutputLongitudinalSection, OutputCrossSection;//计算设置configuration

extern int it_orien;//记录orientation次数##########全局record the iteration. global
extern int MaxIteration;//防止死循环Max Iterations to forbid endless loop

/*不需要设置的变量. parameters, don't need setting*/
//暂时不可修改
//extern bool IfPrestoreInverseVTV;//是否预存储InverseVTV!!!!!!线性问题节省计算量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//其他遗留的一些变量，仅验证算例时有用
extern double MoveCoefficient_in;//内部节点移动系数0.64InternalMoveCoefficient
extern double MoveCoefficient_b;//0.8boundaryMoveCoefficient
extern const bool IfUseRoughStiffness;//1使用粗糙移动刚度!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


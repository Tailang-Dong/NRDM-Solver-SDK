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
//�ƶ���moving rate
extern double alpha_in;
extern double alpha_b;

extern bool UseBinaryModel;//�Ƿ�ʹ�ö�����ģ��if Use Binary Model
extern double MoveCoefficient;//�ƶ�ϵ��

extern double ResiTolerance;//�����ж��ݲ�##########ȫ�� residual tolerance. global
extern double ResiToleranceInternal;//�ڲ��ڵ������ж��ݲ�##########ȫ��internal residual tolerance. global
extern double ResiToleranceBounadry;//�ڲ��ڵ������ж��ݲ�##########ȫ��boundary residual tolerance. global

extern bool IfWeightLeastSquares, IfOutputMoveStiffness, RecordRelaxation, OutputLongitudinalSection, OutputCrossSection;//��������configuration

extern int it_orien;//��¼orientation����##########ȫ��record the iteration. global
extern int MaxIteration;//��ֹ��ѭ��Max Iterations to forbid endless loop

/*����Ҫ���õı���. parameters, don't need setting*/
//��ʱ�����޸�
//extern bool IfPrestoreInverseVTV;//�Ƿ�Ԥ�洢InverseVTV!!!!!!���������ʡ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//����������һЩ����������֤����ʱ����
extern double MoveCoefficient_in;//�ڲ��ڵ��ƶ�ϵ��0.64InternalMoveCoefficient
extern double MoveCoefficient_b;//0.8boundaryMoveCoefficient
extern const bool IfUseRoughStiffness;//1ʹ�ôֲ��ƶ��ն�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


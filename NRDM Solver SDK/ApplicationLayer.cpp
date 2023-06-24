#include "ApplicationLayer.h"
#include <cstddef>

string ProjectName;
path ProjectPath, BinaryModelFilePath, ASCIIModelFilePath, ConfigureFilePath;

int Dimension = 2;
Node2D* m=NULL;//pointer to the Node Array
//�ƶ���moving rate
double alpha_in=1.0e-13;
double alpha_b=1.0e-12;

bool UseBinaryModel = 1;//�Ƿ�ʹ�ö�����ģ��if Use Binary Model
double MoveCoefficient = 0.5;//�ƶ�ϵ��

double ResiTolerance = 100;//�����ж��ݲ�##########ȫ��residual tolerance. global
double ResiToleranceInternal = 1000;//�ڲ��ڵ������ж��ݲ�##########ȫ��internal residual tolerance. global
double ResiToleranceBounadry = 60;//�ڲ��ڵ������ж��ݲ�##########ȫ��boundary residual tolerance. global

bool IfWeightLeastSquares = 0, IfOutputMoveStiffness = 0, RecordRelaxation = 0, OutputLongitudinalSection = 0, OutputCrossSection = 0;//��������configuration

int it_orien = 0;//��¼orientation����##########ȫ��record the iteration. global
int MaxIteration = 10000;//��ֹ��ѭ��Max Iterations to forbid endless loop

/*����Ҫ���õı��� parameters, need not set*/
//��ʱ�����޸�
//bool IfPrestoreInverseVTV = 0;//�Ƿ�Ԥ�洢InverseVTV!!!!!!���������ʡ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//����������һЩ����������֤����ʱ����
double MoveCoefficient_in = 0.64;//�ڲ��ڵ��ƶ�ϵ��0.64InternalMoveCoefficient
double MoveCoefficient_b = 0.8;//0.8boundaryMoveCoefficient
const bool IfUseRoughStiffness = 0;//1ʹ�ôֲ��ƶ��ն�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
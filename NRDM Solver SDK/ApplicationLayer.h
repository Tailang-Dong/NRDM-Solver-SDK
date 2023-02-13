#pragma once
#include "MechanicsLayer.h"
#include <string>
#include <filesystem>
using namespace std;
using namespace std::filesystem;

extern string ProjectName;
extern path ProjectPath, BinaryModelFilePath, ASCIIModelFilePath, ConfigureFilePath;

extern int Dimension;
extern Node2D* m;
//�½���
extern double alpha_in;
extern double alpha_b;

extern bool UseBinaryModel;//�Ƿ�ʹ�ö�����ģ��
extern double MoveCoefficient;//�ƶ�ϵ��

extern double ResiTolerance;//�����ж��ݲ�##########ȫ��
extern double ResiToleranceInternal;//�ڲ��ڵ������ж��ݲ�##########ȫ��
extern double ResiToleranceBounadry;//�ڲ��ڵ������ж��ݲ�##########ȫ��

extern bool IfWeightLeastSquares, IfOutputMoveStiffness, RecordRelaxation, OutputLongitudinalSection, OutputCrossSection;//��������

extern int it_orien;//��¼orientation����##########ȫ��
extern int MaxIteration;//��ֹ��ѭ��

/*����Ҫ���õı���*/
//��ʱ�����޸�
extern bool IfPrestoreInverseVTV;//�Ƿ�Ԥ�洢InverseVTV!!!!!!���������ʡ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//����������һЩ����������֤����ʱ����
extern double MoveCoefficient_in;//�ڲ��ڵ��ƶ�ϵ��0.64InternalMoveCoefficient
extern double MoveCoefficient_b;//0.8boundaryMoveCoefficient
extern const bool IfUseRoughStiffness;//1ʹ�ôֲ��ƶ��ն�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "FileLayer.h"

std::ifstream InputProjectFile;//��Ŀ�ļ�������
std::ifstream InputConfigureFile;//�����ļ�������

std::ifstream mfin;//��ȡģ���ļ�������
std::ofstream fout;//������
std::ofstream LongSecfout;//�ݽ��������Longitudinal section
std::ofstream CrossSecfout;//����������Cross Section
std::ofstream Itfout;//������������¼
std::ofstream RelaxOut;//�ɳڼ�¼

std::ofstream LogFout;//��־��¼

double SpecificResidual_Max=0;//��¼���нڵ�����������������������ֵ
#include "FileLayer.h"

std::ifstream InputProjectFile;//��Ŀ�ļ�������.InputProjectFile ifstream
std::ifstream InputConfigureFile;//�����ļ�������.InputConfigureFile ifstream

std::ifstream mfin;//��ȡģ���ļ�������ifstream to read the model file
std::ofstream fout;//������ofstream to write the result file
std::ofstream LongSecfout;//�ݽ��������ofstream to write the result file. Longitudinal Section
std::ofstream CrossSecfout;//����������ofstream to write the result file. Cross Section
std::ofstream Itfout;//ofstream to write the result file. record the iteration curves.
std::ofstream RelaxOut;//�ɳڼ�¼ofstream to write the result file. record the iterations.

std::ofstream LogFout;//��־��¼

double SpecificResidual_Max=0;//��¼���нڵ�����������������������ֵrecord the R^s_max
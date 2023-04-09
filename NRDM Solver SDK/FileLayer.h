#pragma once
#include <fstream>
//5�ļ���5File Layer
//��ǰ׼���ļ�I/O��. prepare for file I/O
//�ļ��ļ��е�ַ����ȵ�file path, folder path, path objects. etc.

extern std::ifstream InputProjectFile;//��Ŀ�ļ�������.InputProjectFile ifstream
extern std::ifstream InputConfigureFile;//�����ļ�������.InputConfigureFile ifstream

extern std::ifstream mfin;//��ȡģ���ļ�������ifstream to read the model file
extern std::ofstream fout;//������ofstream to write the result file
extern std::ofstream LongSecfout;//�ݽ��������ofstream to write the result file. Longitudinal Section
extern std::ofstream CrossSecfout;//����������ofstream to write the result file. Cross Section
extern std::ofstream Itfout;//������������¼ofstream to write the result file. record the iteration curves.
extern std::ofstream RelaxOut;//�ɳڼ�¼ofstream to write the result file. record the iterations.

extern std::ofstream LogFout;//��־��¼

extern double SpecificResidual_Max;//��¼���нڵ�����������������������ֵrecord the R^s_max


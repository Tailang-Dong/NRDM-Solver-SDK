#pragma once
#include <fstream>
//5�ļ���
//��ǰ׼���ļ�I/O��
//�ļ��ļ��е�ַ����ȵ�

extern std::ifstream InputProjectFile;//��Ŀ�ļ�������
extern std::ifstream InputConfigureFile;//�����ļ�������

extern std::ifstream mfin;//��ȡģ���ļ�������
extern std::ofstream fout;//������
extern std::ofstream LongSecfout;//�ݽ��������Longitudinal Section
extern std::ofstream CrossSecfout;//����������Cross Section
extern std::ofstream Itfout;//������������¼
extern std::ofstream RelaxOut;//�ɳڼ�¼

extern std::ofstream LogFout;//��־��¼

extern double SpecificResidual_Max;//��¼���нڵ�����������������������ֵ


#include "FileLayer.h"

std::ifstream InputProjectFile;//项目文件输入流.InputProjectFile ifstream
std::ifstream InputConfigureFile;//配置文件输入流.InputConfigureFile ifstream

std::ifstream mfin;//读取模型文件输入流ifstream to read the model file
std::ofstream fout;//结果输出ofstream to write the result file
std::ofstream LongSecfout;//纵截面结果输出ofstream to write the result file. Longitudinal Section
std::ofstream CrossSecfout;//横截面结果输出ofstream to write the result file. Cross Section
std::ofstream Itfout;//ofstream to write the result file. record the iteration curves.
std::ofstream RelaxOut;//松弛记录ofstream to write the result file. record the iterations.

std::ofstream LogFout;//日志记录

double SpecificResidual_Max=0;//记录所有节点的无量纲余量各分量的最大值record the R^s_max
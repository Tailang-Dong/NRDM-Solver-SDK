#pragma once
#include <fstream>
//5文件层5File Layer
//提前准备文件I/O流. prepare for file I/O
//文件文件夹地址对象等等file path, folder path, path objects. etc.

extern std::ifstream InputProjectFile;//项目文件输入流.InputProjectFile ifstream
extern std::ifstream InputConfigureFile;//配置文件输入流.InputConfigureFile ifstream

extern std::ifstream mfin;//读取模型文件输入流ifstream to read the model file
extern std::ofstream fout;//结果输出ofstream to write the result file
extern std::ofstream LongSecfout;//纵截面结果输出ofstream to write the result file. Longitudinal Section
extern std::ofstream CrossSecfout;//横截面结果输出ofstream to write the result file. Cross Section
extern std::ofstream Itfout;//迭代至收敛记录ofstream to write the result file. record the iteration curves.
extern std::ofstream RelaxOut;//松弛记录ofstream to write the result file. record the iterations.

extern std::ofstream LogFout;//日志记录

extern double SpecificResidual_Max;//记录所有节点的无量纲余量各分量的最大值record the R^s_max


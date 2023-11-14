#pragma once
//3-1����Node m�����Ա���phi���ݶ�, �����ݶ�Gradphi[]��ָ�����ΪVec[][]���ٽ�������phi_nei[]������ڵ�����phi_m���ٽ�����n_nei��ά��Dim��
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);

//3-2����Node m������ʸ��E���ݶ�, �����ݶ�GradE[][]��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-3������ɢ������Node m������ʸ��E��ɢ��, ����ɢ��DivE��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-4������ɢ������Node m����������T��ɢ��, ����ɢ��GradE��ָ�����ΪVec[][]���ٽ�������T_nei[][][]������ڵ�����T_m[][]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);

/*һ���ع�-��������������Է�������⣬����Ҫ�����ݲ�*/
//3-1����Node m�����Ա���phi���ݶ�, �����ݶ�Gradphi[]��ָ�����ΪVec[][]���ٽ�������phi_nei[]������ڵ�����phi_m���ٽ�����n_nei��ά��Dim��
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim);

//3-2����Node m������ʸ��E���ݶ�, �����ݶ�GradE[][]��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-3������ɢ������Node m������ʸ��E��ɢ��, ����ɢ��DivE��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-4������ɢ������Node m����������T��ɢ��, ����ɢ��GradE��ָ�����ΪVec[][]���ٽ�������T_nei[][][]������ڵ�����T_m[][]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*�����ع�-��������������Է�������⣬����Ҫ�����ݲԤ��VTV����*/
//3-1����Node m�����Ա���phi���ݶ�, �����ݶ�Gradphi[]��ָ�����ΪVec[][]���ٽ�������phi_nei[]������ڵ�����phi_m���ٽ�����n_nei��ά��Dim��
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim);

//3-2����Node m������ʸ��E���ݶ�, �����ݶ�GradE[][]��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-3������ɢ������Node m������ʸ��E��ɢ��, ����ɢ��DivE��ָ�����ΪVec[][]���ٽ�������phi_nei[][]������ڵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);

//3-4������ɢ������Node m����������T��ɢ��, ����ɢ��GradE��ָ�����ΪVec[][]���ٽ�������T_nei[][][]������ڵ�����T_m[][]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*һ-2���ع�-��������������Է�������⣬����Ҫ�����ݲ�-���Ӽ�Ȩ��С����ѡ��*/
//3-1�����ݶ�
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2ʸ��E���ݶ� E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3ʸ��E��ɢ��
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4������ɢ��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);
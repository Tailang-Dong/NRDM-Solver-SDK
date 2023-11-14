#pragma once
/*****Í¨ÓÃÎïÀí²ãÍ·ÎÄ¼þ, º¯ÊýÉùÃ÷. Header of GeneralPhysicsLayer, Founctions Declarations.******/

/***×î³£ÓÃ£¬±»Á¦Ñ§²ãÕæÕýµ÷ÓÃµÄ.  most commonly called, Called by the mechanics layer***/
/*½âÏßÐÔ·½³Ìµ÷ÓÃÇóÄæ·¨, ÊÇ·ñ²ÉÓÃ¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî, Dim=2 or 3. call inverse method to solve matrix equation, with weighted least squares option.*/
//3-1±êÁ¿³¡µÄÌÝ¶ÈËã×Ó. ÌÝ¶ÈGradphi[Dim], Çóµ¼¾ØÕóVec[n_nei][Dim], ÁÙ½üµãÊôÐÔphi_nei[n_nei], ÖÐÐÄµãÊôÐÔphi_m, ÁÙ½üµãÊýn_nei, Î¬¶ÈDim, WeightedLeastSquaresÎª¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî
//3-1gradient of scalar field. gradient=Gradphi[Dim], derivative Matix=Vec[n_nei][Dim], scalar on neighbors=phi_nei[n_nei], scalar on central point=phi_m, number of neighbors=n_nei, Î¬¶ÈDim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2Ê¸Á¿³¡µÄÌÝ¶ÈËã×Ó.ÌÝ¶ÈGradE[Dim][Dim], Çóµ¼¾ØÕóVec[n_nei][Dim], ÁÙ½üµãÊôÐÔE_nei[n_nei][Dim], ÖÐÐÄµãÊôÐÔE_m[Dim], ÁÙ½üµãÊýn_nei, Î¬¶ÈDim, WeightedLeastSquaresÎª¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî
//3-2gradient of vector field. gradient=GradE[Dim][Dim], derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, Î¬¶ÈDim, weighted least squares option=WeightedLeastSquares.
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3Ê¸Á¿³¡µÄÉ¢¶ÈËã×Ó.É¢¶ÈDivE, Çóµ¼¾ØÕóVec[n_nei][Dim], ÁÙ½üµãÊôÐÔE_nei[n_nei][Dim], ÖÐÐÄµãÊôÐÔE_m[Dim], ÁÙ½üµãÊýn_nei, Î¬¶ÈDim, WeightedLeastSquaresÎª¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî
//3-3divergence of vector field. gradient=DivE, derivative Matix=Vec[n_nei][Dim], vector on neighbors=E_nei[n_nei][Dim], vector on central point=E_m[Dim], number of neighbors=n_nei, Î¬¶ÈDim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4ÕÅÁ¿³¡µÄÉ¢¶ÈËã×Ó.É¢¶ÈGradT[Dim], Çóµ¼¾ØÕóVec[n_nei][Dim], ÁÙ½üµãÊôÐÔT_nei[n_nei][Dim][Dim], ÖÐÐÄµãÊôÐÔT_m[Dim][Dim], ÁÙ½üµãÊýn_nei, Î¬¶ÈDim, WeightedLeastSquaresÎª¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî
//3-4divergence of tensor field. gradient=GradT[Dim], derivative Matix=Vec[n_nei][Dim], tensor on neighbors=T_nei[n_nei][Dim][Dim], tensor on central point=T_m[Dim][Dim], number of neighbors=n_nei, Î¬¶ÈDim, weighted least squares option=WeightedLeastSquares.
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);

/***ÖØÔØIII, ²»±»µ÷ÓÃÔçÆÚ°æ±¾. Overload III, early version, not be called***/
/*½âÏßÐÔ·½³Ìµ÷ÓÃµü´ú·¨, ÐèÒªÔ¤Éè¾«¶È. call iterative method to solve matrix equation, need preset tolerance*/
//3-1-III
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);
//3-2-III
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-3-III
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);
//3-4-III
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);

/***ÖØÔØI,²»±»µ÷ÓÃÔçÆÚ°æ±¾. Overload I, early version, not be called***/
/*½âÏßÐÔ·½³Ìµ÷ÓÃÇóÄæ·¨, ÎÞ¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî, Dim=2 or 3. call inverse method to solve matrix equation, No weighted least squares option*/
//3-1-I
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-I
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-I
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-I
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim);

/***ÖØÔØII,ÒÑ²»±»µ÷ÓÃµÄÔçÆÚ°æ±¾. Overload II, early version, not be called***/
/*½âÏßÐÔ·½³Ìµ÷ÓÃÇóÄæ·¨, ÎÞ¼ÓÈ¨×îÐ¡¶þ³ËÑ¡Ïî,Dim=2 or 3.  ÐèÒªInverse_VecTVec[Dim][Dim]
call inverse method to solve matrix equation, No weighted least squares option, need Inverse_VecTVec[Dim][Dim]*/
//3-1-II
void Gradient_of_Scalar(double* Gradphi, double** Vec, double** Inverse_VecTVec, double* phi_nei, double phi_m, int n_nei, int Dim);
//3-2-II
void Gradient_of_Vector(double** GradE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-3-II
void Divergence_of_Vector(double& DivE, double** Vec, double** Inverse_VecTVec, double** E_nei, double* E_m, int n_nei, int Dim);
//3-4-II
void Divergence_of_Tensor(double* GradT, double** Vec, double** Inverse_VecTVec, double*** T_nei, double** T_m, int n_nei, int Dim);

/*ä¸€-2æ¬¡é‡æž?æ–¹é˜µæ±‚é€†è¿›è¡Œçº¿æ€§æ–¹ç¨‹ç»„æ±‚è§£ï¼Œä¸éœ€è¦è¿­ä»£å®¹å·?å¢žåŠ åŠ æƒæœ€å°äºŒä¹˜é€‰é¡¹*/
//3-1æ ‡é‡æ¢¯åº¦
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-2çŸ¢é‡Eçš„æ¢¯åº?E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-3çŸ¢é‡Eçš„æ•£åº?
void Divergence_of_Vector(double& DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, bool WeightedLeastSquares);

//3-4å¼ é‡çš„æ•£åº?
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, bool WeightedLeastSquares);
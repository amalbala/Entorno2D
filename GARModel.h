//*--------------------------------------------------*/
/*					UCLM - ISA						*/
/*				E.T.S.I. INDUSTRIALES				*/
/*													*/
/*	Avdn. Camilo Jose Cela s/n - 13071 Ciudad Real	*/
/*													*/
/*	Proyecto: 3D Target Tool 						*/
/*	Investigador Principal:							*/
/*													*/
/*	Fecha: Diciembre -2004							*/
/*--------------------------------------------------*/



 #ifndef _GAR_MODEL_H_
 #define _GAR_MODEL_H_

 #include "imagutils.h"
 #include "listagen.h"




 
 void InitGARModel(unsigned char* imagen_reg,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,int bInterpol, 
  float weight, int RegionFunc, float ValorReg, int OnlyInit);
 void InitGAMRModel(unsigned char* imagen_reg,int mod,float sigma,int NARROW,float EPSILON,
				  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,int bInterpol, 
				  float weight,int RegionFunc, float ValorReg, int OnlyInit);
 void ReInitGARModel();
 void NextIterGARModel(float inct);
 void NextIterGAMRModel(float inct);
 bool StopCriteriorGARModel();
 void PresentGRAModel(unsigned char *ls);
 void DestroyGARModel();
 void DrawGARModel();
 void FilterBubbles();
 bool EndNarrow();
 bool EndNarrowDist();
 void RecoverCurveGARModel();
 void StoreGRAModel(char* nameFile);
 void SetUpCurveGARModel(lista l);
 int NumLS(unsigned char* imag);

 #endif
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



 #ifndef _LEVEL_SET_H_
 #define _LEVEL_SET_H_

 #include "imagutils.h"
 #include "listagen.h"

//estructura de lista de levelset y regiones contenidas


 
 void inicializar_level_set(unsigned char* imag_reg,float increm,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,bool fuz,
  int num_cl,int cl_int,bool interp);
 
 void inicializar_level_setGVF(unsigned char* imag_reg,float increm,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,bool fuz,
  int num_cl,int cl_int,bool interp,float aporte);
 
 void inicializar_ls_ffm(unsigned char* imag,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral,double th_ki);
 
 void reinicializa_ls();
 void reinicializa_lsGVF();
 
 void almacenar_double(double *mat,char *nombre);
 
 void almacenar_int(int *mat, char* nombre);
 void almacenar_uchar(unsigned char *mat,char *nombre);
 void siguiente_LevelSet(double **MatrizEnerg);
 void siguiente_Geodesica(double **MatrizEnerg);
 void siguiente_Geodesica_GVF(double **MatrizEnerg);
 bool parada(double *kiext);
 bool parada_fmm();
 void actual_level_set(unsigned char *ls);
 void liberar_ls();
 void pintarLevelSet();
 void filtrado_burb();
 bool fin_narrow();
 bool fin_narrow_dist();
 void extraerLeveSet();
 void almacenarLevelSet(char* nameFile);
 void actualizar_selec_ls(lista l);

 #endif
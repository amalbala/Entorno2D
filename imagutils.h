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




#ifndef IMAGUTILS_H
#define IMAGUTILS_H

#define calcSingleSubscript(i,j,w) (((i) * (w)) + (j))
#define calcSingleSubscript3D(i,j,k,h,w) (((h)*(w)*(k)) + ((i) * (w)) + (j))
#define CABECERA_UCLM "(c)  Universidad de Castilla-La Mancha, 2005. Group UCLM-ISA, contact: gloria.bueno@uclm.es\n"

#include <math.h>
#include "listagen.h"
#include "listalista.h"
#define Pi 3.1416

//Definimos el numero de descriptores a utilizar cuanto mayor
//mas precision en la comparacion pero mas memoria
#define N_DESCRIP 124


//estructura de punto 2D
typedef struct{
	    int x; 
        int y;
} punto2D;

typedef struct{
	    float r; 
        float a;
} polar2D;

typedef struct{
	    double x; 
        double y;
} puntof2D;


//Estructura numero complejo (descriptor Fourier)
typedef struct{
        double re; 
        double im;
} descriptores_fourier;

typedef struct StrImage{
	int width;
	int height;
	unsigned char* matrix;
} TImagen;

typedef struct StrMatrizD{
	int width;
	int height;
	double* matrix;
} TMatrizD;


typedef struct StrMatrizI{
	int width;
	int height;
	int* matrix;
} TMatrizI;

//Parte 2D
int cargarImagenPGM(unsigned char **imag, char arch[40]);
void guardarImagenPGM(unsigned char *imag, char arch[40]);
int cargarImagenRAW(unsigned char **imag, char arch[40],int width, int height);
void guardarImagenRAW(unsigned char *imag, char arch[40]);
void guardarImagenMRI(unsigned char *imag, char arch[20]);
void guardarImagenTXT(unsigned char *imag, char arch[20]);
void guardarImagenTXT(double *imag, char arch[20]);
void guardarContornoCNT(lista l, char arch[20]);

//localiza una posicion en un vector unidimensional una posicion bidimiensional
//int calcSingleSubscript(int i, int j, int h);
//Halla la imagen complementaria de una imagen binaria
void imagen_complementaria(double *imag,double* image_comp, int tama);
void Set_X_Y( int w, int h);
//Devuelve las dimensiones de la imagen con la que trabajaremos
void Get_X_Y( int* w, int* h);
//Devuelve las dimensiones de la imagen con la que trabajaremos
void Get_X_Y_Z( int* x, int* y,int *z);
//Gradiente por Prewitt
void grad_diffin(double *imag,double *grad);
void grad_diffin(double *imag,double *gradx, double *grady,double *grad);
void grad_diffin_norm(double *imag,double *gradx, double *grady,double *grad);
//Reescalado a 255
void doubleTounsignedchar(double *imag, unsigned char *imag2);
//Negativizado
void negativo(unsigned char *imag, unsigned char *imagres);
//Gradiente convolucionado con la gaussiana
void grad_conv_gauss(unsigned char* imag, float sigma, float *max,float *min, double* total_gradient);
//Calculo la gaussiana de la imagen
void gaussiana(float sigma, unsigned char *imagen, double *imagsuav);
void mediana(unsigned char *imag, unsigned char* imagsuav);
//Calculo de la curvatura de una imagen
void curvatura(double *matriz,double *curv,void gradiente(double * m, double * gx, double* gy, double* g));
void curv_div(double *matriz,double *curv);
//Calculo las distancias aplicando la mascara de Chamfer
void chamfer_distance4v(double *im_input,double *im_output);
void chamfer_distance8v(double *im_input,double *im_output);
void vecin4ssed_distance(double *im_input, punto2D *im_dist);
//Determinar por recursion 4 conexa los puntos pertenecientes a una region
void evol_punto(int *imagen_reg,int x, int y);
void evol_punto(double *imagen_reg,int x, int y, int z); 
bool determinar_pto_interno(int *imag, punto2D *pto);
//Funcion que cierra los puntos de una curva usando interpolacion lineal
void cierre_lineal(lista ptos,int numptos, lista curva_cerrada);
//Funcion que cierra los puntos de una curva usando interpolacion por splines
void cierre_spline(lista ptos,int numptos, lista curva_cerrada);
//Funcion que genera ruido aleatorio
void ruido_aleatorio(unsigned char* imag, unsigned char *imag_res);
double curvatura(punto2D pto, punto2D pto_ant, punto2D pto_sig);
void crear_circunferencia_centrada(int orig_x, int orig_y, float radio,int muestreo,lista puntos);
void crear_circunferencia_cerrada(int orig_x, int orig_y, float radio,lista curva);
void multiples_circulos(int num_circ,lista l);
void centro_circulos(lista l);
void extraer_region(unsigned char* imag,lista puntos,unsigned char *result);
//Fuzzy
int etiq_region(int* imagen,int *m_etiq);
int etiq_mayoritaria (lista contorno,int *matriz_etiq);
void fuzzyCmeans(unsigned char *imag,int dim_x, int dim_y,unsigned char *cluster_result,
    int numclusters,int epsilon);
void fuzzyCmeansforzado(unsigned char *imag,int dim_x, int dim_y,
	unsigned char *cluster_result, unsigned char *centroides,int numclusters,int epsilon);
void cargar_centros_fuzzy(unsigned char** vec, char *nombre,int *num_c);
void binarizaImag(unsigned char *imag, unsigned char umb);
unsigned char umbralOtsu(int histo[256]);
void laplaciano(unsigned char *imag,double *laplac);
void histograma(unsigned char* imagen, int histo[256]);

float entropia(unsigned char* imagen);

void crecRegion(double * pdIm,punto2D p2DPto,lista lListaC, lista lListaP);
int extraerSegmentacion(unsigned char *imag, TLlista lGeneral, int ValContour);

void CartesianasToPolares2D(lista lCartesianas,lista lPolares, punto2D p2DCentroGravedad);

void CodigoCadena(unsigned char *imag, int iIntensidad, lista l);
void SelectorPuntos(lista l1, lista l2, int numPuntos);

void DoubleMtxToTXT(char *FileName,double* DoubleMtx);

void GradientVF(unsigned char* iImage, puntof2D* mGVF,double alpha);
void GradientVF_Paragios(unsigned char* iImage, puntof2D* mGVF,double alpha);

void DrawVectImag(puntof2D* ImgVectMtx, unsigned char* Img);
void DrawVectImag(punto2D* ImgVectMtx, unsigned char* Img);
void StoreVectImg(puntof2D* P2dfMtx,char *name);
void StoreVectImg(punto2D* P2dfMtx,char *name);
void VectorModule(puntof2D* P2dfMtx, double *ModuleMtx);
void VectorModule(punto2D* P2dfMtx, double *ModuleMtx);
void ScaleTransform(double* iImageOri,double intervA,double intervB,double* iImageFin);
void VectorNormalized(puntof2D* P2dfMtxIni,puntof2D* P2dfMtxFin);
void VectorNormalized(double* MtxIni,double* MtxFin);

TMatrizI* CreateMatrizI(TMatrizI *Matrix);
TMatrizI* CreateMatrizI(int height, int width);
TMatrizI* CreateMatrizI(TMatrizD *Matrix);
TMatrizD* CreateMatrizD(TMatrizD *Matrix);
TMatrizD* CreateMatrizD(int height, int width);
TImagen* CreateImagen(int height, int width);
TImagen* CreateImagen(TImagen *Imag);

void DestroyMatrix(TMatrizI **Matrix);
void DestroyMatrix(TMatrizD **Matrix);
void DestroyImagen(TImagen **Imag);

int numPicos(lista lPicos, int vecindad, double umbral, unsigned char* imag);
void RegionesInteres(lista lPicos, int vecindad, double umbral, unsigned char* imag);

#endif

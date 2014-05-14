


//Librerias utilizadas
#include <GL/glut.h>
#include <GL/glui.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include "comun.h"
#include "imagutils.h"
#include "listagen.h"
#include "levelset.h"
#include "GARModel.h"
#include "snakes.h"
#include "error.h"

//Si vas a compilar en Windows, descomenta esta linea, solo afecta 
//al caracter de separacion entre directorio \ en Win32 y / en Linux
//#define WIN_32



using namespace std;

#define MAX_ITER 50

#ifndef WIN_32
	//#define PATH_IMAG "Imagenes/dapelvisW/dapelvisW2/\0" 	
	#define PATH_IMAG "Imagenes/\0" 	
#else
	#define PATH_IMAG "Imagenes\dapelvisV1\\0"
	
#endif

#define PATH_SELEC "Inic\\"

//Constantes para los menus
enum {SEGMEN,CARGAR,CARGAR_MRI,MRI_X,MRI_Y,MRI_Z,SUAVIZAR,BINAR,BINAR_UMB,GRADIEN,LAPLAC,NEGAT,LEVELSET,GEODES,SNAKES,PARARSEG,
      LS_EXPAND,LS_COMP,GEOD_EXPAND,GEOD_COMP,GAUSS_GRAD,NUEVA_SELEC,GUARDAR,ENSUCIAR,CLUSTERING, DISTAN,
      GUARDAR_SEL_GEOD,GUARDAR_SEL_LS,GUARDAR_SEL_SN,CIRCULOS,CURVAT,FUZZYGEODES,FUZZYLEVELSET,
      FUZZY,FUZZY_FORZ,FUZZYSNAKES,ENTROPIA,DIB_CIRCULO,GUARDA_SELEC,CARGA_SELEC,GARMODEL,GVF,GVFMOD,GEODESGVF,
		GUARDAR_AUTO,GVFPARAG,GVFMODPARAG,PROCESAR_SEC,PROC_SEC,GAMRMODEL};

enum{MTX_CURV,MTX_FI,MTX_KI,MTX_KIEXT, MTX_GRAD};

enum{GAR_BT, GAR_RT};

#define CAB_PARAM_GAR "Dirección: %d\nSigma: %f\nNarrow Band: %d\nEpsilon: %f\nThressHold: %f\nAporte: %f\nFuncion: %d\nValor Region: %f\nSolo Al Inicio: %d\n"
//#define INI_SEC "Inic/iniP1_s21.txt"
//#define INI_SEC "Inic/initial1.txt"
//#define INI_SEC "Inic/ini1679P2a.txt"
#define INI_SEC "Inic/ini1ap2_32"
//#define INI_SEC "Inic/ini1679P3a.txt"
//#define INI_SEC "Inic/iniAtlasProstate4a.txt"
//#define INI_SEC "Inic/iniAtlasProstate5a.txt"
//#define INI_SEC "Inic/iniAtlasProstate6a.txt"
//#define INI_SEC "Inic/iniDapelvisw2_6.txt"
//#define INI_SEC "Inic/inidapelvis2.txt"
//#define INI_SEC "Inic/IniAtlasProstate1-30sec.txt"
//#define INI_SEC "Inic/Ini29AtlasProstate2Sec.txt"
//#define INI_SEC "Inic/iniBC-0000.txt"
//#define INI_SEC "Inic/BC-0000.txt"
//#define INI_SEC "Inic/iniGAMR3.txt"

int extension;
int FuncReg = 1;
int direccion;
int modoExtraer;
int iFondo;

int dib_circulo = 0;
int bInterp = 0;
int bOnlyInitGAR = 1;
int iMuestreoCirc = 8;
//angulo de la camara
GLint angulo = 0;

//Dimensiones de la ventana de vision
GLsizei weight = 800;
GLsizei height = 800;
//Deplazamiento de centrado de la imagen
GLfloat desplaz_x;
GLfloat desplaz_y;
//Seleccion del menu
GLint opcion = 10000;
GLint subopcion = 10000;

GLint opcion_ls = 10000;
GLint opcion_MRI = 10000;

//Posiciones relativas del raton
GLfloat x_m = 0.0;
GLfloat y_m = 0.0;

int winIdMain;
GLUI* glui;
char label[100];

char name[sizeof(GLUI_String)];
char name_selec[sizeof(GLUI_String)];
char name_wind[sizeof(GLUI_String)];
char name_secuen[sizeof(GLUI_String)];
char name_actual[40];

char arch[40];
int i_imag,j_imag,z_imag;
static unsigned char *imag=NULL;
double *MatrizEnerg = NULL;
int iterac=0;
bool seg = false;
lista puntosselec;
int mov_x=0;
int mov_y=0;
int mov_z=0;
//Parametros de las Snakes
float alpha_sn = 1.0;
float beta_sn = 1.0;
float teta_sn = 1.0;
float w_term = 1.0;
float w_line = 1.0;
float w_edge = 1.0;
float w_press = 1.0;
//Parametros del LevelSet
float SIGMA_GAUSS = 1.0;
float THRESS = 0.0001;
float THRESS_K = 0.0;
int THRESS_ITER = 0;
int NARROWBAND = 10;
int NUM_CIRC = 6;
int iRadio = 10;
float INCREMENTO = 0.3;
float APORTE = 0.4;
float VALREG = 0;
int num_clusters=7;
bool fuzzy_seg = false;
int clust_int=3;
int iterac_actual;
double tiempo;


float VELOC_EPSILON = 0.001;
bool CIRC = false;

int Dim_RAW_Y=512;
int Dim_RAW_X=512;

int Dim_SEC_X=512;
int Dim_SEC_Y=512;
int Dim_SEC_Z=36;

int Inf_SEC = 32;
int Sup_SEC = 45;

//Tiempo Computacional
clock_t tiempo_total;
clock_t tiempo1, tiempo2;
int err;
char mensaje[100];
GLUI_StaticText *stEtiquet;
GLUI_StaticText *stTiempo;
GLUI_StaticText *stIterac;
GLUI_StaticText *stEntropia;

unsigned char* imag3D;
FILE *textresult;
FILE *timeresult;
int n;


/*--------------------------------------------------*/
/* FUNCIONES DE DIBUJADO DE OPENGL					*/
/*--------------------------------------------------*/
/*--------Tralacion de pixeles ---------------------*/
void swappixelOGL(unsigned char *image){
  unsigned char *imag_tras;
  imag_tras = (unsigned char*)calloc(sizeof(unsigned char),i_imag*j_imag);

  int k=0;
  int i,j;
  for(i=i_imag; i>0; i--)
    for(j=0; j<j_imag; j++){
      imag_tras[k] = image[calcSingleSubscript(i-1,j,j_imag)];
      k++;
    }
  for(i=0; i<i_imag*j_imag;i++)
    image[i] = imag_tras[i];
  free(imag_tras);
}

void guardarPuntosImag(char arch[40]){
  ofstream desc;
  desc.open(arch);
  punto2D pto;
  
  for(int i=0; i<i_imag; i++){
	  for(int j=0; j<j_imag; j++){
		  if(imag[calcSingleSubscript(i,j,j_imag)] ==255){
			  desc<<i<<" "<<j<<"\n";
		  }
	  }
  }  
  desc.close();
	
}
/*---------------- SUB MENU -------------------------------*/
void guardarPuntos(lista list,char arch[40]){
  ofstream desc;
  desc.open(arch);
  punto2D pto;
  l_empieza(list);
  while(l_quedan(list)){
    l_dame(list,&pto);
    desc<<pto.x<<" "<<pto.y<<"\n";
  }
  desc.close();
}
/*---------------- SUB MENU -------------------------------*/
void cargarPuntos(char arch[40],lista list){
  FILE *archivo;
  archivo = fopen(arch,"r");
  if(!archivo){
	  error(err,mensaje);
	  stEtiquet->set_text(mensaje);
  }else{
	punto2D pto;
	while(fscanf(archivo,"%d %d\n",&pto.x,&pto.y)!= EOF){
		l_meted(list,&pto);
	}
	fclose(archivo);
  }
}
/*------------------------------------------*/
void cargar_parametros(char *fich_param){

	FILE *archivo;
	archivo = fopen(fich_param,"r");
	fscanf(archivo,CAB_PARAM_GAR,&direccion,&SIGMA_GAUSS,&NARROWBAND,&VELOC_EPSILON,
		&THRESS,&APORTE,&FuncReg,&VALREG,&bOnlyInitGAR);
	fclose(archivo);

}

/*------------------- Bucle de LevelSet -----------------------*/
void bucle (){
 switch(subopcion){
 case LEVELSET: case GEODES: case GEODESGVF: 

	 switch(subopcion){
	  case LEVELSET:
		  siguiente_LevelSet(&MatrizEnerg);
		  break;
	  case GEODES:
		  siguiente_Geodesica(&MatrizEnerg);
		  break;
	  case GEODESGVF:
		  siguiente_Geodesica_GVF(&MatrizEnerg);
	  }

      tiempo2 = clock();
      tiempo_total += tiempo2-tiempo1;
      iterac++;
	  actual_level_set(imag);

	if((!parada(MatrizEnerg))/*parada_fmm()*/&&(!((THRESS_ITER > 0)&&(iterac>THRESS_ITER)))){
      tiempo1 = clock();
      if(/*iterac%NARROWBAND==0*/ fin_narrow()){
		  if((subopcion == LEVELSET)||(subopcion == GEODES)){
			  reinicializa_ls();
		  }
		  if(subopcion == GEODESGVF){
			  reinicializa_lsGVF();
			  }
      }
 
	  char imPGM[40];
	  if(iterac>10){
		  if(iterac>100){
			sprintf(imPGM,"Iteracion_%d.pgm",iterac);
		  }else{
			  sprintf(imPGM,"Iteracion_0%d.pgm",iterac);
		  }
	  }else{
		  sprintf(imPGM,"Iteracion_00%d.pgm",iterac);
	  }
	  int a,b;
	  Get_X_Y(&a,&b);
	  Set_X_Y(j_imag,i_imag);
	  swappixelOGL(imag);
	  guardarImagenPGM(imag,imPGM);
	  swappixelOGL(imag);
	  Set_X_Y(a,b);	  

	  	
    }else{
      printf("Fin de la segmentacion.\n");
      stEtiquet->set_text("Fin de la segmentacion.");
      sprintf(mensaje,"Iteraciones: %i",iterac);
      stIterac->set_text(mensaje);
      sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
      stTiempo->set_text(mensaje);
      printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
      seg=false;
      tiempo_total = 0;

      actual_level_set(imag);
	  if(!l_vacia(puntosselec)){
		  l_dest(&puntosselec);
		  puntosselec = l_nuev(sizeof(punto2D));
	  }
	  actualizar_selec_ls(puntosselec);
	  

	  almacenarLevelSet("ptosResLevelSet.txt");
	  pintarLevelSet();	  
      liberar_ls();
      iterac_actual = iterac;
      iterac=0;
    }
    break;
    case SNAKES:
    
		if(!parada_snk()||(iterac==1)){
        tiempo1 = clock();
		if(fuzzy_seg)
			siguienteSnakeFuzzy();
		else
			siguienteSnake();
        iterac++;
//        printf("%d\n",iterac);
        tiempo2 = clock();
        tiempo_total += tiempo2-tiempo1;
        pintaguardaSnake(imag);
        guardaSnakeInterp();
        swappixelOGL(imag);
        char nombre[20];
        sprintf(nombre,"ResultSnake_%d",iterac);
	      strcat(nombre,".pgm");
        guardarImagenPGM(imag,nombre);
        swappixelOGL(imag);

      }else{
		seg = false;
       printf("Fin de la segmentacion.\n");
       stEtiquet->set_text("Fin de la segmentacion.");
       sprintf(mensaje,"Iteraciones: %i",iterac);
       stIterac->set_text(mensaje);
       sprintf(mensaje,"Tiempo: %f.3 seg.",(double)tiempo_total/CLOCKS_PER_SEC);
       stTiempo->set_text(mensaje);
       printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
       tiempo_total = 0;
       pintaguardaSnakeInterp(imag);
     guardarImagenPGM(imag,"Final_sn.pgm");
     guardaSnakeInterp();
       libera_snk();
	iterac_actual = iterac;
      iterac=0;
      }

    break;
	case GAMRMODEL:
		if((!StopCriteriorGARModel())&&(!((THRESS_ITER > 0)&&(iterac>THRESS_ITER)))){
			tiempo1 = clock();
			if(/*iterac%NARROWBAND==0*/ EndNarrowDist()){
				ReInitGARModel();
			}

      NextIterGAMRModel(INCREMENTO);

      tiempo2 = clock();
      tiempo_total += tiempo2-tiempo1;
      iterac++;
	  PresentGRAModel(imag);
	  /*
	  char imPGM[40];
	  if(iterac>10){
		  if(iterac>100){
			sprintf(imPGM,"Iteracion_%d.pgm",iterac);
		  }else{
			  sprintf(imPGM,"Iteracion_0%d.pgm",iterac);
		  }
	  }else{
		  sprintf(imPGM,"Iteracion_00%d.pgm",iterac);
	  }
	  int a,b;
	  Get_X_Y(&a,&b);
	  Set_X_Y(j_imag,i_imag);
	  swappixelOGL(imag);
	  guardarImagenPGM(imag,imPGM);
	  swappixelOGL(imag);
	  Set_X_Y(a,b);
	  */
	
	  
    }else{
      printf("Fin de la segmentacion.\n");
      stEtiquet->set_text("Fin de la segmentacion.");
      sprintf(mensaje,"Iteraciones: %i",iterac);
      stIterac->set_text(mensaje);
      sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
      stTiempo->set_text(mensaje);
      printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
      seg=false;
	  tiempo = (double)tiempo_total/CLOCKS_PER_SEC;
      tiempo_total = 0;
	  PresentGRAModel(imag);
	  if(!l_vacia(puntosselec)){
		  l_dest(&puntosselec);
		  puntosselec = l_nuev(sizeof(punto2D));
	  }
	  SetUpCurveGARModel(puntosselec);
	  Set_X_Y(j_imag,i_imag);
	  //printf("Numero de LS: %d\n",NumLS(imag));
	  //DrawGARModel();
      DestroyGARModel();
	  iterac_actual = iterac;
      iterac=0;
    }
    break;
case GARMODEL:
		if((!StopCriteriorGARModel())&&(!((THRESS_ITER > 0)&&(iterac>THRESS_ITER)))){
			tiempo1 = clock();
			if(/*iterac%NARROWBAND==0*/ EndNarrowDist()){
				
				ReInitGARModel();
			}

      NextIterGARModel(INCREMENTO);

      tiempo2 = clock();
      tiempo_total += tiempo2-tiempo1;
      iterac++;
	  PresentGRAModel(imag);
	  /*
	  char imPGM[40];
	  if(iterac>10){
		  if(iterac>100){
			sprintf(imPGM,"Iteracion_%d.pgm",iterac);
		  }else{
			  sprintf(imPGM,"Iteracion_0%d.pgm",iterac);
		  }
	  }else{
		  sprintf(imPGM,"Iteracion_00%d.pgm",iterac);
	  }
	  int a,b;
	  Get_X_Y(&a,&b);
	  Set_X_Y(j_imag,i_imag);
	  swappixelOGL(imag);
	  guardarImagenPGM(imag,imPGM);
	  swappixelOGL(imag);
	  Set_X_Y(a,b);
	  */
	
	  
    }else{
      printf("Fin de la segmentacion.\n");
      stEtiquet->set_text("Fin de la segmentacion.");
      sprintf(mensaje,"Iteraciones: %i",iterac);
      stIterac->set_text(mensaje);
      sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
      stTiempo->set_text(mensaje);
      printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
      seg=false;
	  tiempo = (double)tiempo_total/CLOCKS_PER_SEC;
      tiempo_total = 0;
	  PresentGRAModel(imag);
	  if(!l_vacia(puntosselec)){
		  l_dest(&puntosselec);
		  puntosselec = l_nuev(sizeof(punto2D));
	  }
	  SetUpCurveGARModel(puntosselec);
	  Set_X_Y(j_imag,i_imag);
	  //printf("Numero de LS: %d\n",NumLS(imag));
	  //DrawGARModel();
      DestroyGARModel();
	  iterac_actual = iterac;
      iterac=0;
    }
    break;

	case PROC_SEC:
		if(n<Sup_SEC){
			char nameres[40];
			if((!StopCriteriorGARModel())&&(!((THRESS_ITER > 0)&&(iterac>THRESS_ITER)))){
				tiempo1 = clock();
				if(/*iterac%NARROWBAND==0*/ EndNarrowDist()){
					ReInitGARModel();
				}
				//NextIterGARModel(INCREMENTO);
				NextIterGAMRModel(INCREMENTO);
				tiempo2 = clock();
				tiempo_total += tiempo2-tiempo1;
				iterac++;					
				PresentGRAModel(imag);	
			}else{
				printf("Fin del procesamiento de la secuencia.\n");
				stEtiquet->set_text("Secuencia Completada.");
				sprintf(mensaje,"Iteraciones: %i",iterac);
				stIterac->set_text(mensaje);
				sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
				stTiempo->set_text(mensaje);
				printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
				tiempo = (double)tiempo_total/CLOCKS_PER_SEC;
				tiempo_total = 0;				
				PresentGRAModel(imag);
				iterac_actual = iterac;
				//sprintf(nameres,"ResultadoGAR_Imagen%d.pgm\0",n);
				sprintf(nameres,"ResGAR-Corte%d-%.2f-%.2f-%i-%.2f-%.4f-%i-%.3f.pgm",n,APORTE,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,tiempo);
				Set_X_Y(j_imag,i_imag);
				swappixelOGL(imag);
				guardarImagenPGM(imag,nameres);
				swappixelOGL(imag);
				
				if(!l_vacia(puntosselec)){
					l_dest(&puntosselec);
					puntosselec = l_nuev(sizeof(punto2D));
				}
			
				SetUpCurveGARModel(puntosselec);
				l_empieza(puntosselec);
				while(l_quedan(puntosselec)){
					punto2D pto;
					l_dame(puntosselec,&pto);
					fprintf(textresult,"%d %d %d\n",pto.x, pto.y, n);
				}
				fprintf(timeresult,"Corte %d\t Iteraciones: %d\t Tiempo: %.3f seg.\n",n,iterac_actual,tiempo);
				DestroyGARModel();
				
				iterac=0;
				int i,j,l;
				if(n+1 < Sup_SEC){
				for(i=0; i<i_imag;i++){		
					for(j=0; j<j_imag; j++){
						imagen[calcSingleSubscript(i,j,j_imag)] = imag3D[calcSingleSubscript3D(i,j,n,i_imag,j_imag)];
					}
				}
				}else{
					swappixelOGL(imagen);
				}
				n++;
				nsecuencia++;
				swappixelOGL(imagen);
				fuzzy_seg = false;
				stEtiquet->set_text("Seg. GAR");
				if (iterac==0){
					if(l_vacia(puntosselec)){
						cargarPuntos(INI_SEC,puntosselec);
					}
				//InitGARModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
				//	puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
				  InitGAMRModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
					puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
				}
				iterac=1;
			}
		}else{
			int i,j,l;
			n--;					
			seg = false;
		   PresentGRAModel(imag);
		   if(!l_vacia(puntosselec)){
			   l_dest(&puntosselec);
			   puntosselec = l_nuev(sizeof(punto2D));
		   }
		   DestroyGARModel();
		   iterac_actual = iterac;
		   iterac=0;					   
		   fclose(textresult);
		   fclose(timeresult);
		   free(imag3D);
		   imag3D = NULL;
		   printf("Fin de la segmentacin.\n");
		}			
	break;	  
  }

}
/*---------------- MENU LEVEL SET -------------------------------*/
void controlmenuLS(int value){//Casos del menu de opciones.
  opcion_ls = value;
  switch(opcion_ls){
      case LS_EXPAND:
        if (iterac==0){
          if(l_vacia(puntosselec)){
			cargarPuntos("ptosResLevelSet.txt",puntosselec);
          }
		  if(subopcion == GEODESGVF){
			  inicializar_level_setGVF(imagen,INCREMENTO,0,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
				puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,fuzzy_seg,num_clusters,clust_int,bInterp,APORTE);
		  }else{
			  inicializar_level_set(imagen,INCREMENTO,0,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
				puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,fuzzy_seg,num_clusters,clust_int,bInterp);
		  }

		   //inicializar_ls_ffm(imagen,0,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
          //  puntosselec,j_imag,i_imag,THRESS,THRESS_K);
        }
        iterac=1;
        seg = true;
      break;
      case LS_COMP:
       if (iterac==0){
          if(l_vacia(puntosselec)){
            cargarPuntos("ptosResLevelSet.txt",puntosselec);
          }
		  if(subopcion == GEODESGVF){
			  inicializar_level_setGVF(imagen,INCREMENTO,1,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
				puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,fuzzy_seg,num_clusters,clust_int,bInterp,APORTE);
		  }else{
			inicializar_level_set(imagen,INCREMENTO,1,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
				puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,fuzzy_seg,num_clusters,clust_int,bInterp);
		  }
        }
        iterac=1;
        seg = true;
      break;
	  case GARMODEL:
		 subopcion = GARMODEL;		 
         fuzzy_seg = false;
		 stEtiquet->set_text("Seg. GAR");
         if (iterac==0){
          if(l_vacia(puntosselec)){
            cargarPuntos("ptosResLevelSet.txt",puntosselec);
          }
          
		  InitGARModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
            puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);
        }
        iterac=1;
        seg = true;
      break;
	  
	  case GAMRMODEL:
		 subopcion = GAMRMODEL;		 
         fuzzy_seg = false;
		 stEtiquet->set_text("Seg. GAR");
         if (iterac==0){
          if(l_vacia(puntosselec)){
            cargarPuntos("ptosResLevelSet.txt",puntosselec);
          }
          
		  InitGAMRModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
            puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);
        }
        iterac=1;
        seg = true;
      break;
  }
  glutPostRedisplay();
 }
/*------------------- CONTROL MENU -----------------------*/
//Para controlar el menu
void controlmenu(int value){//Casos del menu de opciones.
   int i;   
	switch(value){
    case CARGAR:
      sprintf(arch,PATH_IMAG);
      strcat(arch,name);
	  strcpy(name_actual,name);
	  if(imagen != NULL){
        delete(imagen);
		imagen=NULL;
	  }
	  if(imag != NULL){
        delete(imagen);
		imag=NULL;
	  }
      switch(extension){
        case 0:

          strcat(arch,".pgm");
		  /**/printf("%s\n",arch);
		  err = cargarImagenPGM(&imagen,arch);
          if(err!= 0){
			  error(err,mensaje);
			  stEtiquet->set_text(mensaje);
		  }
        break;
        case 1:
          strcat(arch,".raw");
          err = cargarImagenRAW(&imagen,arch,Dim_RAW_X,Dim_RAW_Y);
          if(err!= 0){
			  error(err,mensaje);
			  stEtiquet->set_text(mensaje);
		  }
        break;
      }
	  if(err==0){
		stEtiquet->set_text("Fichero Abierto.");
		Get_X_Y(&j_imag,&i_imag);
		swappixelOGL(imagen);
		imag = (unsigned char *)calloc(sizeof(unsigned char),i_imag * j_imag);
		memcpy(imag,imagen,sizeof(unsigned char)*i_imag*j_imag);
		for(i=0; i<(int)strlen(name);i++)
			name[i] = '\0';
	  }
	break;
	case PROCESAR_SEC:
		char secname[40];
		if(imag3D != NULL) free(imag3D);
		if(imagen != NULL) free(imagen); 
		if(imag != NULL) free(imag);
		imag3D = (unsigned char*) calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y*Dim_SEC_Z);
		imagen = (unsigned char*)calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y);
		imag = (unsigned char*)calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y);
		textresult = fopen("PuntosSecuencia.txt","w");
		timeresult = fopen("TiemposIteraciones.txt","w");
		int l;
		for(l=0;l<40;l++) secname[l] = '\0';
		
		  sprintf(secname,PATH_IMAG);
		  strcat(secname,name_secuen);
		  printf("%s\n",secname);
		  FILE* archivo;
		  archivo = fopen(secname,"rb");
		  fread((char*)imag3D,sizeof(char),Dim_SEC_X*Dim_SEC_Y*Dim_SEC_Z,archivo);
		  fclose(archivo);
		  seg = true;
		  subopcion = PROC_SEC;
		  char fichname[40];
			j_imag = Dim_SEC_X;
			i_imag = Dim_SEC_Y;
			Set_X_Y(j_imag,i_imag);			
		    n = Inf_SEC;
			nsecuencia = 0;
		  int i,j;
	for(i=0; i<i_imag;i++)		
		for(j=0; j<j_imag; j++){
			imagen[calcSingleSubscript(i,j,j_imag)] = imag3D[calcSingleSubscript3D(i,j,n,i_imag,j_imag)];
		}
		subopcion = PROC_SEC;		 
	fuzzy_seg = false;
	//rotacion_vert(imagen);
	swappixelOGL(imagen);
	stEtiquet->set_text("Seg. GAR");
	if (iterac==0){
		if(l_vacia(puntosselec)){
			cargarPuntos(INI_SEC,puntosselec);
		}
	//	InitGARModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
    //       puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
		InitGAMRModel(imagen,direccion,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
           puntosselec,j_imag,i_imag,THRESS,CIRC,THRESS_K,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
    }
	iterac=1;

		
	break;
	   
	}
	if(imagen!=NULL){
 	switch(value){
       case SUAVIZAR:
        double *imagen_suav;
        imagen_suav = (double*)calloc(sizeof(double),i_imag*j_imag);
        gaussiana(SIGMA_GAUSS,imag,imagen_suav);
        doubleTounsignedchar(imagen_suav,imag);
        free(imagen_suav);
       break;
       case CURVAT:
         double *imag_curv;
 		    imag_curv = (double*)calloc(sizeof(double),i_imag*j_imag);
         double *imag_curv2;
 		    imag_curv2 = (double*)calloc(sizeof(double),i_imag*j_imag);
         for(i=0; i<(i_imag*j_imag); i++)
           imag_curv[i] = imag[i];
         curv_div(imag_curv,imag_curv2);
         doubleTounsignedchar(imag_curv2,imag);
         free(imag_curv);
         free(imag_curv2);
       break;
       case ENSUCIAR:
         ruido_aleatorio(imagen,imagen);
         memcpy(imag,imagen,(sizeof(unsigned char) * i_imag * j_imag));
       break;
       case GAUSS_GRAD:
        double *image_gradgauss;
 		    float max,min;
 		    image_gradgauss = (double*)calloc(sizeof(double),i_imag*j_imag);
 		    grad_conv_gauss(imag,0.5,&max,&min,image_gradgauss );
        for(i=0; i<(i_imag*j_imag); i++)
           imag[i] = image_gradgauss[i];        
 		    free(image_gradgauss);
 	    break;
       case GRADIEN:
         double *image;
 		    double *image_prw;
         image = (double*)calloc(sizeof(double),i_imag*j_imag);
 		    image_prw = (double*)calloc(sizeof(double),i_imag*j_imag);
         for(i=0; i<(i_imag*j_imag); i++)
           image[i] = imag[i];
         grad_diffin_norm(image,image_prw,image_prw,image_prw);
         doubleTounsignedchar(image_prw,imag);
         free(image);
       break;

       case NUEVA_SELEC:
		   if(puntosselec){
			   l_dest(&puntosselec);
			   puntosselec = l_nuev(sizeof(punto2D));
			   memcpy(imag,imagen,(sizeof(unsigned char) * i_imag * j_imag));
		   }
       break;
       case GUARDAR:
         int a,b;
         Get_X_Y(&a,&b);
         Set_X_Y(j_imag,i_imag);         
         swappixelOGL(imag);
         guardarImagenPGM(imag,name_wind);
         swappixelOGL(imag);
         Set_X_Y(a,b);                  
       break;
	   case GUARDAR_AUTO:
         int m,n;
		 char nombre_auto[80];
         Get_X_Y(&m,&n);
         Set_X_Y(j_imag,i_imag);         
         swappixelOGL(imag);
		 switch(subopcion){
			 case LEVELSET: 
				 switch(opcion_ls){
					 case LS_EXPAND:
						 sprintf(nombre_auto,"ResLSExp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
					 case LS_COMP:
						 sprintf(nombre_auto,"ResLSComp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
				 };
			 break;
			 case GEODES: 
				 switch(opcion_ls){
					 case LS_EXPAND:
						 sprintf(nombre_auto,"ResGEODExp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
					 case LS_COMP:
						 sprintf(nombre_auto,"ResGEODComp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
				 };
			 break;
			 case GEODESGVF: 
				 switch(opcion_ls){
					 case LS_EXPAND:
						 sprintf(nombre_auto,"ResGVFExp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
					 case LS_COMP:
						 sprintf(nombre_auto,"ResGVFComp-%s-%.2f-%i-%.2f-%.2f-%i-%.2f.pgm",name_actual,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,(double)tiempo_total/CLOCKS_PER_SEC);
					 break;
				 };
			 break;
			 case GARMODEL:
				 sprintf(nombre_auto,"ResGAR-%s-%.2f-%.2f-%i-%.2f-%.4f-%i-%.3f.pgm",name_actual,APORTE,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,tiempo);
			 break;
			 case GAMRMODEL:
				 sprintf(nombre_auto,"ResGAMR-%s-%.2f-%.2f-%i-%.2f-%.4f-%i-%.3f.pgm",name_actual,APORTE,SIGMA_GAUSS,NARROWBAND,
							 INCREMENTO,THRESS,iterac_actual,tiempo);
			 break;
		 };
		 guardarImagenPGM(imag,nombre_auto);
         swappixelOGL(imag);
         Set_X_Y(m,n);                  
       break;
       case LEVELSET:
       subopcion = LEVELSET;
       fuzzy_seg = false;
	   stEtiquet->set_text("Seg. Geometricas");
	   switch(direccion){
           case 0:
             controlmenuLS(LS_EXPAND);
           break;
           case 1:
             controlmenuLS(LS_COMP);
           break;
         }
       break;
       case GEODES:
         subopcion = GEODES;
         fuzzy_seg = false;
		 stEtiquet->set_text("Seg. Geodesicas");
         switch(direccion){
           case 0:
             controlmenuLS(LS_EXPAND);
           break;
           case 1:
             controlmenuLS(LS_COMP);
           break;
         }
       break;
	   case GEODESGVF:
         subopcion = GEODESGVF;
         fuzzy_seg = false;
		 stEtiquet->set_text("Seg. Geodesicas");
         switch(direccion){
           case 0:
             controlmenuLS(LS_EXPAND);
           break;
           case 1:
             controlmenuLS(LS_COMP);
           break;
         }
       break;
	   


       case FUZZYLEVELSET:
       subopcion = LEVELSET;
       fuzzy_seg = true;
	   stEtiquet->set_text("Seg. FuzzyGeometricas");
        switch(direccion){
           case 0:
             controlmenuLS(LS_EXPAND);
           break;
           case 1:
             controlmenuLS(LS_COMP);
           break;
         }
       break;
       case FUZZYGEODES:
         fuzzy_seg = true;
         subopcion = GEODES;
		 stEtiquet->set_text("Seg. FuzzyGeodesicas");
         switch(direccion){
           case 0:
             controlmenuLS(LS_EXPAND);
           break;
           case 1:
             controlmenuLS(LS_COMP);
           break;
         }
       break;
   
       case SNAKES:
		  fuzzy_seg = false;         
		   if(iterac==0){
             if (l_vacia(puntosselec)){
               cargarPuntos("ptosSnake.txt",puntosselec);
             }
             inicializarSnake(imagen,puntosselec,i_imag,j_imag,
             alpha_sn,beta_sn,teta_sn,w_term,w_line,w_edge,w_press,fuzzy_seg,num_clusters,clust_int);			 
 		    }
         iterac=1;
         seg = true;
     		subopcion = SNAKES;
			stEtiquet->set_text("Seg. Snakes");
       break;
       case FUZZYSNAKES:
		  fuzzy_seg = true;         
		   if (iterac==0){
             if (l_vacia(puntosselec)){
               cargarPuntos("ptosSnake.txt",puntosselec);
             }
             inicializarSnake(imagen,puntosselec,i_imag,j_imag,
             alpha_sn,beta_sn,teta_sn,w_term,w_line,w_edge,w_press,fuzzy_seg,num_clusters,clust_int );			
 		    }
         iterac=1;
         seg = true;
     	 subopcion = SNAKES;
		 stEtiquet->set_text("Seg. Snakes");			
       break;
       case PARARSEG:
		   if (seg==true){
			   seg = false;
			   stEtiquet->set_text("Segmentacion Detenida.");
			   switch(subopcion){
			   case LEVELSET: case GEODES: case GEODESGVF:
					   printf("Fin de la segmentacin.\n");
					   printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
					   sprintf(mensaje,"Iteraciones: %i",iterac);
					   stIterac->set_text(mensaje);
					   sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
					   stTiempo->set_text(mensaje);      
					   filtrado_burb();
					   actual_level_set(imag);       
					   //actual_level_set_FFM(imag);       
					   almacenarLevelSet("ptosResLevelSet.txt");
					   liberar_ls();					   
					   iterac_actual = iterac;
						iterac=0;
					break;
					case SNAKES:
						libera_snk();
						iterac_actual = iterac;
						iterac=0;
					break;
					case GARMODEL: case GAMRMODEL:
					   printf("Fin de la segmentacin.\n");
					   printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
					   sprintf(mensaje,"Iteraciones: %i",iterac);
					   stIterac->set_text(mensaje);
					   sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
					   stTiempo->set_text(mensaje);      
					   FilterBubbles();
					   PresentGRAModel(imag);
					   if(!l_vacia(puntosselec)){
						   l_dest(&puntosselec);
						   puntosselec = l_nuev(sizeof(punto2D));
					   } 
					   SetUpCurveGARModel(puntosselec);
					   Set_X_Y(j_imag,i_imag);
						//printf("Numero de LS: %d\n",NumLS(imag));
	  
	  				   DestroyGARModel();
					   iterac_actual = iterac;
						iterac=0;
					break;
					case PROC_SEC:
					   printf("Fin de la segmentacin.\n");
					   printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
					   sprintf(mensaje,"Iteraciones: %i",iterac);
					   stIterac->set_text(mensaje);
					   sprintf(mensaje,"Tiempo: %.3f seg.",(double)tiempo_total/CLOCKS_PER_SEC);
					   stTiempo->set_text(mensaje);      
					   PresentGRAModel(imag);
					   if(!l_vacia(puntosselec)){
						   l_dest(&puntosselec);
						   puntosselec = l_nuev(sizeof(punto2D));
					   }
					   //SetUpCurveGARModel(puntosselec);
	  				   DestroyGARModel();
					   iterac_actual = iterac;
						iterac=0;					   
					   fclose(textresult);
					   free(imag3D);
					   imag3D = NULL;
					break;
			   }
		   }
       break;
       case GUARDAR_SEL_GEOD:
         guardarPuntos(puntosselec,"ptosGeodesic.txt");
       break;
       case GUARDAR_SEL_LS:
         guardarPuntos(puntosselec,"ptosLevelSet.txt");
       break;
 	   case GUARDAR_SEL_SN:
         guardarPuntos(puntosselec,"ptosSnake.txt");
       break;
       case CIRCULOS:
         multiples_circulos(NUM_CIRC,puntosselec);
         CIRC = true;
       break;
       case FUZZY:
 		      fuzzyCmeans(imag,i_imag,j_imag,imag,num_clusters,10);
       break;
       case FUZZY_FORZ:
          int num_cl;
          unsigned char *centros;
          cargar_centros_fuzzy(&centros,"centros.txt",&num_cl);
          fuzzyCmeansforzado(imag,i_imag,j_imag,imag,centros,num_cl,10);          
       break;
	   case BINAR:
		   int histo[256];
		   histograma(imag,histo);
		   binarizaImag(imag, umbralOtsu(histo));
       break;
	    case LAPLAC:
			double *image_lap;
			image_lap = (double*)calloc(sizeof(double),i_imag*j_imag);
			laplaciano(imag,image_lap);
			doubleTounsignedchar(image_lap,imag);
			free(image_lap);
       break;
	   case NEGAT:
		   negativo(imag,imag);
       break;
	   case ENTROPIA:
		   float entro;
		   entro = entropia(imag);
		   sprintf(mensaje," %.4f",entro);
		   stEntropia->set_text(mensaje);      
	   break;
	   case DIB_CIRCULO:
		   //crear_circunferencia_cerrada(int orig_x, int orig_y, float radio,lista curva);
	   break;
	   case GUARDA_SELEC:
		   char nom_sel[50];
		   sprintf(nom_sel, PATH_SELEC);
		   strcat(nom_sel,name_selec);
		   guardarPuntos(puntosselec,nom_sel);
		   //guardarPuntosImag(nom_sel);
	   break;
	   case CARGA_SELEC:
		   if(name_selec[0] != '\0' && (imagen!= NULL)){
			   char nom_sel[40];
			   sprintf(nom_sel,PATH_SELEC);
			   strcat(nom_sel,name_selec);
			   cargarPuntos(nom_sel,puntosselec);
		   }
	   break;
	   case GVF:
		   puntof2D* GVFMtx;
		   GVFMtx = (puntof2D*)malloc(sizeof(puntof2D)*i_imag*j_imag);
		   GradientVF(imag,GVFMtx,0.2);
		   DrawVectImag(GVFMtx,imag);
		   StoreVectImg(GVFMtx,"GVF.txt");
		   free(GVFMtx);
	   break;
	   case GVFPARAG:
		   puntof2D* GVFMtxParag;
		   GVFMtxParag = (puntof2D*)malloc(sizeof(puntof2D)*i_imag*j_imag);
		   GradientVF_Paragios(imag,GVFMtxParag,0.2);
		   DrawVectImag(GVFMtxParag,imag);
		   StoreVectImg(GVFMtxParag,"GVF.txt");
		   free(GVFMtxParag);
	   break;	   
	   case GVFMOD:
		   puntof2D* GVFMtxMod;
		   GVFMtxMod = (puntof2D*)malloc(sizeof(puntof2D)*i_imag*j_imag);
		   double* Module;
		   Module = (double*)malloc(sizeof(double)*i_imag*j_imag);
		   GradientVF(imag,GVFMtxMod,0.2);
		   VectorModule(GVFMtxMod,Module);
		   /**/DoubleMtxToTXT("ModuloGVF.txt",Module);
		   doubleTounsignedchar(Module,imag);
		   free(GVFMtxMod);
		   free(Module);
	   break;
	   case GVFMODPARAG:
		   puntof2D* GVFMtxModParag;
		   GVFMtxModParag = (puntof2D*)malloc(sizeof(puntof2D)*i_imag*j_imag);
		   double* ModuleParag;
		   ModuleParag = (double*)malloc(sizeof(double)*i_imag*j_imag);
		   GradientVF_Paragios(imag,GVFMtxModParag,0.2);
		   VectorModule(GVFMtxModParag,ModuleParag);
		   /**/DoubleMtxToTXT("ModuloGVF.txt",ModuleParag);
		   doubleTounsignedchar(ModuleParag,imag);
		   free(GVFMtxModParag);
		   free(ModuleParag);
	   break;
	   
	   }    
  }
  glutPostRedisplay();
}
/*---------------- VENTANA PRINCIPAL -------------------------------*/
void mainDisplay (void)
{
  /* Clean drawing board */
  glutSetWindow (winIdMain);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  if(seg){
       l_dest(&puntosselec);
       puntosselec = l_nuev(sizeof(punto2D));	  
	   bucle();
  }else{
    if(!l_vacia(puntosselec)){
      l_empieza(puntosselec);
      while(l_quedan(puntosselec)){
        punto2D pto;
        l_dame(puntosselec,&pto);
        imag[calcSingleSubscript(pto.x,pto.y,j_imag)] = 255.0;
      }
    }

  }
  desplaz_x = -1.0 + (((weight - j_imag)/2.0)/(weight/2.0));
  desplaz_y = -1.0 + (((height - i_imag)/2.0)/(height/2.0));;
  glRasterPos2f(desplaz_x,desplaz_y);
  glDrawPixels(j_imag,i_imag,GL_LUMINANCE,GL_UNSIGNED_BYTE,imag);
  glutSwapBuffers();
  
}
/*---------------- INIT -------------------------------*/
//Inicializacion.
void init() {
	//Definicion de la luz
	GLfloat light_ambient[] = { 0.75, 0.75, 0.75, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.1, 0.25, 1.0, 0.0 };

	glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv (GL_LIGHT0, GL_POSITION, light_position);

	//Activamos las luces
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT0);

  //Tipo de modelos
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	//Color de fondo
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
//  glPixelStorei(GL_UNPACK_ALIGNMENT,2);
  puntosselec = l_nuev(sizeof(punto2D));
  tiempo_total = 0;
}

/*------------- RESHAPE -------------------------------------*/
void mainReshape (int w, int h)
{
  int tx,ty,tw,th;
  GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);

  glViewport (tx, ty, tw, th);

  weight = w;
  height = h;
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D (-1.0F, 1.0F, -1.0F, 1.0F);
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();

}

/*------------- idle -------------------------------------*/
void idle (void){

  /* Update main and sub window */
  glutSetWindow (winIdMain);
  glutPostRedisplay ();

}

/*------------- MOTION -------------------------------------*/


//Rutina de control del movimiento del raton
void motion(int x, int y)
 {
	if(!dib_circulo){
		punto2D pto;

		pto.x = (j_imag - y) + (weight - j_imag) - ((height  - i_imag)/2.0);
		pto.y = x  - ((weight - j_imag)/2.0);
		l_meted(puntosselec,&pto);   
		/**/printf("punto %d\t %d\n",pto.x, pto.y);
	}
 }


/*--------------- MOUSE -----------------------------------*/


//Rutina de control del raton
void mouse(int button, int state, int x, int y)
{
  if(dib_circulo){	  
	  int x_center = (j_imag - y) + (weight - j_imag) - ((height  - i_imag)/2.0);
	  int y_center = x  - ((weight - j_imag)/2.0);
	  crear_circunferencia_centrada(x_center,y_center,iRadio,iMuestreoCirc,puntosselec);	  
  }else{   
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		motion(x,y);
	}
  }

}
/*--------------- MAIN ---------------------------------*/
void MenuMatricesGAR(int value){
	int i,j,position ;
	punto2D pto;

	switch(value){
	case GAR_BT:

		double *Fi_new;
		double *Fi_old;
		double *Fi_grad;
		double *Fi_temp;
		int *Fi_int;
		double *Veloc;
		double *Ki;
		double *Ki_grad;
		int *etiquetas;
		punto2D *dist_4sed;
  		int valor;
		float max;
		lista LPuntosInterp;

		Fi_new = (double*)calloc(sizeof(double),j_imag*i_imag);
		Fi_old= (double*)calloc(sizeof(double),j_imag*i_imag);
		Fi_grad= (double*)calloc(sizeof(double),j_imag*i_imag);
		Fi_temp= (double*)calloc(sizeof(double),j_imag*i_imag);
		Veloc= (double*)calloc(sizeof(double),j_imag*i_imag);
		Ki= (double*)calloc(sizeof(double),j_imag*i_imag);
		Ki_grad= (double*)calloc(sizeof(double),j_imag*i_imag);
		etiquetas = (int*)calloc(sizeof(int),j_imag*i_imag);
		dist_4sed = (punto2D*)calloc(sizeof(punto2D),j_imag*i_imag);
		Fi_int = (int*)calloc(sizeof(int),j_imag*i_imag);
		
		grad_conv_gauss(imag,SIGMA_GAUSS,&max,&max,Fi_grad);
		LPuntosInterp = l_nuev(sizeof(punto2D));
		cierre_lineal(puntosselec,l_elementos(puntosselec),LPuntosInterp);
		l_empieza(LPuntosInterp);
		while(l_quedan(LPuntosInterp)){
			l_dame(LPuntosInterp,&pto);
			Fi_temp[calcSingleSubscript(pto.x,pto.y,j_imag)] = Fi_old[calcSingleSubscript(pto.x,pto.y,j_imag)] = 1.0;
			Ki[calcSingleSubscript(pto.x,pto.y,j_imag)] = 1.0/1.0 + Fi_grad[calcSingleSubscript(pto.x,pto.y,j_imag)];
			Fi_int[calcSingleSubscript(pto.x,pto.y,j_imag)] = 1;
		}		
		chamfer_distance4v(Fi_old,Fi_old);
		etiq_region(Fi_int,etiquetas);
		valor = etiquetas[calcSingleSubscript(0,0,j_imag)];
		for(i=0; i<j_imag*i_imag;i++){
			if(etiquetas[i] == valor){
				Fi_old[i] *= -1;
			}
		}
		
		curvatura(Fi_old,Veloc,grad_diffin);
		for(i=0;i<j_imag*i_imag;i++){
			Veloc[i] *= 0.01;
		}
		
		vecin4ssed_distance(Fi_temp,dist_4sed);
		int pos;
		for(i=1;i<i_imag-1;i++){
			for(j=1;j<j_imag-1;j++){
				pos= calcSingleSubscript(i,j,j_imag);
				Ki[pos] = Ki[calcSingleSubscript(i - dist_4sed[pos].x,j - dist_4sed[pos].y,j_imag)];
			}
		}
		grad_diffin(Fi_old,Fi_temp,Fi_temp,Fi_grad);
		grad_diffin(Ki,Fi_temp,Fi_temp,Ki_grad);

		for(i=0; i<j_imag*i_imag; i++){
			Fi_new[i] = (Ki[i]*(1.0 + Veloc[i])*Fi_grad[i])+((Ki_grad[i]*Fi_grad[i])+(Ki_grad[i]*Fi_grad[i]));		
			imag[i] = Fi_new[i];
		}
		free(Fi_new);
		free(Fi_old);
		free(Fi_grad);
		free(Fi_temp);
		free(Veloc);
		free(Ki);
		free(Ki_grad);
		free(etiquetas);
		free(dist_4sed);
		l_dest(&LPuntosInterp);

		break;
	case GAR_RT:

		float Media, Prob;
		lista LPtosInterp;

		int *Binar;
		int *Map;
		int *etiq;
		
		LPtosInterp = l_nuev(sizeof(punto2D));
		Binar = (int*)calloc(sizeof(int),j_imag*i_imag);
		etiq = (int*)calloc(sizeof(int),j_imag*i_imag);
		Map = (int*)malloc(sizeof(int)*j_imag*i_imag);


		cierre_lineal(puntosselec,l_elementos(puntosselec),LPtosInterp);
		l_empieza(LPtosInterp);
		while(l_quedan(LPtosInterp)){
			l_dame(LPtosInterp,&pto);
			Binar[calcSingleSubscript(pto.x,pto.y,j_imag)] = 1;
		}
		for(i=0; i<i_imag*j_imag; i++)
			Map[i]=1;
		
		etiq_region(Binar,etiq);
		valor = etiq[calcSingleSubscript(0,0,j_imag)];
		for(i=0; i<j_imag*i_imag;i++){
			if(etiq[i] == valor){
				Map[i] *= -1;
			}
		}

		double dTotal=0;
		int iNum=0;
		for(i=0;i<i_imag*j_imag;i++){
			if(Map[i]<0){
				iNum++;
				dTotal += imag[i];
			}
		}
		Media = (dTotal/iNum);

		int vHistRegion[255];
		for(i=0; i<255; i++) vHistRegion[i]= 0;
		for(i=0; i<j_imag*i_imag;i++){
			if(Map[i]<0){
				vHistRegion[imag[i]]++;	
			}
		}	
		int iMax = 0;
		int iIntensity;
		for(i=0;i<255; i++){
			if(vHistRegion[i] > iMax){
				iMax = vHistRegion[i];
				iIntensity = i;
			}
		}

		Prob = iIntensity;	
		double dDifMedia, dDifThres;

		for(i=0; i<j_imag*i_imag;i++){
			dDifMedia = imag[i] - Media;
			dDifThres = imag[i] - Prob;
			imag[i] = (sqrt((dDifMedia * dDifMedia)+(dDifThres * dDifThres)));	
		}

		l_dest(&LPtosInterp);
		free(Binar);
		free(etiq);
		free(Map);
		break;
	
	  }

}

/*--------------- MAIN ---------------------------------*/
void MenuMatricesLS(int value){

	int i,j,position ;
	punto2D pto;

	switch(value){
	case MTX_CURV:
		double * curvMtx;
		
		curvMtx = (double*)malloc(sizeof(double)*j_imag*i_imag);
		for(i=0; i<j_imag*i_imag; i++)
			curvMtx[i] = imag[i];
		curvatura(curvMtx,curvMtx,grad_diffin_norm);
		doubleTounsignedchar(curvMtx, imag);
		free(curvMtx);
		break;
	case MTX_FI:
		//Inicializamos la matriz
		double *fiMtx;
		int *fiTemp;
		lista LPuntosInterp;
		int *etiquetas;
		int valor;
		
		fiMtx = (double*)calloc(sizeof(double),j_imag*i_imag);
		fiTemp = (int*)calloc(sizeof(int),j_imag*i_imag);
		etiquetas = (int*)calloc(sizeof(int),j_imag*i_imag);
		LPuntosInterp = l_nuev(sizeof(punto2D));
		//Pintamos los puntos de la curva inicial.
		cierre_lineal(puntosselec,l_elementos(puntosselec),LPuntosInterp);
		l_empieza(LPuntosInterp);
		while(l_quedan(LPuntosInterp)){
			l_dame(LPuntosInterp,&pto);
			position = calcSingleSubscript(pto.x,pto.y,j_imag);
			fiTemp[position] = fiMtx[position] = 1;
		}
		chamfer_distance4v(fiMtx,fiMtx);
		etiq_region(fiTemp,etiquetas);
		valor = etiquetas[calcSingleSubscript(0,0,j_imag)];
		for(i=0; i<j_imag*i_imag;i++){
			if(etiquetas[i] == valor){
				fiMtx[i] *= -1;
			}
		}
		doubleTounsignedchar(fiMtx, imag);
		free(fiMtx);
		free(fiTemp);
		free(etiquetas);
		l_dest(&LPuntosInterp);
		break;	
	case MTX_KI:
		double *kiext;
		float max, min;
		kiext = (double*)malloc(sizeof(double)*i_imag*j_imag);
		grad_conv_gauss(imag,SIGMA_GAUSS,&max,&min,kiext);
		for(i=0; i<i_imag*j_imag; i++){
			kiext[i]= 1.0/(1.0 + kiext[i]);
		}
		doubleTounsignedchar(kiext,imag);
		free(kiext);
		break;
	case MTX_KIEXT:
		double *ki;
		double *grad;
		double *bin;
		punto2D *dist_4sed;
		lista LPtoInterp;

		ki = (double*)calloc(sizeof(double),i_imag*j_imag);
		grad = (double*)calloc(sizeof(double),i_imag*j_imag);
		bin = (double*)calloc(sizeof(double),i_imag*j_imag);
		dist_4sed = (punto2D*)calloc(sizeof(punto2D),i_imag*j_imag);
		grad_conv_gauss(imag,SIGMA_GAUSS,&max,&min,grad);
		LPtoInterp = l_nuev(sizeof(punto2D));
		//Pintamos los puntos de la curva inicial.
		cierre_lineal(puntosselec,l_elementos(puntosselec),LPtoInterp);
		
		l_empieza(LPtoInterp);
		while(l_quedan(LPtoInterp)){
			l_dame(LPtoInterp,&pto);
			position = calcSingleSubscript(pto.x,pto.y,j_imag);
			ki[position]= 1.0/(1.0 + grad[position]);
			bin[position]= 1.0;
		}

		vecin4ssed_distance(bin,dist_4sed);
		for(i=1;i<i_imag-1;i++){
			for(j=1;j<j_imag-1;j++){
				position= calcSingleSubscript(i,j,j_imag);
				ki[position] = ki[calcSingleSubscript(i - dist_4sed[position].x,j - dist_4sed[position].y,j_imag)];
			}
		}
		doubleTounsignedchar(ki,imag);
		free(grad);
		free(bin);
		free(dist_4sed);	
		free(ki);
		l_dest(&LPtoInterp);
    	break;
	}
}


/*--------------- MAIN ---------------------------------*/
int main(int argc,char** argv){
	for(int i=0; i<(int)strlen(name); i++)
    name[i] = '\0';
	glutInit(&argc,argv);
	//inicializa la ventana en el modo Doble Buffer
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	//Tamao de la ventana
	//glutInitWindowSize (1100, 1100);
	glutInitWindowSize (800, 800);
	//Posicion de la ventana
	glutInitWindowPosition (0, 0);
	//Crea la Ventana
	winIdMain = glutCreateWindow ("Ventana de Visualizacion");
	//sentencias de Inicio
	init ();
	//Funcion de dibujo 
	glutDisplayFunc(mainDisplay);
	//Funcion cuando cambiemos el tamao de la ventana
  //GLUI_Master.set_glutReshapeFunc(mainReshape);
  glutReshapeFunc(mainReshape);
	//Control del teclado
 // GLUI_Master.set_glutKeyboardFunc(keyboard);
  //Control del Raton
  GLUI_Master.set_glutMouseFunc(mouse);
	glutMotionFunc(motion);
  //Control de las ventanas
  GLUI_Master.set_glutIdleFunc(idle);
  //Creacion del menu de level set
  

  //Creacion de la subventana
  glui = GLUI_Master.create_glui("Ventana de Control",0,610,0);
  glui->set_main_gfx_window(winIdMain);
 //Panel de Apertura
  GLUI_Rollout *pnAbrir = glui->add_rollout("Abrir",true);
  glui->add_statictext_to_panel(pnAbrir,"Archivo:");
  glui->add_edittext_to_panel(pnAbrir,"",GLUI_EDITTEXT_TEXT,
    &name);
  GLUI_RadioGroup *rgExtension = glui->add_radiogroup_to_panel(pnAbrir,&extension);
  glui->add_radiobutton_to_group(rgExtension,"PGM");
  glui->add_radiobutton_to_group(rgExtension,"RAW");  
  glui->add_spinner_to_panel(pnAbrir,"Dimension X = ",
    GLUI_SPINNER_INT,&Dim_RAW_X);
  glui->add_spinner_to_panel(pnAbrir,"Dimension Y = ",
    GLUI_SPINNER_INT,&Dim_RAW_Y);
  glui->add_button_to_panel(pnAbrir,"Abrir",CARGAR,(GLUI_Update_CB)controlmenu);
  glui->add_separator_to_panel(pnAbrir);
  glui->add_statictext_to_panel(pnAbrir,"Nombre:");
  glui->add_edittext_to_panel(pnAbrir,"",GLUI_EDITTEXT_TEXT,
    &name_wind);
  glui->add_button_to_panel(pnAbrir,"Guardar Pantalla",GUARDAR,controlmenu);
  glui->add_button_to_panel(pnAbrir,"Guardar Automatico",GUARDAR_AUTO,controlmenu);
  
 
  //Funciones Generales
  GLUI_Rollout *roFuncSelecc = glui->add_rollout("Funciones de Seleccion",true);
  glui->add_button_to_panel(roFuncSelecc,"Nueva Seleccion",
     NUEVA_SELEC,(GLUI_Update_CB)controlmenu);
  glui->add_checkbox_to_panel(roFuncSelecc,"Dib. Circulo",&dib_circulo);
  glui->add_spinner_to_panel(roFuncSelecc,"Radio: ",
    GLUI_SPINNER_INT,&iRadio);
  glui->add_spinner_to_panel(roFuncSelecc,"Muestreo: ",
    GLUI_SPINNER_INT,&iMuestreoCirc);
  glui->add_edittext_to_panel(roFuncSelecc,"",GLUI_EDITTEXT_TEXT,
    &name_selec);
  glui->add_button_to_panel(roFuncSelecc,"Guardar Seleccion",
     GUARDA_SELEC,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncSelecc,"Cargar Seleccion",
     CARGA_SELEC,(GLUI_Update_CB)controlmenu);
  
  glui->add_button_to_panel(roFuncSelecc,"Seleccion Circulos",
    CIRCULOS,(GLUI_Update_CB)controlmenu);
  glui->add_spinner_to_panel(roFuncSelecc,"Circulos: ",
    GLUI_SPINNER_INT,&NUM_CIRC);
  GLUI_Rollout *roFuncGenerales = glui->add_rollout("Funciones Generales",false);
  glui->add_button_to_panel(roFuncGenerales,"Binarizar por OTSU",
     BINAR,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Gradiente",
    GRADIEN,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Laplaciano",
    LAPLAC,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Negativo",
    NEGAT,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Ensuciar",
    ENSUCIAR,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"GVF",
    GVF,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Modulo GVF",
    GVFMOD,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"GVF-Parag",
    GVFPARAG,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Mod GVF-Parag",
    GVFMODPARAG,(GLUI_Update_CB)controlmenu);
  
   glui->add_separator_to_panel(roFuncGenerales);
  glui->add_button_to_panel(roFuncGenerales,"Gausiana",
    SUAVIZAR,(GLUI_Update_CB)controlmenu);
  GLUI_Spinner *spGauss = glui->add_spinner_to_panel(roFuncGenerales,"Sigma: ",
    GLUI_SPINNER_FLOAT,&SIGMA_GAUSS);
  spGauss->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spGauss->set_speed(1.0);
  glui->add_button_to_panel(roFuncGenerales,"GradGauss",
    GAUSS_GRAD,(GLUI_Update_CB)controlmenu);
  GLUI_Spinner *spGradGauss = glui->add_spinner_to_panel(roFuncGenerales,"Sigma: ",
    GLUI_SPINNER_FLOAT,&SIGMA_GAUSS);
  spGradGauss->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spGradGauss->set_speed(1.0);
   glui->add_separator_to_panel(roFuncGenerales);
    glui->add_button_to_panel(roFuncGenerales,"Fuzzy Cluster",
    FUZZY,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Fuzzy Forzado",
    FUZZY_FORZ,(GLUI_Update_CB)controlmenu);
  GLUI_Spinner *spFuzzy = glui->add_spinner_to_panel(roFuncGenerales,"Num. Clusters: ",
    GLUI_SPINNER_INT,&num_clusters);
  spFuzzy->set_float_limits(0.0,100.0,GLUI_LIMIT_CLAMP);
  spFuzzy->set_speed(1.0);
 
  glui->add_button("Salir",0,exit);
  glui->add_separator();
  stEtiquet = glui->add_statictext("Mensajes...");
 
  glui->add_column(true);
  //Panel de Segmentacion
  GLUI_Rollout *roSegmen = glui->add_rollout("Segmentacion",true);
  //Panel de Snakes
  GLUI_Rollout *pnSnakes = glui->add_rollout_to_panel(roSegmen,"PDM",false);
  GLUI_Spinner *spAlphaSn = glui->add_spinner_to_panel(pnSnakes,"Alfa: ",
  GLUI_SPINNER_FLOAT,&alpha_sn);
  spAlphaSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spAlphaSn->set_speed(1.0);
  GLUI_Spinner *spBetaSn = glui->add_spinner_to_panel(pnSnakes,"Beta: ",
  GLUI_SPINNER_FLOAT,&beta_sn);
  spBetaSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spBetaSn->set_speed(1.0);
  GLUI_Spinner *spTetaSn = glui->add_spinner_to_panel(pnSnakes,"Teta: ",
  GLUI_SPINNER_FLOAT,&teta_sn);
  spTetaSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spTetaSn->set_speed(1.0);

  GLUI_Spinner *spWlineSn = glui->add_spinner_to_panel(pnSnakes,"Wlineas: ",
  GLUI_SPINNER_FLOAT,&w_line);
  spWlineSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spWlineSn->set_speed(1.0);
  GLUI_Spinner *spWedgeSn = glui->add_spinner_to_panel(pnSnakes,"Wbordes: ",
  GLUI_SPINNER_FLOAT,&w_edge);
  spWedgeSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spWedgeSn->set_speed(1.0);
  GLUI_Spinner *spWtermSn = glui->add_spinner_to_panel(pnSnakes,"Wterm: ",
  GLUI_SPINNER_FLOAT,&w_term);
  spWtermSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spWtermSn->set_speed(1.0);
  GLUI_Spinner *spWpressSn = glui->add_spinner_to_panel(pnSnakes,"Wpresion: ",
  GLUI_SPINNER_FLOAT,&w_press);
  spWpressSn->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spWpressSn->set_speed(1.0);


  glui->add_button_to_panel(pnSnakes,"Segmentar por PDM",SNAKES,controlmenu);
  glui->add_button_to_panel(pnSnakes,"Seg. con Fuzzy PDM",FUZZYSNAKES,controlmenu);
  glui->add_button_to_panel(pnSnakes,"Guardar Seleccion",GUARDAR_SEL_SN,controlmenu);
  GLUI_Spinner *spClustInt = glui->add_spinner_to_panel(roSegmen,"Cluster Interes: ",
  GLUI_SPINNER_INT,&clust_int);
  spClustInt->set_int_limits(0,num_clusters,GLUI_LIMIT_CLAMP);
  spClustInt->set_speed(1.0);
  //Panel de Level Set
  GLUI_Rollout *pnLevelSet = glui->add_rollout_to_panel(roSegmen,"M. Geometrico y Geodesico",true);
  GLUI_Panel *pnOpcionLS = glui->add_panel_to_panel(pnLevelSet,"",GLUI_PANEL_NONE);  
  glui->add_checkbox_to_panel(pnOpcionLS,"Interpolar Entrada",&bInterp);
  glui->add_separator_to_panel(pnOpcionLS);
  GLUI_RadioGroup *rgDireccion = glui->add_radiogroup_to_panel(pnOpcionLS,&direccion);
  glui->add_radiobutton_to_group(rgDireccion,"Expansion");
  glui->add_radiobutton_to_group(rgDireccion,"Contraccion");
  glui->add_separator_to_panel(pnOpcionLS);
  GLUI_Spinner *spGaussLS = glui->add_spinner_to_panel(pnOpcionLS,"Sigma: ",
  GLUI_SPINNER_FLOAT,&SIGMA_GAUSS);
  spGaussLS->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spGaussLS->set_speed(1.0);
  GLUI_Spinner *spNBLS = glui->add_spinner_to_panel(pnOpcionLS,"Narrow Band: ",
  GLUI_SPINNER_INT,&NARROWBAND);
  spNBLS->set_float_limits(0.0,100.0,GLUI_LIMIT_CLAMP);
  spNBLS->set_speed(1.0);
   GLUI_Spinner *spIncremLS = glui->add_spinner_to_panel(pnOpcionLS,"Increm. T: ",
  GLUI_SPINNER_FLOAT,&INCREMENTO);
  spIncremLS->set_float_limits(0.0,10.0,GLUI_LIMIT_CLAMP);
  spIncremLS->set_speed(1.0);
  GLUI_Spinner *spEpsilonLS = glui->add_spinner_to_panel(pnOpcionLS,"Epsilon: ",


  GLUI_SPINNER_FLOAT,&VELOC_EPSILON);
  spEpsilonLS->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spEpsilonLS->set_speed(1.0);

  GLUI_Spinner *spUmbralLS = glui->add_spinner_to_panel(pnOpcionLS,"Umbral Parada: ",
  GLUI_SPINNER_FLOAT,&THRESS);
  spUmbralLS->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spUmbralLS->set_speed(1.0);
  GLUI_Spinner *spUmbralKI = glui->add_spinner_to_panel(pnOpcionLS,"Umbral Ki: ",
  GLUI_SPINNER_FLOAT,&THRESS_K);
  spUmbralKI->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spUmbralKI->set_speed(1.0);
  GLUI_Spinner *spUmbralIter = glui->add_spinner_to_panel(pnOpcionLS,"Num. Iter: ",
  GLUI_SPINNER_INT,&THRESS_ITER);
  
  
  glui->add_separator_to_panel(pnLevelSet);
  glui->add_button_to_panel(pnLevelSet,"Guardar Seleccion",GUARDAR_SEL_LS,controlmenu);
  glui->add_button_to_panel(pnLevelSet,"Segmentar Modelo Geometrico",LEVELSET,controlmenu);
  glui->add_button_to_panel(pnLevelSet,"Segmentar Fuzzy-Geometrico",FUZZYLEVELSET,controlmenu);
  glui->add_button_to_panel(pnLevelSet,"Segmentar Modelo Geodesico",GEODES,controlmenu);
  glui->add_button_to_panel(pnLevelSet,"Segmentar Fuzzy-Geodesico",FUZZYGEODES,controlmenu);
  glui->add_button_to_panel(pnLevelSet,"Geodes + GVF",GEODESGVF,controlmenu);
glui->add_column_to_panel(roSegmen,true);
  GLUI_Rollout *pnGAR = glui->add_rollout_to_panel(roSegmen,"M. GAR y GAMR",true);

  GLUI_Spinner *spPesoGAR = glui->add_spinner_to_panel(pnGAR,"Aporte: ",
  GLUI_SPINNER_FLOAT,&APORTE);
  spPesoGAR->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spPesoGAR->set_speed(1.0);
  glui->add_button_to_panel(pnGAR,"Geodesic Active Region",GARMODEL,controlmenuLS);
  glui->add_button_to_panel(pnGAR,"G. A. Multiple Region",GAMRMODEL,controlmenuLS);
  GLUI_Panel *pnOpcionGAR = glui->add_panel_to_panel(pnGAR,"Funciones GAR",GLUI_PANEL_EMBOSSED);  
  GLUI_RadioGroup *rgTipoFuncRegion = glui->add_radiogroup_to_panel(pnOpcionGAR,&FuncReg);
  glui->add_radiobutton_to_group(rgTipoFuncRegion ,"Media");
  glui->add_radiobutton_to_group(rgTipoFuncRegion ,"Max. Intensidad");  
  GLUI_Spinner *spValorRegGAR = glui->add_spinner_to_panel(pnOpcionGAR,"Valor Region: ",
  GLUI_SPINNER_FLOAT,&VALREG);
  spPesoGAR->set_float_limits(0.0,255.0,GLUI_LIMIT_CLAMP);
  spPesoGAR->set_speed(1.0);
  glui->add_checkbox_to_panel(pnOpcionGAR,"Solo al Inicio",&bOnlyInitGAR);

  glui->add_separator_to_panel(pnGAR);
  glui->add_statictext_to_panel(pnGAR,"Secuencia:");
  glui->add_edittext_to_panel(pnGAR,"",GLUI_EDITTEXT_TEXT,
    &name_secuen);
  glui->add_spinner_to_panel(pnGAR,"Dimension X = ",
    GLUI_SPINNER_INT,&Dim_SEC_X);
  glui->add_spinner_to_panel(pnGAR,"Dimension Y = ",
    GLUI_SPINNER_INT,&Dim_SEC_Y);
  glui->add_spinner_to_panel(pnGAR,"Cortes = ",
    GLUI_SPINNER_INT,&Dim_SEC_Z);
  GLUI_Panel *pnIntervSECGAR = glui->add_panel_to_panel(pnGAR,"Intervalo Cortes",GLUI_PANEL_EMBOSSED);
  glui->add_spinner_to_panel(pnIntervSECGAR ,"",
    GLUI_SPINNER_INT,&Inf_SEC);
  glui->add_spinner_to_panel(pnIntervSECGAR ,"",
    GLUI_SPINNER_INT,&Sup_SEC);
      
  glui->add_button_to_panel(pnGAR,"Procesar Secuencia",PROCESAR_SEC,controlmenu);

  glui->add_button_to_panel(roSegmen,"Parar Segmentacion",PARARSEG,controlmenu);
  stIterac = glui->add_statictext_to_panel(roSegmen,"Tiempo:");
  stTiempo = glui->add_statictext_to_panel(roSegmen,"Iterac:");
  glui->add_column_to_panel(roSegmen,true);
  
  

  
	//inicia el proceso de eventos
	glutMainLoop();

	return 0;
}


 

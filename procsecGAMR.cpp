/*--------------------------------------------------*/
/*					UCLM - ISA						*/
/*				E.T.S.I. INDUSTRIALES				*/
/*													*/
/*	AvnD. Camilo Jose Cela s/n - 13071 Ciudad Real	*/
/*													*/
/*	Proyecto: Funciones Separables					*/
/*	Investigador Principal: Gloria Bueno Garcia		*/
/*  Programadores: Antonio Martinez					*/
/*													*/
/*	Fecha: Junio - 2005								*/
/*--------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "comun.h"
#include "imagutils.h"
#include "snakes.h"
#include "levelset.h"
#include "GARModel.h"


#define CAB_PARAM_GAR "Dim_X: %d\nDim_Y: %d\nDim_Z: %d\n\nInicial: %d\nFinal: %d\n\nSigma: %f\nIncremento: %f\nNarrow Band: %d\nEpsilon: %f\nThressHold: %f\nInterpolar: %d\nAporte: %f\nFuncion: %d\nValor Region: %d\nSolo Al Inicio: %d\n"
#define THRESS_ITER 0


int Dim_SEC_X;
int Dim_SEC_Y;
int Dim_SEC_Z;
int Inf_SEC;
int Sup_SEC;
float SIGMA_GAUSS;
float INCREMENTO;
int NARROWBAND;
float VELOC_EPSILON;
float THRESS;
int bInterp;
float APORTE;
int FuncReg;
int VALREG;
int bOnlyInitGAR;
int iterac_actual;

/*---------------- SUB MENU -------------------------------*/
void cargarPuntos(char arch[40],lista list){
  FILE *contorno_archivo;
  contorno_archivo = fopen(arch,"r");
  punto2D pto;
  int a,b;
  
  Get_X_Y(&a,&b);
  while(fscanf(contorno_archivo,"%d %d\n",&pto.x,&pto.y)!= EOF){
	pto.x = b - pto.x;
	l_meted(list,&pto);
  }
  fclose(contorno_archivo);
}
/*---------------- SUB MENU -------------------------------*/
void cargar_parametros(char *fich_param){

	FILE *param_archivo= fopen(fich_param,"r");

	fscanf(param_archivo,CAB_PARAM_GAR,&Dim_SEC_X,&Dim_SEC_Y,&Dim_SEC_Z,&Inf_SEC,&Sup_SEC,&SIGMA_GAUSS,&INCREMENTO,&NARROWBAND,&VELOC_EPSILON,
		&THRESS,&bInterp,&APORTE,&FuncReg,&VALREG,&bOnlyInitGAR);
	fclose(param_archivo);
	

}
/*---------------- SUB MENU -------------------------------*/
int main(int argc, char *argv[]){
  //Tiempo Computacional
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  int opcion;
  unsigned char *imag3D;
  unsigned char *imagen;
  unsigned char *imag;
  //unsigned char *imag3D;
  FILE *textresult;
  FILE *param;
  int a,b;
  float increm;
  int nb,n,nsecuencia,iterac;
  lista puntosselec;
  puntosselec=l_nuev(sizeof(punto2D));
  tiempo_total = tiempo1 = tiempo2 = 0;
  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("procsecGAMR <imagen_orginal.raw> <contorno.txt> <parametros.txt>\n");
    
  }else{
	  char secname[40];
	  iterac=0;
	  cargar_parametros(argv[3]);
	  imag3D = (unsigned char*) calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y*Dim_SEC_Z);
	  imagen = (unsigned char*)calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y);
	  imag = (unsigned char*)calloc(sizeof(unsigned char),Dim_SEC_X*Dim_SEC_Y);
	  textresult = fopen("PuntosSecuencia.txt","w");
	  
	  int l;
	  for(l=0;l<40;l++) secname[l] = '\0';
	  sprintf(secname,argv[1]);
	  FILE* archivo;
	  archivo = fopen(secname,"rb");
	  fread((char*)imag3D,sizeof(char),Dim_SEC_X*Dim_SEC_Y*Dim_SEC_Z,archivo);
	  fclose(archivo);
	  char fichname[40];
	  int j_imag = Dim_SEC_X;
	  int i_imag = Dim_SEC_Y;
	  Set_X_Y(j_imag,i_imag);			
	  n = Inf_SEC;
	  nsecuencia = 0;
	  int i,j;
	  for(i=0; i<Dim_SEC_X;i++)		
		  for(j=0; j<Dim_SEC_Y; j++){
			  imagen[calcSingleSubscript(i,j,j_imag)] = imag3D[calcSingleSubscript3D(i,j,n,i_imag,j_imag)];
		  }
	  		 
	  //swappixelOGL(imagen);
	  
	  if (iterac==0){
		  if(l_vacia(puntosselec)){
			  cargarPuntos(argv[2],puntosselec);
		  }
		  InitGAMRModel(imagen,0,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
			  puntosselec,j_imag,i_imag,THRESS,0,0,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
	  }
	  iterac=1;
	  while(n<Sup_SEC){
			char nameres[40];
			if((!StopCriteriorGARModel())&&(!((THRESS_ITER > 0)&&(iterac>THRESS_ITER)))){
				tiempo1 = clock();
				if(/*iterac%NARROWBAND==0*/ EndNarrowDist()){
					ReInitGARModel();
				}
				NextIterGAMRModel(INCREMENTO);
				tiempo2 = clock();
				tiempo_total += tiempo2-tiempo1;
				iterac++;					
			//	PresentGRAModel(imag);	
			}else{
				printf("Fin del procesamiento de la secuencia.\n");
				printf("Numero de Iteraciones: %i\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
				double tiempo = (double)tiempo_total/CLOCKS_PER_SEC;
				tiempo_total = 0;				
				PresentGRAModel(imag);
				sprintf(nameres,"ResultadoGAR_Imagen%d.pgm\0",n);
				Set_X_Y( i_imag,j_imag);
			//	swappixelOGL(imag);
				guardarImagenPGM(imag,nameres);
			//	swappixelOGL(imag);
				
				if(!l_vacia(puntosselec)){
					l_dest(&puntosselec);
					puntosselec = l_nuev(sizeof(punto2D));
				}
			
				SetUpCurveGARModel(puntosselec);
				l_empieza(puntosselec);
				while(l_quedan(puntosselec)){
					punto2D pto;
					l_dame(puntosselec,&pto);
					fprintf(textresult,"(%d, %d, %d)\n",pto.x, pto.y, n);
				}
				DestroyGARModel();
				iterac_actual = iterac;
				iterac=0;
				int i,j,l;
				if(n+1 < Sup_SEC){
				for(i=0; i<Dim_SEC_X;i++){		
					for(j=0; j<Dim_SEC_Y; j++){
						imagen[calcSingleSubscript(i,j,j_imag)] = imag3D[calcSingleSubscript3D(i,j,n,i_imag,j_imag)];
					}
				}
				}else{
//					swappixelOGL(imagen);
				}
				n++;
				nsecuencia++;
//				swappixelOGL(imagen);

				
				if (iterac==0){					
				InitGAMRModel(imagen,0,SIGMA_GAUSS,NARROWBAND,VELOC_EPSILON,
					puntosselec,j_imag,i_imag,THRESS,0,0,bInterp,APORTE,FuncReg,VALREG, bOnlyInitGAR);		
				}
				iterac=1;
			}
		}
	 
			n--;					
		
		   PresentGRAModel(imag);
		   if(!l_vacia(puntosselec)){
			   l_dest(&puntosselec);
			   puntosselec = l_nuev(sizeof(punto2D));
		   }
		   DestroyGARModel();
		   iterac_actual = iterac;
		   iterac=0;					   
		   fclose(textresult);
		   free(imag3D);
		   imag3D = NULL;
		   printf("Fin de la segmentacin.\n");
	
  }
}
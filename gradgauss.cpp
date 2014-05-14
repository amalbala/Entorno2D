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


int main(int argc, char *argv[])
{
  unsigned char *imag;
  double *imag3;
  float num;
  float max, min;
  int a,b;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;

  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("gradgauss <imagen_orginal> <imagen_result> <sigma>\n");
  }else{
  if(argc != 4){
    cout<<"Numero incorrecto de parametros el correctos es 3."<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){

    num = atof(argv[3]);
	if(0<num && num<10){
		Get_X_Y(&a,&b);
   		imag3 = (double*)calloc(sizeof(double),a*b);
		tiempo1 = clock();
		grad_conv_gauss(imag, num, &max,&min,imag3);
		tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
		doubleTounsignedchar(imag3,imag);
		free(imag3);
		guardarImagenPGM(imag,argv[2]);
		free(imag);
	}else{
		cout<<"Parametro fuera de los limites [0,10]"<<endl;
	}
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

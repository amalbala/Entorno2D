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
#include "imagutils.h"

int main(int argc, char *argv[])
{
  unsigned char *imag;
  double *imag3;
  float num;
  int a,b;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;

  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("laplaciano <imagen_orginal> <imagen_result>\n");
  }else{


  if(argc != 3){
    cout<<"Numero incorrecto de parametros el correctos es 2.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
    Get_X_Y(&a,&b);
   	imag3 = (double*)calloc(sizeof(double),a*b);
	tiempo1 = clock();
	laplaciano(imag,imag3);
	tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
	doubleTounsignedchar(imag3,imag);
    free(imag3);
    guardarImagenPGM(imag,argv[2]);
	free(imag);
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

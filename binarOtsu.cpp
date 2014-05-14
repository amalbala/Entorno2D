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
#include "levelset.h"
#include "snakes.h"

int main(int argc, char *argv[])
{
  unsigned char *imag;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;
 

  if(argc == 1){
	 printf(CABECERA_UCLM);
    printf("binarOtsu <imagen_orginal> <imagen_result>\n");
  }else{


  if(argc != 3){
    cout<<"Numero incorrecto de parametros el correctos es 2.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
    int hist[256];
    histograma(imag,hist);
	tiempo1 = clock();
	binarizaImag(imag,umbralOtsu(hist));
	tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
    guardarImagenPGM(imag,argv[2]);
	free(imag);
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }  
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

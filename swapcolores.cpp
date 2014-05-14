
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
#include "imagutils.h"


int main(int argc, char *argv[])
{
  unsigned char *imag;
 

  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("swapcolores <imagen_orginal> <imagen_result>\n");
  }else{


  if(argc != 3){
    cout<<"Numero incorrecto de parametros el correctos es 2.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
	  int X,Y;
	Get_X_Y( &X, &Y);
	for(int i=0; i<X*Y; i++)
		if(imag[i] == 0)
			imag[i] = 255;

    guardarImagenPGM(imag,argv[2]);
	free(imag);
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  return EXIT_SUCCESS;
}

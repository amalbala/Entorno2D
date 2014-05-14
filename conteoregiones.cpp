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

#include <iostream.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#include "comun.h"

#include "imagutils.h"


int main(int argc, char *argv[])
{
  unsigned char *imag;
  FILE *arch; 

  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("conteoregiones <imagen_orginal>\n");
  }else{
  if(argc != 2){
    cout<<"Numero incorrecto de parametros el correctos es 2.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
   extraerSegmentacion(imag);
   free(imag);
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }  
  return EXIT_SUCCESS;
}

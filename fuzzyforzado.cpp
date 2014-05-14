
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
  int num,error;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;
  
  int a,b;

  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("fuzzzyforzado <imagen_orginal> <imagen_result> <archivo_centros> <error>\n");
  }else{


  if(argc != 5){
    cout<<"Numero incorrecto de parametros el correctos es 3.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
	  error = atoi(argv[4]);
	  if(0<error && error<1000){
		Get_X_Y(&a,&b);
		unsigned char *vector;
		cargar_centros_fuzzy(&vector,argv[3],&num);
		tiempo1 = clock();
		fuzzyCmeansforzado(imag,a,b,imag,vector,num,error);   
		tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
		guardarImagenPGM(imag,argv[2]);
		free(imag);	
	  }else{
	  	cout<<"Parametro fuera de los limites [0,1000]"<<endl;
	  }
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

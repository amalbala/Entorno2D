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
  unsigned char *imag_res;
  int num,error;
  int a,b;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;


  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("fuzzycmeans <imagen_orginal> <imagen_result> <num_clusters> <error>\n");
  }else{


  if(argc != 5){
    cout<<"Numero incorrecto de parametros el correctos es 3.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
    Get_X_Y(&a,&b);
    imag_res=(unsigned char *)calloc(sizeof(unsigned char),a*b);
    num = atoi(argv[3]);
	error = atoi(argv[4]);
	if(0<num && num<1000 && 0<error && error<1000){
		tiempo1 = clock();
		fuzzyCmeans(imag,a,b,imag_res,num,error);    
		tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
		guardarImagenPGM(imag_res,argv[2]);
		free(imag_res);
	}else{
		cout<<"Parametros fuera de los limites [0,1000]"<<endl;
	}
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

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


void cargarPuntos(char arch[40],lista list){
  FILE *archivo;
  archivo = fopen(arch,"r");
  punto2D pto;

  while(fscanf(archivo,"%d %d\n",&pto.x,&pto.y)!= EOF){
    l_meted(list,&pto);
  }
  fclose(archivo);
}


int main(int argc, char *argv[])
{
	unsigned char *imag;
	TLlista lGen = new (TCab);
	lista lRes;
	clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;

	printf(CABECERA_UCLM);
	if(argc != 2){
		cout<<"Numero incorrecto de parametros el correctos es 1."<<endl;
		printf("extraer <imagen_orginal>\n");
		exit(0);
	}
	if(cargarImagenPGM(&imag,argv[1])==0){
		tiempo1 = clock();
		extraerSegmentacion(imag,lGen,255);
		tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
		free(imag);
	}else{
		cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
	}
	printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}




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
  unsigned char uNum;
  int iNum;
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  tiempo_total = tiempo1 = tiempo2 = 0;
 
  printf(CABECERA_UCLM);
  if(argc == 1){	
    printf("binarizar <imagen_orginal> <imagen_result> <umbral>\n");
  }else{
  if(argc != 4){
    cout<<"Numero incorrecto de parametros el correctos es 3.\n"<<endl;
    exit(0);
  }
  if(cargarImagenPGM(&imag,argv[1])==0){
    iNum = atoi(argv[3]);
	if(0<iNum && iNum<256){
		uNum = iNum;
		tiempo1 = clock();
		binarizaImag(imag,uNum);
		tiempo2 = clock();
		tiempo_total += tiempo2-tiempo1;
		guardarImagenPGM(imag,argv[2]);
		free(imag);
	}else{
		cout<<"Parametro 3 fuera de los limites [0,255]"<<endl;
	}
  }else{
    cout<<"Error de al abrir fichero. Compruebe que existe"<<endl;
  }
  }
  printf("Tiempo: %f",(double)tiempo_total/CLOCKS_PER_SEC);
  return EXIT_SUCCESS;
}

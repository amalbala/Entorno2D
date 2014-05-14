/*--------------------------------------------------*/
/*					UCLM - ISA						*/
/*				E.T.S.I. INDUSTRIALES				*/
/*													*/
/*	Avdn. Camilo Jose Cela s/n - 13071 Ciudad Real	*/
/*													*/
/*	Proyecto: 3D Target Tool 						*/
/*	Investigador Principal:							*/
/*													*/
/*	Fecha: Diciembre -2004							*/
/*--------------------------------------------------*/

#include <string.h>

using namespace std;

void error(int num,char *msg){
	switch(num){
	case 1:
		strcpy(msg,"Error al abrir fichero.");
	break;
	case 2:
		strcpy(msg,"Error al leer fichero.");
	break;

	}
}

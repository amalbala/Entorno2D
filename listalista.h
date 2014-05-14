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

#ifndef _LISTALISTA_H_
#define _LISTALISTA_H_

#include "listagen.h"

struct TRegion{
	    lista lContorno;
		lista lPuntos;
		TRegion* sig;
};



struct TCab{
	TRegion* prim;
	TRegion* ult;
	TRegion* iter;
	int numElem;
};

//Tipo de Dato Lista
typedef struct TCab *TLlista;

void llInicializa(TLlista l);
void llMete(TLlista l, TRegion r);
void llSaca(TLlista l, TRegion r);
void llDame(TLlista l, TRegion *r);
void llEmpieza(TLlista l);
bool llVacia(TLlista l);
bool llQuedan(TLlista l);

#endif


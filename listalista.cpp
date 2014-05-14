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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "listalista.h"
#include "imagutils.h"

using namespace std;

void llInicializa(TLlista l){
	l->prim = l->iter = l->ult = NULL;
	l->numElem = 0;
}
		


void llMete(TLlista l, TRegion r){
	punto2D pto;
	if(l->numElem == 0){
		l->ult = new(TRegion);
		l->ult->sig = NULL;
		l->prim = l->iter = l->ult;
	}else{
		l->ult->sig = new(TRegion);
		l->ult= l->ult->sig;
	}
	
	l->ult->lContorno = l_nuev(sizeof(punto2D));
	l_empieza(r.lContorno);
	while(l_quedan(r.lContorno)){
		l_dame(r.lContorno,&pto);
		l_meted(l->ult->lContorno,&pto);
	}
	l->ult->lPuntos = l_nuev(sizeof(punto2D));
	l_empieza(r.lPuntos);
	while(l_quedan(r.lPuntos)){
		l_dame(r.lPuntos,&pto);
		l_meted(l->ult->lPuntos,&pto);
	}
	l->ult->sig = NULL;
	l->numElem++;
}


void llSaca(TLlista l, TRegion r){
	l_copia(l->prim->lContorno,r.lContorno);
	l_copia(l->prim->lPuntos,r.lPuntos);
	l_dest(&(l->prim->lContorno));
	l_dest(&(l->prim->lPuntos));
	if(l->numElem != 1){	
		l->iter=l->prim->sig;
		l->prim = l->iter;
		delete(l->prim);
	}else{
		delete(l->prim);
		l->prim = l->iter = l->ult = NULL;
	}
	l->numElem--;
}
void llDame(TLlista l, TRegion *r){
	punto2D pto;
	r->lContorno = l_nuev(sizeof(punto2D));
	r->lPuntos = l_nuev(sizeof(punto2D));
	l_empieza(l->iter->lContorno);
	while(l_quedan(l->iter->lContorno)){
		l_dame(l->iter->lContorno,&pto);
		l_meted(r->lContorno,&pto);
	}
	l_empieza(l->iter->lPuntos);
	while(l_quedan(l->iter->lPuntos)){
		l_dame(l->iter->lPuntos,&pto);
		l_meted(r->lPuntos,&pto);
	}
	l->iter = l->iter->sig;	
}

void llEmpieza(TLlista l){
	l->iter = l->prim;
}

bool llVacia(TLlista l){
	return (l->numElem == 0);
}

bool llQuedan(TLlista l){
	return (l->iter != NULL);
}


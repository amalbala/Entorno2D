/*--------------------------------------------------*/
/*					UCLM - ISA						*/
/*				E.T.S.I. INDUSTRIALES				*/
/*													*/
/*	AvnD. Camilo Jose Cela s/n - 13071 Ciudad Real	*/
/*													*/
/*	Proyecto: Descriptores							*/
/*	Investigador Principal: Gloria Bueno Garcia		*/
/*  Programadores: 									*/
/*													*/
/*	Fecha: sepia - 2006								*/
/*--------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "listagen.h"

using namespace std;

struct l_elem{
	struct l_elem *sig;	
};

struct l_rep{
	struct l_elem* prim;
	struct l_elem* ult;
	struct l_elem* iter;
	size_t tamab;
	int num_elem;
};

lista l_nuev(size_t tbase){
	lista list = (lista) malloc (sizeof(struct l_rep ));
	if (!list) //Si no hay memoria
	{
		fprintf(stderr,"l_nuev :  Error de Memoria.\n");
		exit(0);
	}//fin del if
	
	list->prim=list->ult=list->iter=NULL;
	list->tamab=tbase;
	list->num_elem = 0;
	return(list);
}//fin del metodo

int l_vacia(lista l){
	if (!l)//Si no existe
	{
		fprintf(stderr, "l_vacia: La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l->prim) return(0);
	return(1);
}//fin del metodo

void l_metei(lista l, void *n){
	struct l_elem *nuevo;

	if(!l)//Si la lista no existe
	{
		fprintf(stderr,"l_metei :  La lista no existe.\n");
		exit(3);
	}//fin del if

	//Cojo memoria para el puntero a sig y para el dato
	nuevo=(struct l_elem*)malloc(sizeof(struct l_elem) + l->tamab);
	if(!nuevo)//No hay memoria
	{
		fprintf(stderr,"l_metei :  Error de Memoria.\n");
		exit(0);
	}//fin del if
	
	nuevo->sig=NULL;
	memcpy(nuevo + 1, n, l->tamab);

	if(!l->prim)//Si es el primero 
	{
		l->prim=l->ult=l->iter=nuevo;
	}else{
		nuevo->sig=l->prim;
		l->prim=nuevo;
	}//fin del if-else
	l->num_elem++;

}//fin del metodo

void l_meted(lista l, void *n)
{
	struct l_elem *nuevo;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_meted :  La lista no existe.\n");
		exit(3);
	}//fin del if
	nuevo=(struct l_elem*) malloc(sizeof(struct l_elem)+l->tamab);

	if(!nuevo)//No hay memoria
	{
		fprintf(stderr,"l_meted :  Error de Memoria.\n");
		exit(0);
	}//fin del if
	
	nuevo->sig=NULL;
	memcpy(nuevo+1, n, l->tamab);

	if(!l->prim)//Es el primero
	{
		l->prim=l->ult=l->iter=nuevo;
	}else{
		l->ult->sig=nuevo;
		l->ult=nuevo;
	}//fin del if-else
	l->num_elem++;

}//fin del metodo

void l_sacai(lista l, void *n){
	struct l_elem *aux;
	
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_sacai :  La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l_vacia(l))//La lista esta vacia
	{
		fprintf(stderr,"l_sacai :  Lista Vacia.\n");
		exit(1);
	}//fin del if
	memcpy(n,l->prim+1,l->tamab);
	if(l->prim==l->iter)//El iterador esta al principio
		l->iter=l->prim->sig;

	aux=l->prim;
	if(l->prim==l->ult)//solo hay un elemento
	{
		l->prim=l->iter=l->ult=NULL;
	}else//Hay mas de un elemento
		l->prim=l->prim->sig;
	free(aux);
	l->num_elem--;
}//fin del metodo

void l_sacad(lista l, void *n){
	struct l_elem *aux;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_sacad :  La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l_vacia(l))//Esta vacia la lista
	{
		fprintf(stderr,"l_sacad :  Lista Vacia.\n");
		exit(1);
	}//fin del if
	memcpy(n,l->ult+1,l->tamab);
	aux=l->prim;
	while((aux->sig!=l->ult)&&(aux->sig!=NULL))
	{//Me situa en el nodo anterior al ultimo
		aux=aux->sig;
	}//fin del while
	if(l->iter==l->ult)//el iterador esta al final
		l->iter=aux;
	if(l->prim==l->ult)//si solo hay un elemento
	{
		free(aux);
		l->prim=l->ult=l->iter=NULL;
	}else //Hay mas de un elemento
	{
		l->ult=aux;
		aux=aux->sig;
		free(aux);
		l->ult->sig=NULL;
	}//fin del if-else
	l->num_elem--;
}//fin del metodo

void l_sacax(lista l, void *n, int posic){
	int i;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_sacax :  La lista no existe.\n");
		exit(3);
	}//fin del if

	if(l_vacia(l))//No hay elementos
	{
		fprintf(stderr,"l_sacax :  Lista Vacia.\n");
		exit(1);
	}//fin del if
	if(l->num_elem<posic){
		fprintf(stderr,"l_sacax :  Posicion No Valida.");
		exit(1);
	}//fin del if

	l_empieza(l);
	for(i = 0; i< posic; i++){
		l->iter = l->iter->sig;
	}//fin del for
	memcpy(n,l->iter+1,l->tamab);
			
}//fin del metodo

void l_metex(lista l, void *n, int posic){
	int i;
	struct l_elem *nuevo;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_metex :  La lista no existe.\n");
		exit(3);
	}//fin del if
	nuevo=(struct l_elem*) malloc(sizeof(struct l_elem)+l->tamab);

	if(!nuevo)//No hay memoria
	{
		fprintf(stderr,"l_metex :  Error de Memoria.\n");
		exit(0);
	}//fin del if
	if(l->num_elem<posic){
		fprintf(stderr,"l_metex :  Posicion No Valida.");
		exit(1);
	}//fin del if
	if(posic == 1){
			l_metei(l,n);
	}else{
		if(l->num_elem == posic){
			l_meted(l,n);
		}else{
			l_empieza(l);
			for(i = 0; i< posic; i++)
			{
				l->iter = l->iter->sig;
			}
			memcpy(n,l->iter+1,l->tamab);
		}//fin del if-else
	}//fin del if-else			
}//fin del if

void l_empieza(lista l){
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_empieza :  La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l_vacia(l))//Esta vacia la lista
	{
		fprintf(stderr,"l_empieza:  Lista Vacia.\n");
		exit(1);
	}//fin del if
	l->iter=l->prim;
}//fin del metodo

int l_quedan(lista l){
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_quedan :  La lista no existe.\n");
		exit(3);
	}//fin del if
	return(!(l->iter==NULL));
}//fin del metodo

void l_dame(lista l, void *n){
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_dame :  La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l_vacia(l))//Esta vacia la lista
	{
		fprintf(stderr,"l_dame :  Lista Vacia.\n");
		exit(1);
	}//fin del if
	memcpy(n, l->iter + 1, l->tamab);
	l->iter = l->iter->sig;
}//fin del metodo

void l_dest(lista *l){
	struct l_elem *temp;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_dest :  La lista no existe.\n");
		exit(3);
	}//fin del if
	while((*l)->prim){
		temp=(*l)->prim;
		(*l)->prim=(*l)->prim->sig;
		free(temp);
	}//fin del while
	free(*l);
	(*l)=NULL;
}//fin del metodo

int l_busca(lista l, void *n, int fcmp(const void *a,const void *b)){
	struct l_elem *aux;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_busca :  La lista no existe.\n");
		exit(3);
	}//fin del if
	if(l_vacia(l))//Esta vacia la lista
	{
		fprintf(stderr,"l_busca :  Lista Vacia.\n");
		exit(1);
	}//fin del if
	aux=l->prim;
	while(aux!=NULL){
		if(fcmp(aux+1,n)==0)//Si se ha encontrado
		{
			return(1);
		}//fin del if
		aux=aux->sig;
	}//fin del while
	return(0);
}//fin del metodo

int l_elementos(lista l){
	return l->num_elem;
}//fin del if

void l_copia(lista l1, lista l2){
	void *aux;
	l_empieza(l1);
	aux = malloc(sizeof(l1->tamab));
	while(l_quedan(l1)){
		l_dame(l1,aux);
		l_meted(l2,aux);
	}//fin del while

}//fin del metodo

lista l_copiaD(lista l){
	lista nueval;
	struct l_elem *aux;
	void *temp;
	if(!l)//No existe la lista
	{
		fprintf(stderr,"l_copia :  La lista no existe.\n");
		exit(3);
	}//fin del if

	nueval=l_nuev(l->tamab);
	
	aux=l->prim;
	temp = (void*) malloc (l->tamab);
	while(aux!=NULL){
		memcpy(temp,aux+1,l->tamab);
		l_meted(nueval,&temp);
		aux=aux->sig;
	}//fin del while
	nueval->num_elem = l->num_elem;
	return(nueval);
}//fin del metodo


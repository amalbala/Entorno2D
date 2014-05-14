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


#ifndef LISTAGEN_H
#define LISTAGEN_H
#include <stdio.h>
#include <stdlib.h>

//Tipo de Dato Lista
typedef struct l_rep *lista;


//Función de comparación para poder buscar dentro de una lista
//habrá que implementarla en el caso de que queramos realizar
//busquedas.
int fcmp(const void *a,const void *b); 

//Devuelve una lista diseñada para contener datos del tamaño 
//especificado
lista l_nuev(size_t tbase);
//Devuelve el numero de elementos de la lista
int l_elementos(lista l);
//Devuelve 0 si la lista esta vacia
int l_vacia(lista l);
//Mete un dato por el principio
void l_metei(lista l, void *n);
//Mete un dato por el final
void l_meted(lista l, void *n);
//Saca un dato por el principio
void l_sacai(lista l, void *n);
//Saca un dato por el final
void l_sacad(lista l, void *n);
//Saca un dato en la posicion indicada
void l_sacax(lista l, void *n, int posic);
//Mete un dato en la posicion indicada
void l_metex(lista l, void *n, int posic);
//Inicializa el iterador al principio de la lista
void l_empieza(lista l);
//Indica si el iterador no ha llegado al final
int l_quedan(lista l);
//Devuelve en n una copia del elemento actual al que apunta el
//iterador y aumenta el iterador
void l_dame(lista l, void *n);
//Devuelve una copia de la lista
lista l_copiaD(lista l1);
void l_copia(lista l1, lista l2);
//Destruye la lista
void l_dest(lista *l);
//Busca el elemento en la lista aplicando la funcion de comparacion
int l_busca(lista l, void *n, int fcmp(const void *a,const void *b));

#endif

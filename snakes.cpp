//*--------------------------------------------------*/
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "imagutils.h"
#include "matematics.h"
#include "listagen.h"
#include "comun.h"
#include "snakes.h"


#define FILE_FUZZY "centros.txt"


//Dimensiones de la Imagen
int IMAGEN_HEIGHT;
int IMAGEN_WIDTH;

//Numero de puntos de la Snake

//Umbral para el criterio de parada
#define UMBRAL 1.0

//Matriz W line
double *matrizLine;
//Matriz W edge
double *matrizEdge;
//Matriz W term
double *matrizTerm;
//Matriz de Presion
double *matrizPress;
//Matriz de Distancia
double *matrizDist;
int *matrizEtiq;
//Puntos de la Snake Actual
lista ptosSnake;
//Nuevos puntos de la Snake
lista nuevosptosSnake;
//Iteracion
int num_iter;
//Vector con las energias de los puntos de la snake
double *ult_energia;
//Energia total de los punto para el criterio de parada
double total_energia;
//Energia total de los punto para el criterio de parada
double energia_vib;
//Coeficientes de elasticidad y rigidez
double alpha,beta,teta;
int centro_x, centro_y;
double wt=0.0;
double wd=1.0;
int etiqueta_mayor;
unsigned char *centros_s;

int num_clust;
int cluster_int_s;
/*----------------------------------------------*/
int pertenece_cluster(punto2D p){
  if(((matrizLine[calcSingleSubscript(p.x,p.y,IMAGEN_HEIGHT)]== 130)||
  (matrizLine[calcSingleSubscript(p.x,p.y,IMAGEN_HEIGHT)]== 137))&&
    (matrizEtiq[calcSingleSubscript(p.x,p.y,IMAGEN_HEIGHT)]==etiqueta_mayor))
    return 0;
  return -1;
}
/*----------------------------------------------*/
int pertenece_cluster(){
  l_empieza(ptosSnake);
  punto2D pto;
  while(l_quedan(ptosSnake)){
    l_dame(ptosSnake,&pto);
    if((matrizLine[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)]!= 137)&&
      (matrizLine[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)]!= 130))        
    return -1;
  }
  return 0;
}
/*----------------------------------------------*/
int pertenece_region(){
  l_empieza(ptosSnake);
  punto2D pto;
  int etiq=0;
  while(l_quedan(ptosSnake)){
    l_dame(ptosSnake,&pto);
    if(etiq==0){
      etiq = matrizEtiq[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)];
    }else{
      if(matrizEtiq[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)]!=etiq){
        return -1;
      }
    }
  }
  return 0;
}
 
/*----------------------------------------------*/

//Calculo la Energia Wline * Eline siendo Eline la intensidad
//del pixel en cada punto de la imagen.
void WlineEline(double wl, unsigned char *image){
  //Imagen
  for(int i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++)
    matrizLine[i] = wl * image[i];  
}
/*----------------------------------------------*/
void WlineElineFuzzy(double wl, unsigned char *image){
  //Gradiente convolucionado
  unsigned char *imag_temp;
  int *imag_etiq;
  int i;
  imag_temp = (unsigned char*)malloc(sizeof(unsigned char)*IMAGEN_HEIGHT*IMAGEN_WIDTH);
  imag_etiq = (int*)malloc(sizeof(int)*IMAGEN_HEIGHT*IMAGEN_WIDTH);
  lista contorno;
  contorno = l_nuev(sizeof(punto2D));

  cargar_centros_fuzzy(&centros_s,FILE_FUZZY,&num_clust);         
  fuzzyCmeansforzado(image,IMAGEN_HEIGHT,IMAGEN_WIDTH,imag_temp,centros_s,15,10000);
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++) imag_etiq[i] = imag_temp[i];
  etiq_region(imag_etiq,matrizEtiq);
  cierre_lineal(ptosSnake,l_elementos(ptosSnake), contorno);
  etiqueta_mayor = etiq_mayoritaria(contorno,matrizEtiq);
  /**/printf("Etiqueta mayoritaria: %d\n",etiqueta_mayor);
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++)
    matrizLine[i] = wl * imag_temp[i];
  free(imag_temp);
  free(imag_etiq);
  l_dest(&contorno);
}
/*----------------------------------------------*/
//Calculo la Energia , asociada a los bordes, Wedge:
//        Wedge = we * G(I(x,y)')
//Es el gradiente convolucionado por la gaussiana de la imagen
void WedgeEedge(double we, unsigned char *image){


  int i;
  float max,min;
  grad_conv_gauss(image, 0.5, &max,&min,matrizEdge);
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++){    
    matrizEdge[i] *= we;
  }
}

/*----------------------------------------------*/
//Calcula la Energia de Presion que se rige por:
/*
           Ep = wp * G(I(x,y));

 siendo:
                             |I(x,y) - media|                            
            G(I(x,y)) = 1 - ------------------
                                   K * desv
*/
void WpressEpress(double wp, unsigned char *image){
  double desv;
  int i;
  unsigned int media,temp,tama;
  int *imagen_reg;
  lista contorno;
  punto2D pto;
  
  media = temp = tama = 0;
  
  contorno = l_nuev(sizeof(punto2D));
  cierre_lineal(ptosSnake,l_elementos(ptosSnake),contorno);
  imagen_reg = (int*) calloc (sizeof(int),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  l_empieza(contorno);
  while(l_quedan(contorno)){
    l_dame(contorno,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)] = 1.0;

  }
  determinar_pto_interno(imagen_reg,&pto);
  evol_punto(imagen_reg,pto.x, pto.y);
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH;i++){
    if(imagen_reg[i]==2){
      media += image[i];
      temp += (image[i] * image[i]);
      tama++;
    }
  }
  media /=tama;
  if((((1.0/tama)*temp)-(media*media))<0)
    desv = 0;
  else
    desv = sqrt(((1.0/tama)*temp)-(media*media));
  //Calculo la energia de presion
  if(desv <= 0){
      for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++){
        matrizPress[i] = wp * (1.0-((fabs(image[i] - media))));
      }
  } else {
    for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++){
      matrizPress[i] = wp * (1.0-((abs(image[i] - media))/(2.0*desv)));
    }
  }
  l_dest(&contorno);
  free(imagen_reg);
}                     
/*----------------------------------------------*/
//Calculo la energia asociada a las lineas, Wterm 
//Wterm = wt * (- G(I))
//Es el coeficiento por la gaussiana de la imagen
void WtermEterm(double wt, unsigned char *image){
  double *imag_d = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  int i;
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++)
    imag_d[i] = imagen[i];
  curvatura(imag_d,matrizTerm,grad_diffin);
  for(i=0; i<IMAGEN_HEIGHT*IMAGEN_WIDTH; i++)
    matrizTerm[i]*=wt; 
  free(imag_d);
}
/*----------------------------------------------*/
/*----------------------------------------------*/

//Calculo la siguiente Snake
void siguienteSnake(){
  
  double min;
  int first = 1;
  num_iter++;
  double valorenerg;
  double alpha_v, beta_v;
  double curv;
  //Primer Punto
  int pos,i,j,k,n = 0;
  valorenerg=0;
	first = 1;
  int NUMPTOS;
  //Hago una copia de los puntos para irlos actualizando
  punto2D pto;
  
//Para resolver la ecuacion:
  //                    2                 2
  //        (alpha * |Vs| ) + (beta * |Vss| )
  //Eint = -----------------------------------
  //                      2
  //Aplicando diferencias finitas para calcular las derivadas
  //Asi la energia total de la Snake sera:
  //Etot = Eint + Eimag - Epress = Eint + (Eline + Eedge + Eterm) - Epress
  
  //Primero estudiamos el primer punto, que habra que compararlo con el anterior, el
  //siguiente y consigo mismo
  
  punto2D *vecptosSnake;
  punto2D *vecnuevosptosSnake; 
  vecptosSnake =(punto2D*) malloc (sizeof(punto2D)*l_elementos(ptosSnake));
  vecnuevosptosSnake =(punto2D*) malloc (sizeof(punto2D)*l_elementos(ptosSnake));
  NUMPTOS = l_elementos(ptosSnake);
  l_empieza(ptosSnake);
  for(i=0; i<l_elementos(ptosSnake); i++){
    l_dame(ptosSnake,&pto);
    vecptosSnake[i] = pto;
    vecnuevosptosSnake[i] = pto;
  }
  for(i=0; i<3; i++){
		for(j=0; j<3; j++){
      //Posicion dentro de las matrices
      pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);
      alpha_v = (alpha *(pow(((vecptosSnake[n].x + (i-1)) -
							(vecptosSnake[NUMPTOS-1].x)),2) +
						pow(( (vecptosSnake[n].y + (j-1)) -
							(vecptosSnake[NUMPTOS-1].y )),2)));
      beta_v = beta * (pow( ((vecptosSnake[n+1].x) -
					(2.0 * (vecptosSnake[n].x + (i-1))) +
					(vecptosSnake[NUMPTOS-1].x)),2) +
					(pow( ((vecptosSnake[n+1].y) -
					(2.0 * (vecptosSnake[n].y + (j-1))) +
					(vecptosSnake[NUMPTOS-1].y )),2)));
     pto.x = vecptosSnake[n].x + (i-1);
     pto.y = vecptosSnake[n].y + (j-1);
     curv=curvatura(pto,vecptosSnake[NUMPTOS-1],vecptosSnake[n+1]);
      curv *= -wt;    
    

      valorenerg =
				matrizLine[pos] +
				matrizEdge[pos] +
				matrizTerm[pos] +
        curv +
        ((-teta/2.0) * matrizPress[pos]) +
        ((1.0/2.0)*(alpha_v+beta_v)); 
			if (first == 1){
				min = valorenerg;   
				first = 2;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
			}else{
				if(valorenerg<=min){
					min = valorenerg;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
				}
			}
		}
	}
  //almaceno la energia
	ult_energia[n] = min;
  

 
	//Puntos con antecesor y siguiente
	for(n=1; n<(NUMPTOS-1); n++){
   	first = 1;
		vecnuevosptosSnake[n].x = vecptosSnake[n].x;
		vecnuevosptosSnake[n].y = vecptosSnake[n].y;
   	for(i=0; i<3; i++){
			for(j=0; j<3; j++){
        pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);
        alpha_v = (alpha *(pow(((vecptosSnake[n].x + (i-1)) -
								(vecptosSnake[n-1].x )),2) +
							pow(((vecptosSnake[n].y + (j-1)) -
								(vecptosSnake[n-1].y )),2)));
        beta_v = beta *(pow (((vecptosSnake[n+1].x ) -
						(2.0 * (vecptosSnake[n].x + (i-1))) +
						(vecptosSnake[n-1].x )),2) +
						(pow(( (vecptosSnake[n+1].y ) -
						(2.0 * (vecptosSnake[n].y + (j-1))) +
						(vecptosSnake[n-1].y )),2)));
       pto.x = vecptosSnake[n].x + (i-1);
       pto.y = vecptosSnake[n].y + (j-1);
       curv=curvatura(pto,vecptosSnake[n-1],vecptosSnake[n+1]);
       curv *= -wt;

				valorenerg =
					matrizLine[pos] +
					matrizEdge[pos] +
  				matrizTerm[pos] +
          curv +
          ((-teta/2.0) * matrizPress[pos]) +
					((1.0/2.0)*(alpha_v+beta_v));
					
				if (first == 1){
					min = valorenerg;
					first = 2;
						vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1); 
						vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
				}else{                   
					if(valorenerg <= min){
						min = valorenerg;
						vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
						vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
					}
				}
			}
		}
      //almaceno la energia
		ult_energia[n] = min;
	}  
	//Calculo el ultimo punto

	n = NUMPTOS-1;
	first = 1;
	vecnuevosptosSnake[n].x = vecptosSnake[n].x;
	vecnuevosptosSnake[n].y = vecptosSnake[n].y;
 
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
      pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);

      alpha_v = (alpha * (pow(((vecptosSnake[n].x + (i-1)) -
							(vecptosSnake[n-1].x )),2) +

						pow(((vecptosSnake[n].y + (j-1)) -
							(vecptosSnake[n-1].y )),2)));
        beta_v = beta *(pow(((vecptosSnake[0].x ) -
					(2.0 * (vecptosSnake[n].x + (i-1))) +
					(vecptosSnake[n-1].x )),2) +
					pow(( (vecptosSnake[0].y ) -
					(2.0 * (vecptosSnake[n].y + (j-1))) +
					(vecptosSnake[n-1].y )),2));
       pto.x = vecptosSnake[n].x + (i-1);
       pto.y = vecptosSnake[n].y + (j-1);
       curv=curvatura(pto,vecptosSnake[n-1],vecptosSnake[0]);
       curv *= -wt;

     valorenerg =

				matrizLine[pos] +
				matrizEdge[pos] +
        ((-teta/2.0) * matrizPress[pos]) +
				matrizTerm[pos] +
        curv+
        ((1.0/2.0)*(alpha_v+beta_v));    				
			if (first == 1){
				min = valorenerg;

				first = 2;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);

			}else{
				if(valorenerg <= min){
					min = valorenerg;
					vecnuevosptosSnake[n].x = (vecptosSnake[n].x) + (i-1);
					vecnuevosptosSnake[n].y = (vecptosSnake[n].y) + (j-1);
				}
			}
		}
	}
  
	ult_energia[n] = min;
  
	for(k=0; k<NUMPTOS; k++){
		vecptosSnake[k].x = vecnuevosptosSnake[k].x;
		vecptosSnake[k].y = vecnuevosptosSnake[k].y;

	}
  l_dest(&ptosSnake);  
  ptosSnake = l_nuev(sizeof(punto2D));
  
  for(i=0; i<NUMPTOS; i++){
    pto = vecnuevosptosSnake[i];
    l_meted(ptosSnake,&pto);
  }
  free(vecptosSnake);
  free(vecnuevosptosSnake);  
}
/*----------------------------------------------*/

//Calculo la siguiente Snake
void siguienteSnakeFuzzy(){
  
  double min=0;
  double distan;

  int first = 1;
  num_iter++;
  double valorenerg;
  double alpha_v, beta_v;
  double curv;
  //Primer Punto
  int pos,i,j,k,n = 0;
  valorenerg=0;
	first = 1;
  int NUMPTOS;
  //Hago una copia de los puntos para irlos actualizando
  punto2D pto;
  
//Para resolver la ecuacion:
  //                    2                 2
  //        (alpha * |Vs| ) + (beta * |Vss| )
  //Eint = -----------------------------------
  //                      2
  //Aplicando diferencias finitas para calcular las derivadas
  //Asi la energia total de la Snake sera:
  //Etot = Eint + Eimag - Epress = Eint + (Eline + Eedge + Eterm) - Epress
  
  //Primero estudiamos el primer punto, que habra que compararlo con el anterior, el
  //siguiente y consigo mismo
  
  punto2D *vecptosSnake;
  punto2D *vecnuevosptosSnake; 
  vecptosSnake =(punto2D*) malloc (sizeof(punto2D)*l_elementos(ptosSnake));
  vecnuevosptosSnake =(punto2D*) malloc (sizeof(punto2D)*l_elementos(ptosSnake));
  NUMPTOS = l_elementos(ptosSnake);
  l_empieza(ptosSnake);
  for(i=0; i<l_elementos(ptosSnake); i++){
    l_dame(ptosSnake,&pto);
    vecptosSnake[i] = pto;
    vecnuevosptosSnake[i] = pto;
  }
  if(!pertenece_cluster(vecnuevosptosSnake[n])){
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
      //Posicion dentro de las matrices
      pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);
      alpha_v = (alpha *(pow(((vecptosSnake[n].x + (i-1)) -
							(vecptosSnake[NUMPTOS-1].x)),2) +
						pow(( (vecptosSnake[n].y + (j-1)) -
							(vecptosSnake[NUMPTOS-1].y )),2)));
      beta_v = beta * (pow( ((vecptosSnake[n+1].x) -
					(2.0 * (vecptosSnake[n].x + (i-1))) +
					(vecptosSnake[NUMPTOS-1].x)),2) +
					(pow( ((vecptosSnake[n+1].y) -
					(2.0 * (vecptosSnake[n].y + (j-1))) +
					(vecptosSnake[NUMPTOS-1].y )),2)));
     pto.x = vecptosSnake[n].x + (i-1);
     pto.y = vecptosSnake[n].y + (j-1);
     curv=curvatura(pto,vecptosSnake[NUMPTOS-1],vecptosSnake[n+1]);
      curv *= -wt;    
     distan = sqrt(pow(vecptosSnake[n].x + (i-1) - centro_x,2) +
      pow(vecptosSnake[n].y + (j-1) - centro_y,2));
     distan *= wd;


      valorenerg =
/**/    distan+
				matrizLine[pos] +
				matrizEdge[pos] +
//				matrizTerm[pos] +
        curv +
        ((-teta/2.0) * matrizPress[pos]) +
        ((1.0/2.0)*(alpha_v+beta_v)); 
			if (first == 1){
				min = valorenerg;   
				first = 2;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
			}else{
				if(valorenerg<=min){
					min = valorenerg;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
				}
			}
		}
	}
  //almaceno la energia
	ult_energia[n] = min;
  }

 
	//Puntos con antecesor y siguiente
	for(n=1; n<(NUMPTOS-1); n++){
    
		first = 1;
		vecnuevosptosSnake[n].x = vecptosSnake[n].x;
		vecnuevosptosSnake[n].y = vecptosSnake[n].y;
    if(!pertenece_cluster(vecnuevosptosSnake[n])){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
        pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);
        alpha_v = (alpha *(pow(((vecptosSnake[n].x + (i-1)) -
								(vecptosSnake[n-1].x )),2) +
							pow(((vecptosSnake[n].y + (j-1)) -
								(vecptosSnake[n-1].y )),2)));
        beta_v = beta *(pow (((vecptosSnake[n+1].x ) -
						(2.0 * (vecptosSnake[n].x + (i-1))) +
						(vecptosSnake[n-1].x )),2) +
						(pow(( (vecptosSnake[n+1].y ) -
						(2.0 * (vecptosSnake[n].y + (j-1))) +
						(vecptosSnake[n-1].y )),2)));
       pto.x = vecptosSnake[n].x + (i-1);
       pto.y = vecptosSnake[n].y + (j-1);
       curv=curvatura(pto,vecptosSnake[n-1],vecptosSnake[n+1]);
       distan = sqrt(pow(vecptosSnake[n].x + (i-1) - centro_x,2) +
          pow(vecptosSnake[n].y + (j-1) - centro_y,2));
       distan *= wd;
       curv *= -wt;

				valorenerg =/**/      distan+
					matrizLine[pos] +
					matrizEdge[pos] +
//  				matrizTerm[pos] +
          curv +
          ((-teta/2.0) * matrizPress[pos]) +
					((1.0/2.0)*(alpha_v+beta_v));
					
				if (first == 1){
					min = valorenerg;
					first = 2;
						vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1); 
						vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
				}else{                   
					if(valorenerg <= min){
						min = valorenerg;
						vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
						vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);
					}
				}
			}
		}
   }   //almaceno la energia
		ult_energia[n] = min;
	}  
	//Calculo el ultimo punto

	n = NUMPTOS-1;
	first = 1;
	vecnuevosptosSnake[n].x = vecptosSnake[n].x;
	vecnuevosptosSnake[n].y = vecptosSnake[n].y;
  if(!pertenece_cluster(vecnuevosptosSnake[n])){
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
      pos= calcSingleSubscript((vecptosSnake[n].x)+(i-1),(vecptosSnake[n].y)+(j-1),IMAGEN_HEIGHT);

      alpha_v = (alpha * (pow(((vecptosSnake[n].x + (i-1)) -
							(vecptosSnake[n-1].x )),2) +
						pow(((vecptosSnake[n].y + (j-1)) -
							(vecptosSnake[n-1].y )),2)));
        beta_v = beta *(pow(((vecptosSnake[0].x ) -
					(2.0 * (vecptosSnake[n].x + (i-1))) +
					(vecptosSnake[n-1].x )),2) +
					pow(( (vecptosSnake[0].y ) -
					(2.0 * (vecptosSnake[n].y + (j-1))) +
					(vecptosSnake[n-1].y )),2));
       pto.x = vecptosSnake[n].x + (i-1);
       pto.y = vecptosSnake[n].y + (j-1);
       curv=curvatura(pto,vecptosSnake[n-1],vecptosSnake[0]);
       curv *= -wt;
      distan = sqrt(pow(vecptosSnake[n].x + (i-1) - centro_x,2) +
      pow(vecptosSnake[n].y + (j-1) - centro_y,2));
     distan *= wd;


      valorenerg =
        distan+
				matrizLine[pos] +
				matrizEdge[pos] +
        ((-teta/2.0) * matrizPress[pos]) +
//				matrizTerm[pos] +
        curv+
        ((1.0/2.0)*(alpha_v+beta_v));    				
			if (first == 1){
				min = valorenerg;

				first = 2;
					vecnuevosptosSnake[n].x = vecptosSnake[n].x + (i-1);
					vecnuevosptosSnake[n].y = vecptosSnake[n].y + (j-1);

			}else{
				if(valorenerg <= min){
					min = valorenerg;
					vecnuevosptosSnake[n].x = (vecptosSnake[n].x) + (i-1);
					vecnuevosptosSnake[n].y = (vecptosSnake[n].y) + (j-1);
				}
			}
		}
	}
  }
	ult_energia[n] = min;
  
	for(k=0; k<NUMPTOS; k++){
		vecptosSnake[k].x = vecnuevosptosSnake[k].x;
		vecptosSnake[k].y = vecnuevosptosSnake[k].y;

	}
  l_dest(&ptosSnake);  
  ptosSnake = l_nuev(sizeof(punto2D));
  
  for(i=0; i<NUMPTOS; i++){
    pto = vecnuevosptosSnake[i];
    l_meted(ptosSnake,&pto);
  }
  free(vecptosSnake);
  free(vecnuevosptosSnake);  
}
/*----------------------------------------------*/
//Inicializacion de parametros
void inicializarSnake(unsigned char *imag,lista puntos,int x_imag, int y_imag,
float alpha_s,float beta_s,float teta_s,float wt_s,float wl_s,float we_s,float wp_s,
int fuzzy,int num_c,int cl_int){
  double wl,we,wp;
  //Dimensiones
  IMAGEN_WIDTH = x_imag;
  IMAGEN_HEIGHT = y_imag;
  num_clust = num_c;
  cluster_int_s = cl_int;

  //Matrices de energias
 
  matrizLine = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  matrizEdge = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  matrizTerm = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  matrizPress = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  matrizDist = (double*)calloc(sizeof(double),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  matrizEtiq = (int*)calloc(sizeof(int),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  //Imagen
  imagen = (unsigned char*)calloc(sizeof(unsigned char),IMAGEN_HEIGHT*IMAGEN_WIDTH);
  memcpy(imagen,imag,sizeof(unsigned char)*IMAGEN_HEIGHT*IMAGEN_WIDTH);
   //Inicializacion de los punto iniciales de la snake en forma de circulo
  punto2D pto;
  ptosSnake = l_nuev(sizeof(punto2D));
  ult_energia =(double*)calloc(sizeof(double),l_elementos(puntos));
  l_empieza(puntos);
  while(l_quedan(puntos)){
    l_dame(puntos,&pto);
    l_meted(ptosSnake,&pto);
 }
 //Coeficientes de energias
 //Coeficiente de elasticidad
	alpha = alpha_s;
 //Coeficiente de rigidez
 	beta = beta_s;
 //Coeficiente de presion
  teta = teta_s;
  //Coeficiente Eline
	wl = wl_s;
  //Coeficiente Eedge
	we = we_s;
  //Coeficiente Eterm
	wt = wt_s;
  //Coeficiente de Epress
  wp = wp_s;
  num_iter = 1;    
  
  printf("Valor de los parametros:\n\tAlpha: %f\n\tBeta: %f\n\tTeta: %f\n",alpha,beta,teta);
  printf("\tWline: %f\n\tWedge: %f\n\tWterm: %f\n\tWpress: %f\n",wl,we,wt,wp);
  //Calculo de energias inciales de la snake
  if(fuzzy)
	WlineElineFuzzy(wl,imag);	
  else
	WlineEline(wl,imag);
	WedgeEedge(we,imag);
	WtermEterm(wt,imag);
  WpressEpress(wp,imag); 
  
}

/*----------------------------------------------*/
//Criterio de parada
int parada_snk(){
  double energia_past = total_energia;
  total_energia = 0;
  for(int i=0; i<l_elementos(ptosSnake); i++){
    total_energia += ult_energia[i];    
  }
  if(energia_vib == total_energia) return 0;
  energia_vib = energia_past;
  printf("Energia Vibracion: %f\tEnergia Pasada: %f\t Energia Actual: %f\n",energia_vib,energia_past,total_energia);
  //Si la diferencia entre energia acutal y la anterior es menor que un umbral, paramos
  if((energia_past!= 0)&&((energia_past-total_energia) < UMBRAL)&&
    ((energia_past-total_energia) > - UMBRAL))
    return 0;                                                
  else
    return -1;                                               
}
/*----------------------------------------------*/
void libera_snk(){
	
	free(matrizLine);
	free(matrizEdge);
	free(matrizTerm);
	free(matrizPress);
	free(matrizDist);
	free(matrizEtiq);
	free(imagen);
	free(ult_energia);
	l_dest(&ptosSnake);
}
/*----------------------------------------------*/
void pintar_energias(){
  int i,j,pos;
  FILE* arch;
  arch = fopen("Energias.txt","w");
  for(i=0; i< IMAGEN_HEIGHT; i++){
    for(j=0; j< IMAGEN_WIDTH; j++){
      pos= calcSingleSubscript(i,j,IMAGEN_WIDTH);
    	fprintf(arch,"%d\t",
     (int)(matrizLine[pos]+matrizEdge[pos]+matrizTerm[pos]-
        ((teta/2.0) * matrizPress[pos])));;
    }
    fprintf(arch,"\n");
   }
   fclose(arch);
}
/*----------------------------------------------*/
//Actualiza la imagen sobreimpresionando los puntos de la snake
void pintaguardaSnakeInterp(unsigned char *imag){
  memcpy(imag,imagen,sizeof(unsigned char)*IMAGEN_WIDTH*IMAGEN_HEIGHT);
  punto2D pto;
  lista contorno;
  contorno = l_nuev(sizeof(punto2D));
  cierre_lineal(ptosSnake,l_elementos(ptosSnake), contorno);
  l_empieza(contorno);
  while(l_quedan(contorno)){
    l_dame(contorno,&pto);
		imag[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)] = 255;
     //(-1)*(imagen[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)]-255);
	}
}
/*----------------------------------------------*/
//Actualiza la imagen sobreimpresionando los puntos de la snake
void pintaguardaSnake(unsigned char *imag){
  memcpy(imag,imagen,sizeof(unsigned char)*IMAGEN_WIDTH*IMAGEN_HEIGHT); 
  punto2D pto;
  l_empieza(ptosSnake);
  while(l_quedan(ptosSnake)){
    l_dame(ptosSnake,&pto);
		imag[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)] = 255;
      //(-1)*(imag[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)]-255.0);
	}
}
/*----------------------------------------------*/
//Pinta solo el contorno interpolado de la snake
void guardaSnakeInterp(){
  unsigned char *imag;
  char nombre[20];
  punto2D pto;
  sprintf(nombre,"SnakeEXT_%d",num_iter+100);
	strcat(nombre,".pgm");
  imag = (unsigned char*)calloc(sizeof(unsigned char) ,(IMAGEN_WIDTH * IMAGEN_HEIGHT));
  for(int i=0; i<IMAGEN_WIDTH*IMAGEN_HEIGHT; i++)
      imag[i] = 0.0;

  l_empieza(ptosSnake);
  while(l_quedan(ptosSnake)){
    l_dame(ptosSnake,&pto);
    imag[calcSingleSubscript(pto.x,pto.y,IMAGEN_HEIGHT)] = 255;
  }
  guardarImagenPGM(imag, nombre);
  free(imag);
}

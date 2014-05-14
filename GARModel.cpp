





 #include <math.h>
 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <fstream>
 #include "GARModel.h"
 #include "matematics.h" 
 #include "comun.h"
 #include "listalista.h"

using namespace std;

 //GAR_fa: Factor constante en el caso de la velocidad no asociada a la curvatura
 #define GAR_fa 1.0
 #define APORTATION	0.1

double GAR_THRESS_KI = 0.4;
 
//Minima diferencia entre dos levelset para el criterio de parada
double GAR_MIN_DIF = 0.01;

unsigned char* GAR_Image;
unsigned char* GAR_ImageOri;

//Dimensiones de la Imagen
int GAR_IMAGE_HEIGHT;
int GAR_IMAGE_WIDTH;
int GAR_IMAGE_HEIGHT_old;
int GAR_IMAGE_WIDTH_old;
int GAR_PresentImageWidth;
int GAR_PresentImageHeight;
int GAR_W_imag_ori;
int GAR_H_imag_ori;
int GAR_RegFunction;

//Puntos del Level Set
lista GAR_CurvePoints;
//Level Set anterior
lista GAR_PreviousCurvePoints;
//GAR_NARROW Band
int GAR_NARROW;
//Factor de la curvatura
float GAR_EPSILON;
float GAR_SIGMA_G;
//Matriz con la imagen convolucionada con la gaussiana
double *GAR_GradGaussMtx;
//Matriz del gradiente
double *GAR_GradMtx;
//Matriz del gradiente de la funcion de parada
double *GAR_GradkiMtx;
//Matriz de la velocidad asociada a cada punto de la imagen
double *GAR_VelocMtx;
//Matriz del criterio de parada GAR_kiMtx extendido
double *GAR_kiextMtx;
double *GAR_kiMtx;

//Matrices del levelset
//Matriz de evolucion del level set
double *GAR_fiMtx;
double *GAR_fiOldMtx;
//Seleccion con circulos
bool GAR_circul = false;

//Numero de iteraciones
int GAR_iter=1;
//GAR_modo de evolucion: 0 = expansion, 1=compresion
int GAR_modo;
//Criterio de Parada
float GAR_kant;
//Numero de repeticiones del mismo valor del level set
int GAR_NumRept=0;
//Para localizar un punto interno del level set
punto2D GAR_InternalPto;
punto2D GAR_MinNBPto;
punto2D GAR_MaxNBPto;
punto2D GAR_OldMinNBPto;

//Etiqueta Mayoritaria
int GAR_BestLabel;
int *GAR_LabelMtx;
//Media de intensidades de una region
double fIntMin;
double fDesviation;
int fIntMax;
//Coeficiente de peso entre el modulo de region y el modulo de contorno

double GAR_BETA = 0.5;
int ValorRegion;
lista lROI;


int *GAR_dist;
/*-----------------------------------------------*/

double ValorRegionCL(int intensidad,float distancia){
	double coef;
	if((intensidad < 120) && (intensidad > 0)/*&&(intensidad > 10))*/){
		coef = distancia/(double)intensidad;
	}else{
		return(255);
	}

	return (coef);
}

/*-----------------------------------------------*/
/*-----------------------------------------------*/
/*    LEVEL SET                                  */
/*-----------------------------------------------*/

//Función que determina una nueva Narrow Band
void GAR_NewNB(lista l,int tama_nb,punto2D *pto_min,punto2D *pto_max){
  int x1,x2,y1,y2;
  punto2D pto;

  l_empieza(l);
  l_dame(l,&pto);
  x1 = x2 = pto.x;
  y1 = y2 = pto.y;
  while(l_quedan(l)){
    l_dame(l,&pto);
    if(pto.x < x1)  x1 = pto.x;
    if(pto.x > x2)  x2 = pto.x;
    if(pto.y < y1)  y1 = pto.y;
    if(pto.y > y2)  y2 = pto.y;
  }

  (*pto_min).x = (x1-tama_nb);
  if((*pto_min).x < 0) (*pto_min).x = 0;
  (*pto_min).y = (y1-tama_nb);
  if((*pto_min).y < 0) (*pto_min).y = 0;
  (*pto_max).x = (x2+tama_nb);
  if((*pto_max).x >= GAR_H_imag_ori) (*pto_max).x = GAR_H_imag_ori-1;
  (*pto_max).y = (y2+tama_nb);
  if((*pto_max).y >= GAR_W_imag_ori) (*pto_max).y = GAR_W_imag_ori-1;

}


/*-----------------------------------------------*/
//Genera una imagen binaria en la que el fondo el 0.0 y los puntos del level set 1.0
void GAR_BinaryImage(double *imagen_reg){
  punto2D pto;
  int i;
  for(i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++)
    imagen_reg[i] = 0.0;
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)] = 1.0;
  }

}
/*-----------------------------------------------*/
void GAR_BinaryImage(int *imagen_reg){
  punto2D pto;
  int i;
  for(i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++)
    imagen_reg[i] = 0;
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)] = 1;
  }

}

/*-----------------------------------------------*/
//Crecimiento de regiones para determinar los signos de las zonas del level set
void GAR_FillLabelRegion(){
  int valor,i;
  
  int *imag2;
  int *etiquetas;
  
  imag2 = (int*)calloc(sizeof(int),GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  etiquetas = (int*)calloc(sizeof(int),GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);

  GAR_BinaryImage(imag2);
  etiq_region(imag2,etiquetas);
  free(imag2);
  valor = etiquetas[calcSingleSubscript(0, 0,GAR_IMAGE_WIDTH)];
  if(GAR_modo == 1){
      for(i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
    //Para los puntos localizados cambio el signo.
      if(etiquetas[i] == valor){
        GAR_fiMtx[i] *= -1;
      }
    }
    }else{
      for(i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
      if(etiquetas[i] != valor){
        GAR_fiMtx[i] *= -1;
      }
    }
  }
   free(etiquetas);
  

}
  
  
  
  
/*-----------------------------------------------------------------*/
//Determina si un punto pertenece o no al level set
bool GAR_CurvePoint(int i,int j){
  //Si el valor es menos que cero
  if((GAR_fiMtx[calcSingleSubscript(i,j,GAR_IMAGE_WIDTH)] < 0)&&
  //y alguno de sus 4 vecinos es positivo, entonces es un punto del level set.
        ((GAR_fiMtx[calcSingleSubscript(i-1,j,GAR_IMAGE_WIDTH)] >= 0)||  
         (GAR_fiMtx[calcSingleSubscript(i+1,j,GAR_IMAGE_WIDTH)] >= 0)||
         (GAR_fiMtx[calcSingleSubscript(i,j-1,GAR_IMAGE_WIDTH)] >= 0)||
         (GAR_fiMtx[calcSingleSubscript(i,j+1,GAR_IMAGE_WIDTH)] >= 0)))
          return true;
     return false;
}
/*-----------------------------------------------------------------*/
//Determina los puntos del nuevo level set
void GAR_NewCurve(){
  int i,j;
  punto2D pto;
  //Reinicializamos la lista de puntos del LevelSet
  l_dest(&GAR_CurvePoints);
  GAR_CurvePoints = l_nuev(sizeof(punto2D));
  //para todos los puntos de la imagen determino aquellos que pertenecen al level set
  for(i=1; i<GAR_IMAGE_HEIGHT-1; i++){
    for(j =1; j<GAR_IMAGE_WIDTH-1; j++){
      //Si pertenece lo inserto en la nueva lista
       if(GAR_CurvePoint(i,j)){
           pto.x = i;
           pto.y = j;
           l_meted(GAR_CurvePoints,&pto);

       }
     }
  }

}
/*--------------------------------------------------------------*/
//CALCULO DEL MAPA DE DISTANCIAS
/*--------------------------------------------------------------*/
/*-----------------------------------------------*/
//Calculo el mapa de distancias con signo aplicando un proceso incremental.
void GAR_DistanceMap(){
  int i,j,a,b;
   
  free(GAR_dist);
  //Calculo la Chamfer Distance
  chamfer_distance4v(GAR_fiMtx,GAR_fiMtx);
  //Inicializo la matriz que luego controlará la narrow band proporcional
  GAR_dist = (int*)malloc(sizeof(int) * GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH);
  for(int n=0; n<GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH; n++)
	  GAR_dist[n] = GAR_fiMtx[n];
  //Calculo el signo de las distancias a partir del valor anterior
  if (GAR_modo == 0){
    for(i=0; i<GAR_IMAGE_HEIGHT_old;i++){
        for(j=0; j<GAR_IMAGE_WIDTH_old;j++){
          if(GAR_fiOldMtx[calcSingleSubscript(i,j,GAR_IMAGE_WIDTH_old)]< 0.0){
            a = i+GAR_OldMinNBPto.x-GAR_MinNBPto.x;
            b = j+GAR_OldMinNBPto.y-GAR_MinNBPto.y;
            if(a>=0 && b>=0 && a<GAR_IMAGE_HEIGHT && b<GAR_IMAGE_WIDTH){
              GAR_fiMtx[calcSingleSubscript(a,b,GAR_IMAGE_WIDTH)]=
                -GAR_fiMtx[calcSingleSubscript(a,b,GAR_IMAGE_WIDTH)];
            }
          }
        }
    }
    for(i=0;i<GAR_IMAGE_HEIGHT;i++){
      if(GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)]<0)
        GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)];
      if(GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)]<0)
        GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)];

    }
    for(j=0;j<GAR_IMAGE_WIDTH;j++){
      if(GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)]<0)
        GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)];
      if(GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)]<0)
        GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)];
    }
 
  }else{


    for(i=0; i<GAR_IMAGE_HEIGHT_old;i++){
      for(j=0; j<GAR_IMAGE_WIDTH_old;j++){         
        if(GAR_fiOldMtx[calcSingleSubscript(i,j,GAR_IMAGE_WIDTH_old)]< 0.0){
         a = i+GAR_OldMinNBPto.x-GAR_MinNBPto.x;
         b = j+GAR_OldMinNBPto.y-GAR_MinNBPto.y;
         if(a>=0 && b>=0 && a<GAR_IMAGE_HEIGHT && b<GAR_IMAGE_WIDTH){
           GAR_fiMtx[calcSingleSubscript(a,b,GAR_IMAGE_WIDTH)]=
            -GAR_fiMtx[calcSingleSubscript(a,b,GAR_IMAGE_WIDTH)];
         }
        }
      }
    }
    for(i=0;i<GAR_IMAGE_HEIGHT;i++){
      if(GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)]>0)
        GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(i,0,GAR_IMAGE_WIDTH)];
      if(GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)]>0)
        GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(i,GAR_IMAGE_WIDTH-1,GAR_IMAGE_WIDTH)];

    }
    for(j=0;j<GAR_IMAGE_WIDTH;j++){
      if(GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)]>0)
        GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(0,j,GAR_IMAGE_WIDTH)];
      if(GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)]>0)
        GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)]=
          -GAR_fiMtx[calcSingleSubscript(GAR_IMAGE_HEIGHT-1,j,GAR_IMAGE_WIDTH)];
    } 
  }
  delete(GAR_fiOldMtx);  
}
/*-----------------------------------------------*/
//calculo las distancias basadas en funciones recursivas
void GAR_DistanceMapLabel(){
  //Calculo la Chamfer Distance
  chamfer_distance4v(GAR_fiMtx,GAR_fiMtx);
  //Inicializo la matriz que luego controlará la narrow band proporcional
  GAR_dist = (int*)malloc(sizeof(int) * GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH);
  for(int n=0; n<GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH; n++)
	  GAR_dist[n] = GAR_fiMtx[n];
  //Aplicamos el algoritmo de relleno para determinar los signos
  GAR_FillLabelRegion();

}
/*-----------------------------------------------*/
//Caclculo la velocidad dependiente de la curvatura o GAR_VelocMtx
void GAR_VelocityMapCurvature(){
	curvatura(GAR_fiMtx,GAR_VelocMtx,grad_diffin);
	for(int i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
		GAR_VelocMtx[i] *= (-GAR_EPSILON);
	}
}
/*-----------------------------------*/
//Calculo el factor de parada GAR_kiMtx
void GAR_StopCriterion(){
  int i;
  float a;
  GAR_GradGaussMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  grad_conv_gauss(GAR_Image,GAR_SIGMA_G,&a,&a,GAR_GradGaussMtx);
  for(i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT; i++){
    GAR_kiextMtx[i]= 1.0/(1.0 + GAR_GradGaussMtx[i]);
  }
  free(GAR_GradGaussMtx);
}
/*-----------------------------------*/
//Calculo el factor de parada extendido GAR_kiMtx
void GAR_StopCriterionExtension(){
  int i,j;
  double *level;
  punto2D *dist_4sed;
  punto2D pto;
  int tama =GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH;
  float a;
  double v_grad;
  level = (double*)calloc(sizeof(double),tama);
  
  GAR_GradGaussMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  grad_conv_gauss(GAR_Image,GAR_SIGMA_G,&a,&a,GAR_GradGaussMtx);
  for(i=0;i<tama;i++)
    GAR_kiextMtx[i] = 0.0;
  l_empieza(GAR_CurvePoints);
  //La velocidad para los puntos del Level set
   while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    //El valor es igual a:
    /*
                  1
    GAR_kiMtx =  ------------------
          (1 + |G * I(x,y)|)
    */
    v_grad = 1.0/(1.0 + GAR_GradGaussMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]);
    if(v_grad<GAR_THRESS_KI)
      GAR_kiextMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= 0.0;
    else
      GAR_kiextMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= v_grad;
      
    level[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= 1.0;
  }
   //Para hallar los puntos de level set mas cercano a cada punto de la imagen
  //aplico un mapa de distancias 4SSED
   free(GAR_GradGaussMtx);
   dist_4sed = (punto2D*)calloc(sizeof(punto2D),tama);
   vecin4ssed_distance(level,dist_4sed);
   //Los puntos que no pertenecen al level set tienen una valor igual a la velocidad del
  //punto de la level set mas cercana a ese punto.
  int pos;
  for(i=1;i<GAR_IMAGE_HEIGHT-1;i++){
    for(j=1;j<GAR_IMAGE_WIDTH-1;j++){
        pos= calcSingleSubscript(i,j,GAR_IMAGE_WIDTH);
        GAR_kiextMtx[pos] =
          GAR_kiextMtx[calcSingleSubscript(i - dist_4sed[pos].x,j - dist_4sed[pos].y,GAR_IMAGE_WIDTH)];
    }
  }

  free(level);
  free(dist_4sed);
}
/*-----------------------------------------------*/
void ReInitGARModel(){
  punto2D pto;
  //Almacenamos la matriz actual de GAR_fiMtx
  
  //free(GAR_GradGaussMtx);
  //Matriz del gradiente
  free(GAR_GradMtx);
  free(GAR_GradkiMtx);
  //Matriz de la velocidad asociada a cada punto de la imagen
  free(GAR_VelocMtx);
  free(GAR_kiextMtx);
  free(GAR_Image);
  
  //Guarda el level set antiguo 
  lista GAR_CurvePoints_old;
  GAR_CurvePoints_old = l_nuev(sizeof(punto2D));
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
     l_dame(GAR_CurvePoints, &pto);
     pto.x += GAR_MinNBPto.x;

     pto.y += GAR_MinNBPto.y;
     l_meted(GAR_CurvePoints_old,&pto);
   }
  l_dest(&GAR_CurvePoints);
  GAR_CurvePoints = l_nuev(sizeof(punto2D));
  //Guarda los valores antiguos de Fi y de la Narrow Band
  GAR_OldMinNBPto=GAR_MinNBPto;  
  GAR_fiOldMtx=(double*)malloc(sizeof(double)*(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  memcpy(GAR_fiOldMtx,GAR_fiMtx,(sizeof(double)*GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  free(GAR_fiMtx);
  //Determina una nueva narrow band
  GAR_NewNB(GAR_CurvePoints_old,GAR_NARROW,&GAR_MinNBPto,&GAR_MaxNBPto);
  //Actualiza los puntos del level set a la nueva narrow band
  l_empieza(GAR_CurvePoints_old);
  while(l_quedan(GAR_CurvePoints_old)){
     l_dame(GAR_CurvePoints_old, &pto);
     pto.x -= GAR_MinNBPto.x;
     pto.y -= GAR_MinNBPto.y;
     l_meted(GAR_CurvePoints,&pto);
   }
  l_dest(&GAR_CurvePoints_old);

  GAR_IMAGE_HEIGHT_old=GAR_IMAGE_HEIGHT;
  GAR_IMAGE_WIDTH_old=GAR_IMAGE_WIDTH;
  GAR_IMAGE_HEIGHT = GAR_MaxNBPto.x -GAR_MinNBPto.x;
  GAR_IMAGE_WIDTH = GAR_MaxNBPto.y -GAR_MinNBPto.y;
  GAR_Image = (unsigned char*)calloc(sizeof(unsigned char),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  int i,j,k=0;
  for(i=GAR_MinNBPto.x; i<GAR_MaxNBPto.x; i++){
    for(j=GAR_MinNBPto.y; j<GAR_MaxNBPto.y; j++){
     GAR_Image[k]=GAR_ImageOri[calcSingleSubscript(i,j,GAR_W_imag_ori)];
     k++;
    }
  }
  
  Set_X_Y(GAR_IMAGE_WIDTH,GAR_IMAGE_HEIGHT);
  //Matriz del gradiente
  GAR_GradMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz del gradiente del criterio de parada
  GAR_GradkiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz de la velocidad asociada a cada punto de la imagen
  GAR_VelocMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz del criterio de parada extendido
  GAR_kiextMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz del levelset
  GAR_fiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
 
  l_empieza(GAR_CurvePoints);
   while(l_quedan(GAR_CurvePoints)){
     l_dame(GAR_CurvePoints, &pto);
     GAR_fiMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= 1.0;
   }
  //Calculamos la matriz de distancias con signo aplicando recursion
  GAR_DistanceMap();  
  //Calculamos la matriz de velocidad
  GAR_VelocityMapCurvature();
  //Calculamos la matriz asociada al criterio de parada
  GAR_StopCriterionExtension();
}
/*-----------------------------------------------------------------*/
//Filtra puntos repetidos en la lista del level set
void GAR_CurveFilter(){
  punto2D pto,ptoant;
  l_empieza(GAR_CurvePoints);
  lista filtrols;
  filtrols = l_nuev(sizeof(punto2D));
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    l_meted(filtrols,&pto);
  }

  l_dest(&GAR_CurvePoints);
  GAR_CurvePoints = l_nuev(sizeof(punto2D));
  l_empieza(filtrols);
  ptoant.x = ptoant.y = 0;
  while(l_quedan(filtrols)){
    l_dame(filtrols,&pto);
    if((pto.x>0)&&(pto.y>0)&&((pto.x != ptoant.x)||(pto.y != ptoant.y)))
      l_meted(GAR_CurvePoints,&pto);
    ptoant.x = pto.x;
    ptoant.y = pto.y;
  }
  l_dest(&filtrols);
}
/*-----------------------------------------------------------------*/
//Guarda el level set anterior
void GAR_StorePreviousCurve(){
  punto2D pto;
  l_empieza(GAR_CurvePoints);
  l_dest(&GAR_PreviousCurvePoints);
  GAR_PreviousCurvePoints= l_nuev(sizeof(punto2D));
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    l_meted(GAR_PreviousCurvePoints,&pto);
  }
}
/*-----------------------------------------------------------------*/
//Calculo del Termino de la Region 
double GAR_RegionTerm(int i){
	double dDifMin, dDifMax;

	//Establece la distancia de una intensidad de la imagen al intervalo de interés, estudiando
	//la distancia a los extremo del intervalo.
	dDifMin = GAR_Image[i] - fIntMin;
	dDifMax = GAR_Image[i] - fIntMax;
	//El valor será la distancia euclidea a ambos extremos
	return (sqrt((dDifMin * dDifMin)+(dDifMax * dDifMax)));	

}
/*-----------------------------------------------------------------*/
//Calculo del Teremino del Contorno
double GAR_BoundaryTerm(int i){
	//Es la funcion tipica de las geodesicas
	//(g(1.0 + K) * N) - (g · N)
	return ((GAR_kiextMtx[i]*(GAR_fa + GAR_VelocMtx[i])*GAR_GradMtx[i]) - 
		(GAR_GradkiMtx[i]*GAR_GradMtx[i]));		
}
/*-----------------------------------------------------------------*/
//Calcula la media de intensidad de la region contenida por el level set
double RegionIntensityAverage(){
	double dTotal=0;
	int iNum=0;
	for(int i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
		if(GAR_fiMtx[i]<0){
			iNum++;
			dTotal += GAR_Image[i];
		}
	}
	return(dTotal/iNum);
}
/*-----------------------------------------------------------------*/
/*
void DeterminarRegionesInteres(lista lRegiones){
	
	numPicos(lRegiones,3,10,GAR_);
}
*/


/*-----------------------------------------------------------------*/
//Calcula el máximo del histograma de la region contenida por el level set
int RegionIntensityMax(){
	int vHistRegion[255];
	int i;
	
	for(i=0; i<255; i++) vHistRegion[i]= 0;

	for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
		if(GAR_fiMtx[i]<0){
			vHistRegion[GAR_Image[i]]++;	
		}
	}
	int iMax = 0;
	int iIntensity;
	for(i=0;i<255; i++){
		if(vHistRegion[i] > iMax){
			iMax = vHistRegion[i];
			iIntensity = i;
		}
	}
	return iIntensity;	
}
/*-----------------------------------------------------------------*/
//Calcula la desviación estandar de la region contenida por le level set
double RegionStandardDesviation(double iAverage){	
	int iNum=0;
	double dSum = 0;
	if(iAverage==0){
		iAverage = RegionIntensityAverage();		
	}
	for(int i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){
			if(GAR_fiMtx[i]<0){
				iNum++;
				dSum += iAverage + GAR_Image[i];
			}
	}
	return(sqrt(dSum/iNum));
	
}
/*-----------------------------------------------------------------*/
//Siguiente iteracion del modelo
void NextIterGARModel(float inct){
  float valor,nuevovalor,media;
  int i;
  int nReg = 0;
  
  //Mapa de Velocidades seguna la curvatura
  GAR_VelocityMapCurvature();
  //Mapa del Criterio de Parada
  GAR_StopCriterionExtension();
  //Calculo del gradiente de Fi
  grad_diffin(GAR_fiMtx,GAR_GradMtx);
   //Calculo el gradiente de Criterio de Parada
  grad_diffin(GAR_kiextMtx,GAR_GradkiMtx);
  
  switch(GAR_RegFunction){
	  //Segun el valor, se tomará como centro de intervalo de interes
	  //la media de la región o el máximo del histograma de la region.
	  case 0:
	ValorRegion = RegionIntensityAverage();
	break;
  case 1:
	  ValorRegion = RegionIntensityMax();
	  break;
  }


  //Estudiamos la primera region de interés que será la referncia de las demas
  //Calculo la desviacion estandar que determina la amplitud del intervalo
  fDesviation = RegionStandardDesviation(ValorRegion);
  fIntMin = ValorRegion - fDesviation;
  fIntMax = ValorRegion + fDesviation;
  
  
  //Matrices que contendran el valor de los terminos de contorno y region
  double *BoundaryMtx,*RegionMtx;
  BoundaryMtx = (double*)malloc(sizeof(double)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  RegionMtx = (double*)malloc(sizeof(double)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
 
  //Valores para el reescalado de las matrices
  double MaxBoundary, MinBoundary;
  double MaxRegion, MinRegion;
 
  MaxBoundary = MinBoundary = GAR_BoundaryTerm(0);
  MaxRegion = MinRegion = GAR_RegionTerm(0);
  
  
  //Vamos calculando los valores de la matrices y acutalizando los maximos y minimos 
  //para poder reescalar
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
	  BoundaryMtx[i] = GAR_BoundaryTerm(i);
	  if(BoundaryMtx[i] > MaxBoundary) MaxBoundary = BoundaryMtx[i];
	  if(BoundaryMtx[i] < MinBoundary) MinBoundary = BoundaryMtx[i];
	  
	  RegionMtx[i] = GAR_RegionTerm(i);
	  if(RegionMtx[i] > MaxRegion) MaxRegion = RegionMtx[i];
	  if(RegionMtx[i] < MinRegion) MinRegion = RegionMtx[i];
  }
    
  //Normalización
  double BoundaryInterv = MaxBoundary - MinBoundary;
  double RegionInterv = MaxRegion - MinRegion;
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
	  BoundaryMtx[i] = (BoundaryMtx[i] - MinBoundary)/BoundaryInterv;
	  RegionMtx[i] = (RegionMtx[i] - MinRegion)/RegionInterv;
  }

  //Ecuación general de la GAR donde se calcula el nuevo valor de Fi teniendo en 
  //cuenta ambos terminos
  // Fi_nuevo = Fi_anterior - t * ((Aporte * TerminoContorno) - ((1-Aporte) * TerminoRegion))
  
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
		  GAR_fiMtx[i] = GAR_fiMtx[i] - (inct *((GAR_BETA*BoundaryMtx[i]) - ((1.0 - GAR_BETA)*RegionMtx[i])));
	}
	
  free(BoundaryMtx);
  free(RegionMtx);
   
  GAR_NewCurve();
  GAR_iter++;
  
}
/*-----------------------------------------------------------------*/
//Siguiente iteracion del modelo
void NextIterGAMRModel(float inct){
  float valor,nuevovalor,media;
  int i;
  int nReg = 0;

  /*
  unsigned char* temp = (unsigned char*)malloc(sizeof(unsigned char)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  char nombreI[40];
  */
  
  //Mapa de Velocidades seguna la curvatura
  GAR_VelocityMapCurvature();
  //Mapa del Criterio de Parada
  GAR_StopCriterionExtension();
  //Calculo del gradiente de Fi
  grad_diffin(GAR_fiMtx,GAR_GradMtx);
   //Calculo el gradiente de Criterio de Parada
  grad_diffin(GAR_kiextMtx,GAR_GradkiMtx);

  //Estudiamos la primera region de interés que será la referncia de las demas
  //Calculo la desviacion estandar que determina la amplitud del intervalo

  switch(GAR_RegFunction){
	  //Segun el valor, se tomará como centro de intervalo de interes
	  //la media de la región o el máximo del histograma de la region.
	  case 0:
		  if(lROI){
			   l_dest(&lROI);
		  }
	    lROI = l_nuev(sizeof(int));
	    RegionesInteres(lROI,2,0.2,GAR_Image);
	  break;
  }


  l_empieza(lROI);
  
	  l_dame(lROI,&ValorRegion);
	  //fDesviation = RegionStandardDesviation(ValorRegion);
	  /**/printf("Valor Region %d: %d\n",nReg,ValorRegion);
	  nReg++;
	  fDesviation = 3;
    
  fIntMin = ValorRegion - fDesviation;
  fIntMax = ValorRegion + fDesviation;
  
  
  //Matrices que contendran el valor de los terminos de contorno y region
  double *BoundaryMtx,*RegionMtx;
  BoundaryMtx = (double*)malloc(sizeof(double)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  RegionMtx = (double*)malloc(sizeof(double)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
 
  //Valores para el reescalado de las matrices
  double MaxBoundary, MinBoundary;
  double MaxRegion, MinRegion;
 
  MaxBoundary = MinBoundary = GAR_BoundaryTerm(0);
  //MaxRegion = MinRegion = ValorRegionCL(GAR_Image[0],GAR_RegionTerm(0));
  MaxRegion = MinRegion = GAR_RegionTerm(0);
  
  
  //Vamos calculando los valores de la matrices y acutalizando los maximos y minimos 
  //para poder reescalar
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	 
	  BoundaryMtx[i] =GAR_BoundaryTerm(i);
	  	  if(BoundaryMtx[i] > MaxBoundary) MaxBoundary = BoundaryMtx[i];
	  if(BoundaryMtx[i] < MinBoundary) MinBoundary = BoundaryMtx[i];
	  
	  RegionMtx[i] = GAR_RegionTerm(i);
	  //RegionMtx[i] = ValorRegionCL(GAR_Image[i],GAR_RegionTerm(i));
	  if(RegionMtx[i] > MaxRegion) MaxRegion = RegionMtx[i];
	  if(RegionMtx[i] < MinRegion) MinRegion = RegionMtx[i];
  }
  /*
  sprintf(nombreI,"MatrizReg_%da_%d.txt",nReg,GAR_iter);
  guardarImagenTXT(RegionMtx, nombreI);  
  */
  //Normalización
  double BoundaryInterv = MaxBoundary - MinBoundary;
  double RegionInterv = MaxRegion - MinRegion;
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
	  BoundaryMtx[i] = (BoundaryMtx[i] - MinBoundary)/BoundaryInterv;
	  RegionMtx[i] = (RegionMtx[i] - MinRegion)/RegionInterv;
  }

  
  
  /*
  sprintf(nombreI,"MatrizReg_%db_%d.txt",nReg,GAR_iter);
  guardarImagenTXT(RegionMtx, nombreI);
  */
  //Vamos estudiando las distintas regiones de interes y actualizando la matriz del GAR
  double *RegionMtx2;
  RegionMtx2 = (double*)calloc(sizeof(double),GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  //l_empieza(lROI);
  while(l_quedan(lROI)){
	  l_dame(lROI,&ValorRegion);
	  /**/printf("Valor Region %d: %d\n",nReg,ValorRegion);
	  nReg++;
	  
	  fIntMin = ValorRegion - fDesviation;
	  fIntMax = ValorRegion + fDesviation;
	  MaxRegion = MinRegion = GAR_RegionTerm(0);
	  //MaxRegion = MinRegion = ValorRegionCL(GAR_Image[0],GAR_RegionTerm(0));

		  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
			  RegionMtx2[i] = GAR_RegionTerm(i);
			  //RegionMtx2[i] = ValorRegionCL(GAR_Image[i],GAR_RegionTerm(i));
			  if(RegionMtx2[i] > MaxRegion) MaxRegion = RegionMtx2[i];
			  if(RegionMtx2[i] < MinRegion) MinRegion = RegionMtx2[i];
		  }
		  double RegionInterv = MaxRegion - MinRegion;
		  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
			  RegionMtx2[i] = (RegionMtx2[i] - MinRegion)/RegionInterv;
		  }
	  
	  /**/
	  //Intentamos que se acerque a ambas regiones
	  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  

		  if(RegionMtx2[i]<RegionMtx[i]) RegionMtx[i] = RegionMtx2[i];
	  }
/*	
  sprintf(nombreI,"MatrizReg_%d_%d.txt",nReg,GAR_iter);
  guardarImagenTXT(RegionMtx, nombreI);
  */
  
  }
 
  //Ecuación general de la GAR donde se calcula el nuevo valor de Fi teniendo en 
  //cuenta ambos terminos
  // Fi_nuevo = Fi_anterior - t * ((Aporte * TerminoContorno) - ((1-Aporte) * TerminoRegion))
  
  for(i=0;i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT;i++){	  
		  GAR_fiMtx[i] = GAR_fiMtx[i] - (inct *((GAR_BETA*BoundaryMtx[i]) - ((1.0 - GAR_BETA)*RegionMtx[i])));
	}
	
  
  
  //free(temp);
  
 free(BoundaryMtx);
  free(RegionMtx);
  free(RegionMtx2);
  
  GAR_NewCurve();
  GAR_iter++;
  

  
}

/*-----------------------------------------------------------------*/
//Inicialización del modelo
void InitGARModel(unsigned char* imagen_reg,int mod,float sigma,int NARROW,float EPSILON,
				  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,int bInterpol, 
				  float weight,int RegionFunc, float ValorReg, int OnlyInit){
  punto2D pto;
  int i,j,k;
  GAR_modo = mod;
  GAR_NARROW = NARROW;
  GAR_EPSILON = EPSILON;
  GAR_MIN_DIF = umbral;
  GAR_SIGMA_G = sigma;
  GAR_THRESS_KI = th_ki;
  GAR_W_imag_ori = width_imag;
  GAR_H_imag_ori = height_imag;
  GAR_iter = 1;
  GAR_NumRept=0;
  GAR_BETA = weight;
  
  //Determina la nueva narrow band  
  GAR_NewNB(selec,GAR_NARROW,&GAR_MinNBPto,&GAR_MaxNBPto);
  GAR_IMAGE_HEIGHT = GAR_MaxNBPto.x -GAR_MinNBPto.x;
  GAR_IMAGE_WIDTH = GAR_MaxNBPto.y -GAR_MinNBPto.y;
  //Inicializa matrices
  GAR_Image = (unsigned char*)calloc(sizeof(unsigned char),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  GAR_ImageOri = (unsigned char*)calloc(sizeof(unsigned char),(GAR_W_imag_ori*GAR_H_imag_ori));
  k=0;
  memcpy(GAR_ImageOri,imagen_reg,sizeof(unsigned char)*GAR_W_imag_ori*GAR_H_imag_ori);
  for(i=GAR_MinNBPto.x; i<GAR_MaxNBPto.x; i++){
    for(j=GAR_MinNBPto.y; j<GAR_MaxNBPto.y; j++){
     GAR_Image[k]=imagen_reg[calcSingleSubscript(i,j,GAR_W_imag_ori)];
     k++;
    }
  }
  Set_X_Y(GAR_IMAGE_WIDTH,GAR_IMAGE_HEIGHT);
  
  //Matriz del gradiente
  GAR_GradMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //MAtriz del gradiente del criterio de parada Ki
  GAR_GradkiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz de la velocidad asociada a cada punto de la imagen
  GAR_VelocMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz del criterio de parada
  GAR_kiextMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matrices del levelset
  GAR_fiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));

  GAR_CurvePoints = l_nuev(sizeof(punto2D));
  lista pto_temp= l_nuev(sizeof(punto2D));
  //Interpolación de la entrada si no son circulos
  if(!circ){
  //Como puede que no sean puntos 4 conexos interpolamos
	  if(bInterpol){
		  cierre_lineal(selec,l_elementos(selec),pto_temp);
		  l_empieza(pto_temp);
		  while(l_quedan(pto_temp)){
			  l_dame(pto_temp,&pto);
			  pto.x -= GAR_MinNBPto.x;
			  pto.y -= GAR_MinNBPto.y;
			  l_meted(GAR_CurvePoints,&pto);
		  }
		  
	  }else{
		  l_empieza(selec);
		  while(l_quedan(selec)){
			  l_dame(selec,&pto);
			  pto.x -= GAR_MinNBPto.x;
			  pto.y -= GAR_MinNBPto.y;
			  l_meted(GAR_CurvePoints,&pto);
		  }
	  }
   //Filtramos si hay puntos repetidos
	GAR_CurveFilter();
   }else{
     GAR_circul = true;
     l_empieza(selec);
	 while(l_quedan(selec)){
		 l_dame(selec, &pto);
		 pto.x -= GAR_MinNBPto.x;
		 pto.y -= GAR_MinNBPto.y;
		 l_meted(GAR_CurvePoints,&pto);
	 }
   }

   l_dest(&pto_temp);
   l_empieza(GAR_CurvePoints);
   //Almacena los puntos del level set 
   while(l_quedan(GAR_CurvePoints)){
     l_dame(GAR_CurvePoints, &pto);
     GAR_fiMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= 1.0;
   }
  //Calculamos la matriz de distancias con signo aplicando recursion
  GAR_DistanceMapLabel();
  //Calculamos la matriz de velocidad
  GAR_VelocityMapCurvature();
  //Calculamos la matriz asociada al criterio de parada
  GAR_StopCriterionExtension();
  //El tipo de función de region seleccionada
  GAR_RegFunction = RegionFunc;
    //Si solo queremos que se calcule en la primera iteracion

  if(OnlyInit){
	   switch(GAR_RegFunction){
		   case 0:
			   ValorRegion = RegionIntensityAverage();

		   break;
		   case 1:
			   ValorRegion = RegionIntensityMax();
		   break;
	   }

	  GAR_RegFunction = -1;
  }

}

//Inicialización del modelo
void InitGAMRModel(unsigned char* imagen_reg,int mod,float sigma,int NARROW,float EPSILON,
				  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,int bInterpol, 
				  float weight,int RegionFunc, float ValorReg, int OnlyInit){
  punto2D pto;
  int i,j,k;
  GAR_modo = mod;
  GAR_NARROW = NARROW;
  GAR_EPSILON = EPSILON;
  GAR_MIN_DIF = umbral;
  GAR_SIGMA_G = sigma;
  GAR_THRESS_KI = th_ki;
  GAR_W_imag_ori = width_imag;
  GAR_H_imag_ori = height_imag;
  GAR_iter = 1;
  GAR_NumRept=0;
  GAR_BETA = weight;
  
  //Determina la nueva narrow band  
  GAR_NewNB(selec,GAR_NARROW,&GAR_MinNBPto,&GAR_MaxNBPto);
  GAR_IMAGE_HEIGHT = GAR_MaxNBPto.x -GAR_MinNBPto.x;
  GAR_IMAGE_WIDTH = GAR_MaxNBPto.y -GAR_MinNBPto.y;
  //Inicializa matrices
  GAR_Image = (unsigned char*)calloc(sizeof(unsigned char),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  GAR_ImageOri = (unsigned char*)calloc(sizeof(unsigned char),(GAR_W_imag_ori*GAR_H_imag_ori));
  k=0;
  memcpy(GAR_ImageOri,imagen_reg,sizeof(unsigned char)*GAR_W_imag_ori*GAR_H_imag_ori);
  for(i=GAR_MinNBPto.x; i<GAR_MaxNBPto.x; i++){
    for(j=GAR_MinNBPto.y; j<GAR_MaxNBPto.y; j++){
     GAR_Image[k]=imagen_reg[calcSingleSubscript(i,j,GAR_W_imag_ori)];
     k++;
    }
  }
  Set_X_Y(GAR_IMAGE_WIDTH,GAR_IMAGE_HEIGHT);
  
  //Matriz del gradiente
  GAR_GradMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //MAtriz del gradiente del criterio de parada Ki
  GAR_GradkiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz de la velocidad asociada a cada punto de la imagen
  GAR_VelocMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matriz del criterio de parada
  GAR_kiextMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));
  //Matrices del levelset
  GAR_fiMtx=(double*)calloc(sizeof(double),(GAR_IMAGE_HEIGHT*GAR_IMAGE_WIDTH));

  GAR_CurvePoints = l_nuev(sizeof(punto2D));
  lista pto_temp= l_nuev(sizeof(punto2D));
  //Interpolación de la entrada si no son circulos
  if(!circ){
  //Como puede que no sean puntos 4 conexos interpolamos
	  if(bInterpol){
		  cierre_lineal(selec,l_elementos(selec),pto_temp);
		  l_empieza(pto_temp);
		  while(l_quedan(pto_temp)){
			  l_dame(pto_temp,&pto);
			  pto.x -= GAR_MinNBPto.x;
			  pto.y -= GAR_MinNBPto.y;
			  l_meted(GAR_CurvePoints,&pto);
		  }
		  
	  }else{
		  l_empieza(selec);
		  while(l_quedan(selec)){
			  l_dame(selec,&pto);
			  pto.x -= GAR_MinNBPto.x;
			  pto.y -= GAR_MinNBPto.y;
			  l_meted(GAR_CurvePoints,&pto);
		  }
	  }
   //Filtramos si hay puntos repetidos
	GAR_CurveFilter();
   }else{
     GAR_circul = true;
     l_empieza(selec);
	 while(l_quedan(selec)){
		 l_dame(selec, &pto);
		 pto.x -= GAR_MinNBPto.x;
		 pto.y -= GAR_MinNBPto.y;
		 l_meted(GAR_CurvePoints,&pto);
	 }
   }

   l_dest(&pto_temp);
   l_empieza(GAR_CurvePoints);
   //Almacena los puntos del level set 
   while(l_quedan(GAR_CurvePoints)){
     l_dame(GAR_CurvePoints, &pto);
     GAR_fiMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)]= 1.0;
   }
  //Calculamos la matriz de distancias con signo aplicando recursion
  GAR_DistanceMapLabel();
  //Calculamos la matriz de velocidad
  GAR_VelocityMapCurvature();
  //Calculamos la matriz asociada al criterio de parada
  GAR_StopCriterionExtension();
  //El tipo de función de region seleccionada
  GAR_RegFunction = RegionFunc;
    //Si solo queremos que se calcule en la primera iteracion

  if(OnlyInit){
/*	   switch(GAR_RegFunction){
		   case 0:
			   vIntRegion[0] = RegionIntensityAverage();
			   //ValorRegion = RegionIntensityAverage();

		   break;
		   case 1:
			   vIntRegion[0] = RegionIntensityMax();
			   //ValorRegion = RegionIntensityMax();
		   break;
	   }
*/
	  lROI = l_nuev(sizeof(int));
	  RegionesInteres(lROI,2,0.2,GAR_Image);
	  GAR_RegFunction = -1;
  }else{
	  GAR_RegFunction = 0;
  }

}
/*-----------------------------------------------------------------*/
//Criterio de parada de la segmentación
bool StopCriteriorGARModel(){
  double kexttot = 0;
  float valor = GAR_MIN_DIF;
  l_empieza(GAR_CurvePoints);
  punto2D pto;
	while(l_quedan(GAR_CurvePoints)){
		l_dame(GAR_CurvePoints,&pto);

    kexttot+= GAR_kiextMtx[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)];
  }
	kexttot /= l_elementos(GAR_CurvePoints);
  /**/printf("GAR_kant= %f Kextt= %f Diferencia: %f\n",GAR_kant,kexttot,(GAR_kant - kexttot));
  if ((kexttot < 0.9)&&((GAR_kant - kexttot) < valor)&& ((GAR_kant - kexttot)>(-valor)))
  		GAR_NumRept++;
	else
		GAR_NumRept = 0;
	GAR_kant = kexttot;
  if (GAR_NumRept >= (GAR_NARROW-3)){
    return(true);
	}else{
		return(false);
  }
}
/*-----------------------------------------------------------------*/
//Función que determina el final de la narrow band, basada en un cuadrado
bool EndNarrow(){
  punto2D pto;
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    if((pto.x == 2)||(pto.y == 2)||
       (pto.x == GAR_IMAGE_HEIGHT-2)||(pto.y == GAR_IMAGE_WIDTH-2))
        return true;
  }
  return false;
}
/*-----------------------------------------------------------------*/
//Función que determina el final de la narrow band, basada en distancias
bool EndNarrowDist(){
  punto2D pto;
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
	if(GAR_dist[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)] >= GAR_NARROW)
        return true;
  }
  return false;
}

/*-----------------------------------------------------------------*/
//Libera las matrices del modelo
void DestroyGARModel(){
	//Lista de puntos
  l_dest(&GAR_CurvePoints);
  //Matriz del gradiente
  free(GAR_GradMtx);
   //Matriz del gradiente de Ki
  free(GAR_GradkiMtx);
  //Matriz de la velocidad asociada a cada punto de la imagen
  free(GAR_VelocMtx);
  //MAtriz de Ki
  free(GAR_kiextMtx);
  //Matrices del levelset
  free(GAR_fiMtx);
  //Imagen con la que vamos a tratar que puede ser una subimagen dentro de la original
  free(GAR_Image);
  free(GAR_ImageOri);
  free(GAR_dist);

  
}
/*-----------------------------------------------------------------*/
/*  FUNCIONES PARA DIBUJAR MATRICES ASOCIADAS EL LEVEL SET         */
/*-----------------------------------------------------------------*/
//Almacena el resultado actual en una imagen PGM
void DrawGARModel(){
  unsigned char *imagen_reg;
  char nombre[20];
  punto2D pto;
  sprintf(nombre,"GAR_%d",GAR_iter+100);
	strcat(nombre,".pgm");
  imagen_reg = (unsigned char*)calloc(sizeof(unsigned char) ,(GAR_IMAGE_WIDTH * GAR_IMAGE_HEIGHT));
  for(int i=0; i<GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT; i++)
      imagen_reg[i] = 0.0;
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
    l_dame(GAR_CurvePoints,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)] = 255;
  }
  guardarImagenPGM(imagen_reg, nombre);
  free(imagen_reg);
}
/*-----------------------------------------------------------------*/
//Devuelve en ls una imagen con el contorno sobreimpresionado
void PresentGRAModel(unsigned char *ls){
  punto2D pto;
  memcpy(ls,GAR_ImageOri,sizeof(unsigned char)*GAR_H_imag_ori*GAR_W_imag_ori);
  l_empieza(GAR_CurvePoints);
  while(l_quedan(GAR_CurvePoints)){
   l_dame(GAR_CurvePoints,&pto);
   pto.x = pto.x + GAR_MinNBPto.x;
   pto.y = pto.y + GAR_MinNBPto.y;
   ls[calcSingleSubscript(pto.x,pto.y,GAR_W_imag_ori)] = 255.0;
  }  
 }
/*-----------------------------------------------------------------*/
//Filtra las burbujas 
void FilterBubbles(){
  int i,j;
  punto2D pto;
  lista l=l_nuev(sizeof(punto2D));
  double *imagen_reg;
  imagen_reg = (double*)malloc(sizeof(double)*GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
  //binarizo la imagen
  GAR_BinaryImage(imagen_reg);
  for(i=0; i<GAR_IMAGE_HEIGHT-5; i++){
    for(j=0; j<GAR_IMAGE_WIDTH-5; j++){
      if((imagen_reg[calcSingleSubscript(i,j,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+1,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+2,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+3,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+4,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+1,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+3,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+4,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+2,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+4,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+1,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+3,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+4,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+1,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+2,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+3,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+4,GAR_IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+2,GAR_IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+1,GAR_IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+3,GAR_IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+2,GAR_IMAGE_WIDTH)]==1.0)){
           pto.x = i+1;
           pto.y = j+2;           
           l_meted(l,&pto);

           pto.x = i+2;
           pto.y = j+1;
           l_meted(l,&pto);
           pto.x = i+2;
           pto.y = j+3;
           l_meted(l,&pto);
           pto.x = i+3;
           pto.y = j+2;
           l_meted(l,&pto);
         }
      }
    }
    if(!l_vacia(l)){
    bool borrar=false;
    lista l_temp = l_nuev(sizeof(punto2D));
    
    l_empieza(GAR_CurvePoints);
    punto2D pto2;
    while(l_quedan(GAR_CurvePoints)){
      borrar=false;
      l_dame(GAR_CurvePoints,&pto);
      l_empieza(l);
      while(l_quedan(l)){
        l_dame(l,&pto2);
        if(pto.x == pto2.x && pto.y == pto2.y)
          borrar = true;
      }
      if(!borrar)
        l_meted(l_temp,&pto);
    }
    l_dest(&GAR_CurvePoints);
    GAR_CurvePoints=l_nuev(sizeof(punto2D));
    l_empieza(l_temp);
    while(l_quedan(l_temp)){
      l_dame(l_temp, &pto);
      l_meted(GAR_CurvePoints,&pto);
    }
	l_dest(&l_temp);
    }
	l_dest(&l);
	
    free(imagen_reg);
}
/*----------------------------------------------*/
//Determina un pixel interno al level set
bool GAR_InternalCurvePoint(double* imag,punto2D *pto){
   for(int i=1; i<GAR_IMAGE_HEIGHT-1; i++){
      for(int j=1; j<GAR_IMAGE_WIDTH-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,GAR_IMAGE_WIDTH)] < 0)&&
          (imag[calcSingleSubscript(i-1,j,GAR_IMAGE_WIDTH)] == 0)&&
          (imag[calcSingleSubscript(i,j-1,GAR_IMAGE_WIDTH)] == 0)){
            (*pto).x = i;
            (*pto).y = j;
           return true;
        }
      }
    }
   return false;
}

/*----------------------------------------------*/
//Almacen en una lista los level set con los puntos del contorno y los puntos internos
void RecoverCurveGARModel(){	
	unsigned char* pdImag;
	punto2D pto;
	TLlista lGeneral = new (TCab);
	TRegion trReg,trReg2;
	int i;
	
	pdImag=(unsigned char*)calloc(sizeof(unsigned char),GAR_IMAGE_WIDTH*GAR_IMAGE_HEIGHT);
	trReg.lContorno = l_nuev(sizeof(punto2D));
	trReg.lPuntos = l_nuev(sizeof(punto2D));
	llInicializa(lGeneral);

	l_empieza(GAR_CurvePoints);
	while(l_quedan(GAR_CurvePoints)){
		l_dame(GAR_CurvePoints,&pto);
		pdImag[calcSingleSubscript(pto.x,pto.y,GAR_IMAGE_WIDTH)] = 255;
	}
	extraerSegmentacion(pdImag,lGeneral,255);
	l_dest(&(trReg.lContorno));
	l_dest(&(trReg.lPuntos));
	delete(lGeneral);
	free(pdImag);	
}
/*----------------------------------------------*/
//Almacena los puntos del contorno en un txt
void StoreGRAModel(char* nameFile){
	punto2D pto;
	FILE *arch;

	arch = fopen(nameFile,"w");

	l_empieza(GAR_CurvePoints);
	while(l_quedan(GAR_CurvePoints)){
		l_dame(GAR_CurvePoints,&pto);
		pto.x = pto.x + GAR_MinNBPto.x;
		pto.y = pto.y + GAR_MinNBPto.y;
		fprintf(arch,"%d %d\n",pto.x,pto.y);
	}	
	fclose(arch);

}
/*----------------------------------------------*/
//Reposiciona los puntos para poder dibujarlos sobre la imagen completa
void SetUpCurveGARModel(lista l){
	punto2D pto;
	
	l_empieza(GAR_CurvePoints);
	while(l_quedan(GAR_CurvePoints)){
		l_dame(GAR_CurvePoints,&pto);
		pto.x = pto.x + GAR_MinNBPto.x;
		pto.y = pto.y + GAR_MinNBPto.y;
		l_meted(l,&pto);
	}	
	
}
/*----------------------------------------------*/
//Inidca el número de level sets
int NumLS(unsigned char* imag){
	TLlista lGen = new (TCab);
	lista lRes;
	int num=extraerSegmentacion(imag,lGen,255);
	return(num);  
}


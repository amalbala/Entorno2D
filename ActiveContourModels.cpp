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





 #include <math.h>
 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <fstream.h>
 #include "GARModel.h"
 #include "matematics.h" 
 #include "comun.h"
 #include "listalista.h"

 //ACM_fa: Factor constante en el caso de la velocidad no asociada a la curvatura
 #define ACM_fa 1.0
 #define ACM_aportation	0.1

double ACM_THRESS_KI = 0.4;
 
//Minima diferencia entre dos levelset para el criterio de parada
double ACM_MIN_DIF = 0.01;

TImagen* ACM_Image;
TImagen* ACM_ImageOri;

//Dimensiones de la Imagen

int ACM_RegFunction;

//Puntos del Level Set
lista ACM_CurvePoints;
//Level Set anterior
lista ACM_PreviousCurvePoints;
//ACM_NARROW Band
int ACM_NARROW;
//Factor de la curvatura
float ACM_EPSILON;
float ACM_SIGMA_G;
//Matriz con la imagen convolucionada con la gaussiana
TMatrizD* ACM_GradGaussMtx;
//Matriz del gradiente
TMatrizD* ACM_GradMtx;
//Matriz del gradiente de la funcion de parada
TMatrizD* ACM_GradkiMtx;
//Matriz de la velocidad asociada a cada punto de la imagen
TMatrizD* ACM_VelocMtx;
//Matriz del criterio de parada ACM_kiMtx extendido
TMatrizD* ACM_kiextMtx;
TMatrizD* ACM_kiMtx;

//Matrices del levelset
//Matriz de evolucion del level set
TMatrizD* ACM_fiMtx;
TMatrizD* ACM_fiOldMtx;
//Seleccion con circulos
bool ACM_circul = false;

//Numero de iteraciones
int ACM_iter=1;
//ACM_modo de evolucion: 0 = expansion, 1=compresion
int ACM_modo;
//Criterio de Parada
float ACM_kant;
//Numero de repeticiones del mismo valor del level set
int ACM_NumRept=0;
//Para localizar un punto interno del level set
punto2D ACM_InternalPto;
punto2D ACM_MinNBPto;
punto2D ACM_MaxNBPto;
punto2D ACM_OldMinNBPto;

//Etiqueta Mayoritaria
int ACM_BestLabel;
TMatrizI* ACM_LabelMtx;
//Media de intensidades de una region
double ACM_GAR_fIntMin;
double ACM_GAR_fDesviation;
int ACM_GAR_fIntMax;
//Coeficiente de peso entre el modulo de region y el modulo de contorno

double ACM_BETA = 0.5;
double ACM_GAR_ValorRegion;


TMatrizI *ACM_dist;
/*--------------------------------------------------*/
/*            ACTIVE CONTOUR MODELS					*/
/*													*/
/*	Utility functions. Common to all the models		*/
/*--------------------------------------------------*/

//New narrow band function
void ACM_NewNB(lista l,int nb_size,punto2D *pto_min,punto2D *pto_max){
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

  (*pto_min).x = (x1-nb_size);
  if((*pto_min).x < 0) (*pto_min).x = 0;
  (*pto_min).y = (y1-nb_size);
  if((*pto_min).y < 0) (*pto_min).y = 0;
  (*pto_max).x = (x2+nb_size);
  if((*pto_max).x >= ACM_ImageOri->height) (*pto_max).x = ACM_ImageOri->height-1;
  (*pto_max).y = (y2+nb_size);
  if((*pto_max).y >= ACM_ImageOri->width) (*pto_max).y = ACM_ImageOri->width-1;

}


/*-----------------------------------------------*/
//Genera una imagen binaria en la que el fondo el 0.0 y los puntos del level set 1.0
void ACM_BinaryImage(TMatrizI *imag_bin, lista lPuntos){
  punto2D pto;
  int i;
  /*for(i=0; i<ACM_Image->width*ACM_Image->height;i++)
    imagen_reg[i] = 0.0;
	*/
  memset(imag_bin->matrix,0,sizeof(int)*imag_bin->height * imag_bin->width);
  l_empieza(lPuntos);
  while(l_quedan(lPuntos)){
    l_dame(lPuntos,&pto);
    imag_bin->matrix[calcSingleSubscript(pto.x,pto.y,imag_bin->width)] = 1;
  }
}
/*-----------------------------------------------*/
//Genera una imagen binaria en la que el fondo el 0.0 y los puntos del level set 1.0
void ACM_BinaryImage(TMatrizD *imag_bin, lista lPuntos){
  punto2D pto;
  int i;
  /*for(i=0; i<ACM_Image->width*ACM_Image->height;i++)
    imagen_reg[i] = 0.0;
	*/
  memset(imag_bin->matrix,0,sizeof(double)*imag_bin->height * imag_bin->width);
  l_empieza(lPuntos);
  while(l_quedan(lPuntos)){
    l_dame(lPuntos,&pto);
    imag_bin->matrix[calcSingleSubscript(pto.x,pto.y,imag_bin->width)] = 1.0;
  }
}
/*-----------------------------------------------*/
//Crecimiento de regiones para determinar los signos de las zonas del level set
void ACM_SignDetermination(){
  int valor,i;
  TMatrizI *imag2;
  int *etiquetas;
  
  imag2 = new(TMatrizI);
  imag2->matrix = (int*)calloc(sizeof(int),ACM_Image->width*ACM_Image->height);
  etiquetas = (int*)calloc(sizeof(int),ACM_Image->width*ACM_Image->height);

  ACM_BinaryImage(imag2,ACM_CurvePoints);
  etiq_region(imag2->matrix,etiquetas);
  free(imag2->matrix);
  delete(imag2);
  valor = etiquetas[calcSingleSubscript(0, 0,ACM_Image->width)];
  if(ACM_modo == 1){
      for(i=0; i<ACM_fiMtx->width*ACM_fiMtx->height;i++){
    //Para los puntos localizados cambio el signo.
      if(etiquetas[i] == valor){
        ACM_fiMtx->matrix[i] *= -1;
      }
    }
    }else{
      for(i=0; i<ACM_fiMtx->width*ACM_fiMtx->height;i++){
      if(etiquetas[i] != valor){
        ACM_fiMtx->matrix[i] *= -1;
      }
    }
  }
  free(etiquetas);
}
  
/*-----------------------------------------------------------------*/
//Determina si un punto pertenece o no al level set
bool ACM_CurvePoint(int i,int j){
  //Si el valor es menos que cero
  if((ACM_fiMtx->matrix[calcSingleSubscript(i,j,ACM_fiMtx->width)] < 0)&&
  //y alguno de sus 4 vecinos es positivo, entonces es un punto del level set.
        ((ACM_fiMtx->matrix[calcSingleSubscript(i-1,j,ACM_fiMtx->width)] >= 0)||  
         (ACM_fiMtx->matrix[calcSingleSubscript(i+1,j,ACM_fiMtx->width)] >= 0)||
         (ACM_fiMtx->matrix[calcSingleSubscript(i,j-1,ACM_fiMtx->width)] >= 0)||
         (ACM_fiMtx->matrix[calcSingleSubscript(i,j+1,ACM_fiMtx->width)] >= 0)))
          return true;
     return false;
}
/*-----------------------------------------------------------------*/
//Determina los puntos del nuevo level set
void ACM_NewCurve(){
  int i,j;
  punto2D pto;
  //Reinicializamos la lista de puntos del LevelSet
  l_dest(&ACM_CurvePoints);
  ACM_CurvePoints = l_nuev(sizeof(punto2D));
  //para todos los puntos de la imagen determino aquellos que pertenecen al level set
  for(i=1; i<ACM_Image->height-1; i++){
    for(j =1; j<ACM_Image->width-1; j++){
      //Si pertenece lo inserto en la nueva lista
       if(ACM_CurvePoint(i,j)){
           pto.x = i;
           pto.y = j;
           l_meted(ACM_CurvePoints,&pto);

       }
     }
  }

}
/*-----------------------------------------------------------------*/
//Calculo el mapa de distancias con signo aplicando un proceso incremental.
void ACM_DistanceMap(){
  int i,j,a,b;   
  free(ACM_dist);
  //Calculo la Chamfer Distance
  chamfer_distance4v(ACM_fiMtx->matrix,ACM_fiMtx->matrix);
  //Inicializo la matriz que luego controlará la narrow band proporcional
  ACM_dist->height = ACM_fiMtx->height;
  ACM_dist->width = ACM_fiMtx->width;
  ACM_dist = (TMatrizI*)malloc(sizeof(TMatrizI)*ACM_dist->height * ACM_dist->width);
  ACM_dist->matrix = (int*)malloc(sizeof(int) * ACM_dist->width*ACM_dist->height);
  for(int n=0; n<ACM_dist->width*ACM_dist->height; n++)
	  ACM_dist->matrix[n] = ACM_fiMtx->matrix[n];
  //Calculo el signo de las distancias a partir del valor anterior
  if (ACM_modo == 0){
    for(i=0; i<ACM_fiOldMtx->height;i++){
        for(j=0; j<ACM_fiOldMtx->width;j++){
          if(ACM_fiOldMtx->matrix[calcSingleSubscript(i,j,ACM_fiOldMtx->width)]< 0.0){
            a = i+ACM_OldMinNBPto.x-ACM_MinNBPto.x;
            b = j+ACM_OldMinNBPto.y-ACM_MinNBPto.y;
            if(a>=0 && b>=0 && a<ACM_Image->height && b<ACM_Image->width){
              ACM_fiMtx->matrix[calcSingleSubscript(a,b,ACM_Image->width)]=
                -ACM_fiMtx->matrix[calcSingleSubscript(a,b,ACM_Image->width)];
            }
          }
        }
    }
    for(i=0;i<ACM_Image->height;i++){
      if(ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)]<0)
        ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)];
      if(ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)]<0)
        ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)];

    }
    for(j=0;j<ACM_Image->width;j++){
      if(ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)]<0)
        ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)];
      if(ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)]<0)
        ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)];
    }
 
  }else{


    for(i=0; i<ACM_fiOldMtx->height;i++){
      for(j=0; j<ACM_fiOldMtx->width;j++){         
        if(ACM_fiOldMtx->matrix[calcSingleSubscript(i,j,ACM_fiOldMtx->width)]< 0.0){
         a = i+ACM_OldMinNBPto.x-ACM_MinNBPto.x;
         b = j+ACM_OldMinNBPto.y-ACM_MinNBPto.y;
         if(a>=0 && b>=0 && a<ACM_Image->height && b<ACM_Image->width){
           ACM_fiMtx->matrix[calcSingleSubscript(a,b,ACM_Image->width)]=
            -ACM_fiMtx->matrix[calcSingleSubscript(a,b,ACM_Image->width)];
         }
        }
      }
    }
    for(i=0;i<ACM_Image->height;i++){
      if(ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)]>0)
        ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(i,0,ACM_Image->width)];
      if(ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)]>0)
        ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(i,ACM_Image->width-1,ACM_Image->width)];

    }
    for(j=0;j<ACM_Image->width;j++){
      if(ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)]>0)
        ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(0,j,ACM_Image->width)];
      if(ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)]>0)
        ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)]=
          -ACM_fiMtx->matrix[calcSingleSubscript(ACM_Image->height-1,j,ACM_Image->width)];
    } 
  }
  delete(ACM_fiOldMtx);  
}
/*-----------------------------------------------*/
//calculo las distancias basadas en funciones recursivas
void ACM_DistanceMapLabel(){
  //Calculo la Chamfer Distance
  chamfer_distance4v(ACM_fiMtx->matrix,ACM_fiMtx->matrix);
  //Inicializo la matriz que luego controlará la narrow band proporcional
  ACM_dist=(TMatrizI*)malloc(sizeof(TMatrizI));
  ACM_dist->height = ACM_fiMtx->height;
  ACM_dist->width = ACM_fiMtx->width;
  ACM_dist->matrix = (int*)malloc(sizeof(int) * ACM_Image->height*ACM_Image->width);
  for(int n=0; n<ACM_Image->height*ACM_Image->width; n++)
	  ACM_dist->matrix[n] = ACM_fiMtx->matrix[n];
  //Aplicamos el algoritmo de relleno para determinar los signos
  ACM_SignDetermination();

}
/*-----------------------------------------------*/
//Caclculo la velocidad dependiente de la curvatura o ACM_VelocMtx
void ACM_VelocityMapCurvature(){
	float valor;
	double *curv;
	curv=(double*)calloc(sizeof(double),(ACM_Image->height*ACM_Image->width));  
	curvatura(ACM_fiMtx->matrix,curv,grad_diffin);
	for(int i=0;i<ACM_Image->width*ACM_Image->height;i++){
		valor = curv[i];
		ACM_VelocMtx->matrix[i] = (-ACM_EPSILON) * valor;
	}
	free(curv);
}


/*-----------------------------------*/
//Calculo el factor de parada ACM_kiMtx
void ACM_StopCriterion(){
  int i;
  float a;
  
  ACM_GradGaussMtx=CreateMatrizD(ACM_kiextMtx);
  grad_conv_gauss(ACM_Image->matrix,ACM_SIGMA_G,&a,&a,ACM_GradGaussMtx->matrix);
  for(i=0; i<ACM_Image->width*ACM_Image->height; i++){
    ACM_kiextMtx->matrix[i]= 1.0/(1.0 + ACM_GradGaussMtx->matrix[i]);
  }
  DestroyMatrix(&ACM_GradGaussMtx);
}
/*-----------------------------------*/
//Calculo el factor de parada extendido ACM_kiMtx
void ACM_StopCriterionExtension(){
  int i,j;
  double *level;
  punto2D *dist_4sed;
  punto2D pto;
  int tama =ACM_Image->height*ACM_Image->width;
  float a;
  double v_grad;
  level = (double*)calloc(sizeof(double),tama);
  
  ACM_GradGaussMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  grad_conv_gauss(ACM_Image->matrix,ACM_SIGMA_G,&a,&a,ACM_GradGaussMtx->matrix);
  for(i=0;i<tama;i++)
    ACM_kiextMtx->matrix[i] = 0.0;
  l_empieza(ACM_CurvePoints);
  //La velocidad para los puntos del Level set
   while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
    //El valor es igual a:
    /*
                  1
    ACM_kiMtx =  ------------------
          (1 + |G * I(x,y)|)
    */
    v_grad = 1.0/(1.0 + ACM_GradGaussMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]);
    if(v_grad<ACM_THRESS_KI)
      ACM_kiextMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]= 0.0;
    else
      ACM_kiextMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]= v_grad;
      
    level[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]= 1.0;
  }
   //Para hallar los puntos de level set mas cercano a cada punto de la imagen
  //aplico un mapa de distancias 4SSED
   free(ACM_GradGaussMtx);
   dist_4sed = (punto2D*)calloc(sizeof(punto2D),tama);
   vecin4ssed_distance(level,dist_4sed);
   //Los puntos que no pertenecen al level set tienen una valor igual a la velocidad del
  //punto de la level set mas cercana a ese punto.
  int pos;
  for(i=1;i<ACM_Image->height-1;i++){
    for(j=1;j<ACM_Image->width-1;j++){
        pos= calcSingleSubscript(i,j,ACM_Image->width);
        ACM_kiextMtx[pos] =
          ACM_kiextMtx[calcSingleSubscript(i - dist_4sed[pos].x,j - dist_4sed[pos].y,ACM_Image->width)];
    }
  }

  free(level);
  free(dist_4sed);
}
/*-----------------------------------------------*/
void ReInitACMModel(){
  punto2D pto;
  //Almacenamos la matriz actual de ACM_fiMtx
  
  //free(ACM_GradGaussMtx);
  //Matriz del gradiente
  free(ACM_GradMtx);
  free(ACM_GradkiMtx);
  //Matriz de la velocidad asociada a cada punto de la imagen
  free(ACM_VelocMtx);
  free(ACM_kiextMtx);
  free(ACM_Image);
  
  //Guarda el level set antiguo 
  lista ACM_CurvePoints_old;
  ACM_CurvePoints_old = l_nuev(sizeof(punto2D));
  l_empieza(ACM_CurvePoints);
  while(l_quedan(ACM_CurvePoints)){
     l_dame(ACM_CurvePoints, &pto);
     pto.x += ACM_MinNBPto.x;

     pto.y += ACM_MinNBPto.y;
     l_meted(ACM_CurvePoints_old,&pto);
   }
  l_dest(&ACM_CurvePoints);
  ACM_CurvePoints = l_nuev(sizeof(punto2D));
  //Guarda los valores antiguos de Fi y de la Narrow Band
  ACM_OldMinNBPto=ACM_MinNBPto;  
  ACM_fiOldMtx = CreateMatrizD(ACM_fiMtx);
  DestroyMatrix(&ACM_fiMtx);
  //Determina una nueva narrow band
  ACM_NewNB(ACM_CurvePoints_old,ACM_NARROW,&ACM_MinNBPto,&ACM_MaxNBPto);
  //Actualiza los puntos del level set a la nueva narrow band
  l_empieza(ACM_CurvePoints_old);
  while(l_quedan(ACM_CurvePoints_old)){
     l_dame(ACM_CurvePoints_old, &pto);
     pto.x -= ACM_MinNBPto.x;
     pto.y -= ACM_MinNBPto.y;
     l_meted(ACM_CurvePoints,&pto);
   }
  l_dest(&ACM_CurvePoints_old);

  ACM_Image = CreateImagen(ACM_MaxNBPto.x -ACM_MinNBPto.x,ACM_MaxNBPto.y -ACM_MinNBPto.y);
  int i,j,k=0;
  for(i=ACM_MinNBPto.x; i<ACM_MaxNBPto.x; i++){
    for(j=ACM_MinNBPto.y; j<ACM_MaxNBPto.y; j++){
     ACM_Image[k]=ACM_ImageOri[calcSingleSubscript(i,j,ACM_ImageOri->width)];
     k++;
    }
  }
  
  Set_X_Y(ACM_Image->width,ACM_Image->height);
  //Matriz del gradiente
  ACM_GradMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz del gradiente del criterio de parada
  ACM_GradkiMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz de la velocidad asociada a cada punto de la imagen
  ACM_VelocMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz del criterio de parada extendido
  ACM_kiextMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz del levelset
  ACM_fiMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
 
  l_empieza(ACM_CurvePoints);
   while(l_quedan(ACM_CurvePoints)){
     l_dame(ACM_CurvePoints, &pto);
     ACM_fiMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]= 1.0;
   }
  //Calculamos la matriz de distancias con signo aplicando recursion
  ACM_DistanceMap();  
  //Calculamos la matriz de velocidad
  ACM_VelocityMapCurvature();
  //Calculamos la matriz asociada al criterio de parada
  ACM_StopCriterionExtension();
}
/*-----------------------------------------------------------------*/
//Filtra puntos repetidos en la lista del level set
void ACM_CurveFilter(){
  punto2D pto,ptoant;
  l_empieza(ACM_CurvePoints);
  lista filtrols;
  filtrols = l_nuev(sizeof(punto2D));
  while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
    l_meted(filtrols,&pto);
  }

  l_dest(&ACM_CurvePoints);
  ACM_CurvePoints = l_nuev(sizeof(punto2D));
  l_empieza(filtrols);
  ptoant.x = ptoant.y = 0;
  while(l_quedan(filtrols)){
    l_dame(filtrols,&pto);
    if((pto.x>0)&&(pto.y>0)&&((pto.x != ptoant.x)||(pto.y != ptoant.y)))
      l_meted(ACM_CurvePoints,&pto);
    ptoant.x = pto.x;
    ptoant.y = pto.y;
  }
  l_dest(&filtrols);
}
/*-----------------------------------------------------------------*/
//Guarda el level set anterior
void ACM_StorePreviousCurve(){
  punto2D pto;
  l_empieza(ACM_CurvePoints);
  l_dest(&ACM_PreviousCurvePoints);
  ACM_PreviousCurvePoints= l_nuev(sizeof(punto2D));
  while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
    l_meted(ACM_PreviousCurvePoints,&pto);
  }
}
/*-----------------------------------------------------------------*/
//Calculo del Termino de la Region 
double ACM_RegionTerm(int i){
	double dDifMin, dDifMax;

	//Establece la distancia de una intensidad de la imagen al intervalo de interés, estudiando
	//la distancia a los extremo del intervalo.
	dDifMin = ACM_Image->matrix[i] - ACM_GAR_fIntMin;
	dDifMax = ACM_Image->matrix[i] - ACM_GAR_fIntMax;
	//El valor será la distancia euclidea a ambos extremos
	return (sqrt((dDifMin * dDifMin)+(dDifMax * dDifMax)));	

}
/*-----------------------------------------------------------------*/
//Calculo del Teremino del Contorno
double ACM_BoundaryTerm(int i){
	//Es la funcion tipica de las geodesicas
	//(g(1.0 + K) * N) - (g · N)
	return ((ACM_kiextMtx->matrix[i]*(ACM_fa + ACM_VelocMtx->matrix[i])*ACM_GradMtx->matrix[i]) - 
		(ACM_GradkiMtx->matrix[i]*ACM_GradMtx->matrix[i]));		
}
/*-----------------------------------------------------------------*/
//Calcula la media de intensidad de la region contenida por el level set
double ACM_GAR_RegionIntensityAverage(){
	double dTotal=0;
	int iNum=0;
	for(int i=0;i<ACM_Image->width*ACM_Image->height;i++){
		if(ACM_fiMtx->matrix[i]<0){
			iNum++;
			dTotal += ACM_Image->matrix[i];
		}
	}
	return(dTotal/iNum);
}
/*-----------------------------------------------------------------*/
//Calcula el máximo del histograma de la region contenida por el level set
int ACM_GAR_RegionIntensityMax(){
	int vHistRegion[255];
	int i;
	
	for(i=0; i<255; i++) vHistRegion[i]= 0;

	for(i=0;i<ACM_Image->width*ACM_Image->height;i++){
		if(ACM_fiMtx->matrix[i]<0){
			vHistRegion[ACM_Image->matrix[i]]++;	
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
double ACM_GAR_RegionStandardDesviation(double iAverage){	
	int iNum=0;
	double dSum = 0;
	if(iAverage==0){
		iAverage = ACM_GAR_RegionIntensityAverage();		
	}
	for(int i=0;i<ACM_Image->width*ACM_Image->height;i++){
			if(ACM_fiMtx->matrix[i]<0){
				iNum++;
				dSum += iAverage + ACM_Image->matrix[i];
			}
	}
	return(sqrt(dSum/iNum));
	
}
/*-----------------------------------------------------------------*/
//Siguiente iteracion del modelo
void NextIterACMModel(float inct){
  float valor,nuevovalor,media;
  int i;
  
  //Mapa de Velocidades seguna la curvatura
  ACM_VelocityMapCurvature();
  //Mapa del Criterio de Parada
  ACM_StopCriterionExtension();
  //Calculo del gradiente de Fi
  grad_diffin(ACM_fiMtx->matrix,ACM_GradMtx->matrix);
   //Calculo el gradiente de Criterio de Parada
  grad_diffin(ACM_kiextMtx->matrix,ACM_GradkiMtx->matrix);
  switch(ACM_RegFunction){
	  //Segun el valor, se tomará como centro de intervalo de interes
	  //la media de la región o el máximo del histograma de la region.
  case 0:
	ACM_GAR_ValorRegion = ACM_GAR_RegionIntensityAverage();
	break;
  case 1:
	  ACM_GAR_ValorRegion = ACM_GAR_RegionIntensityMax();
	  break;
  }

  //Calculo la desviacion estandar que determina la amplitud del intervalo
  ACM_GAR_fDesviation = ACM_GAR_RegionStandardDesviation(ACM_GAR_ValorRegion);
  ACM_GAR_fIntMin = ACM_GAR_ValorRegion - ACM_GAR_fDesviation;
  ACM_GAR_fIntMax = ACM_GAR_ValorRegion + ACM_GAR_fDesviation;
  
  //Matrices que contendran el valor de los terminos de contorno y region
  double *BoundaryMtx,*RegionMtx;
  BoundaryMtx = (double*)malloc(sizeof(double)*ACM_Image->width*ACM_Image->height);
  RegionMtx = (double*)malloc(sizeof(double)*ACM_Image->width*ACM_Image->height);
 
  //Valores para el reescalado de las matrices
  double MaxBoundary, MinBoundary;
  double MaxRegion, MinRegion;
 
  MaxBoundary = MinBoundary = ACM_BoundaryTerm(0);
  MaxRegion = MinRegion = ACM_RegionTerm(0);
  
  
  //Vamos calculando los valores de la matrices y acutalizando los maximos y minimos 
  //para poder reescalar
  for(i=0;i<ACM_Image->width*ACM_Image->height;i++){	  
	  BoundaryMtx[i] = ACM_BoundaryTerm(i);
	  if(BoundaryMtx[i] > MaxBoundary) MaxBoundary = BoundaryMtx[i];
	  if(BoundaryMtx[i] < MinBoundary) MinBoundary = BoundaryMtx[i];
	  
	  RegionMtx[i] = ACM_RegionTerm(i);
	  if(RegionMtx[i] > MaxRegion) MaxRegion = RegionMtx[i];
	  if(RegionMtx[i] < MinRegion) MinRegion = RegionMtx[i];
  }
  
  //Normalización
  double BoundaryInterv = MaxBoundary - MinBoundary;
  double RegionInterv = MaxRegion - MinRegion;
  for(i=0;i<ACM_Image->width*ACM_Image->height;i++){	  
	  BoundaryMtx[i] = (BoundaryMtx[i] - MinBoundary)/BoundaryInterv;
	  RegionMtx[i] = (RegionMtx[i] - MinRegion)/RegionInterv;
  }
 
  //Ecuación general de la ACM donde se calcula el nuevo valor de Fi teniendo en 
  //cuenta ambos terminos
  // Fi_nuevo = Fi_anterior - t * ((Aporte * TerminoContorno) - ((1-Aporte) * TerminoRegion))
  for(i=0;i<ACM_Image->width*ACM_Image->height;i++){	  
		  ACM_fiMtx->matrix[i] = ACM_fiMtx->matrix[i] - (inct *((ACM_BETA*BoundaryMtx[i]) - ((1.0 - ACM_BETA)*RegionMtx[i])));
	}
	
  free(BoundaryMtx);
  free(RegionMtx);
  ACM_NewCurve();
  ACM_iter++;
  
}
/*-----------------------------------------------------------------*/
//Inicialización del modelo
void InitACMModel(unsigned char* imagen_reg,int mod,float sigma,int NARROW,float EPSILON,
				  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,int bInterpol, 
				  float weight,int RegionFunc, float ValorReg, int OnlyInit){
  punto2D pto;
  int i,j,k;
  ACM_modo = mod;
  ACM_NARROW = NARROW;
  ACM_EPSILON = EPSILON;
  ACM_MIN_DIF = umbral;
  ACM_SIGMA_G = sigma;
  ACM_THRESS_KI = th_ki;
  ACM_ImageOri->width = width_imag;
  ACM_ImageOri->height = height_imag;
  ACM_iter = 1;
  ACM_NumRept=0;
  ACM_BETA = weight;
  
  //Determina la nueva narrow band  
  ACM_NewNB(selec,ACM_NARROW,&ACM_MinNBPto,&ACM_MaxNBPto);
  ACM_Image->height = ACM_MaxNBPto.x -ACM_MinNBPto.x;
  ACM_Image->width = ACM_MaxNBPto.y -ACM_MinNBPto.y;
  //Inicializa matrices
  ACM_Image = CreateImagen(ACM_Image->height,ACM_Image->width);
  ACM_ImageOri = CreateImagen(ACM_ImageOri->width,ACM_ImageOri->height);
  k=0;
  memcpy(ACM_ImageOri,imagen_reg,sizeof(unsigned char)*ACM_ImageOri->width*ACM_ImageOri->height);
  for(i=ACM_MinNBPto.x; i<ACM_MaxNBPto.x; i++){
    for(j=ACM_MinNBPto.y; j<ACM_MaxNBPto.y; j++){
     ACM_Image->matrix[k]=imagen_reg[calcSingleSubscript(i,j,ACM_ImageOri->width)];
     k++;
    }
  }
  Set_X_Y(ACM_Image->width,ACM_Image->height);
  
  //Matriz del gradiente
  ACM_GradMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //MAtriz del gradiente del criterio de parada Ki
  ACM_GradkiMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz de la velocidad asociada a cada punto de la imagen
  ACM_VelocMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matriz del criterio de parada
  ACM_kiextMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //Matrices del levelset
  ACM_fiMtx=CreateMatrizD(ACM_Image->height,ACM_Image->width);

  ACM_CurvePoints = l_nuev(sizeof(punto2D));
  lista pto_temp= l_nuev(sizeof(punto2D));
  //Interpolación de la entrada si no son circulos
  if(!circ){
  //Como puede que no sean puntos 4 conexos interpolamos
	  if(bInterpol){
		  cierre_lineal(selec,l_elementos(selec),pto_temp);
		  l_empieza(pto_temp);
		  while(l_quedan(pto_temp)){
			  l_dame(pto_temp,&pto);
			  pto.x -= ACM_MinNBPto.x;
			  pto.y -= ACM_MinNBPto.y;
			  l_meted(ACM_CurvePoints,&pto);
		  }
		  
	  }else{
		  l_empieza(selec);
		  while(l_quedan(selec)){
			  l_dame(selec,&pto);
			  pto.x -= ACM_MinNBPto.x;
			  pto.y -= ACM_MinNBPto.y;
			  l_meted(ACM_CurvePoints,&pto);
		  }
	  }
   //Filtramos si hay puntos repetidos
	ACM_CurveFilter();
   }else{
     ACM_circul = true;
     l_empieza(selec);
	 while(l_quedan(selec)){
		 l_dame(selec, &pto);
		 pto.x -= ACM_MinNBPto.x;
		 pto.y -= ACM_MinNBPto.y;
		 l_meted(ACM_CurvePoints,&pto);
	 }
   }

   l_dest(&pto_temp);
   l_empieza(ACM_CurvePoints);
   //Almacena los puntos del level set 
   while(l_quedan(ACM_CurvePoints)){
     l_dame(ACM_CurvePoints, &pto);
     ACM_fiMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)]= 1.0;
   }
  //Calculamos la matriz de distancias con signo aplicando recursion
  ACM_DistanceMapLabel();
  //Calculamos la matriz de velocidad
  ACM_VelocityMapCurvature();
  //Calculamos la matriz asociada al criterio de parada
  ACM_StopCriterionExtension();
  //El tipo de función de region seleccionada
  ACM_RegFunction = RegionFunc;
  if (ValorReg != 0){
	  ACM_RegFunction = -1;	  
  }
  //Si solo queremos que se calcule en la primera iteracion
  if(OnlyInit){
	   switch(ACM_RegFunction){
		   case 0:
			   ACM_GAR_ValorRegion = ACM_GAR_RegionIntensityAverage();
		   break;
		   case 1:
			   ACM_GAR_ValorRegion = ACM_GAR_RegionIntensityMax();
		   break;
	   }
	   ACM_RegFunction = -1;
  }

}
/*-----------------------------------------------------------------*/
//Criterio de parada de la segmentación
bool StopCriteriorACMModel(){
  double kexttot = 0;
  float valor = ACM_MIN_DIF;
  l_empieza(ACM_CurvePoints);
  punto2D pto;
	while(l_quedan(ACM_CurvePoints)){
		l_dame(ACM_CurvePoints,&pto);

    kexttot+= ACM_kiextMtx->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)];
  }
	kexttot /= l_elementos(ACM_CurvePoints);
  /**/printf("ACM_kant= %f Kextt= %f Diferencia: %f\n",ACM_kant,kexttot,(ACM_kant - kexttot));
  if ((kexttot < 0.9)&&((ACM_kant - kexttot) < valor)&& ((ACM_kant - kexttot)>(-valor)))
  		ACM_NumRept++;
	else
		ACM_NumRept = 0;
	ACM_kant = kexttot;
  if (ACM_NumRept >= (ACM_NARROW-3)){
    return(true);
	}else{
		return(false);
  }
}
/*-----------------------------------------------------------------*/
//Función que determina el final de la narrow band, basada en un cuadrado
bool ACM_EndNarrow(){
  punto2D pto;
  l_empieza(ACM_CurvePoints);
  while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
    if((pto.x == 2)||(pto.y == 2)||
       (pto.x == ACM_Image->height-2)||(pto.y == ACM_Image->width-2))
        return true;
  }
  return false;
}
/*-----------------------------------------------------------------*/
//Función que determina el final de la narrow band, basada en distancias
bool ACM_EndNarrowDist(){
  punto2D pto;
  l_empieza(ACM_CurvePoints);
  while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
	if(ACM_dist->matrix[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)] >= ACM_NARROW)
        return true;
  }
  return false;
}

/*-----------------------------------------------------------------*/
//Libera las matrices del modelo
void ACM_DestroyModel(){
	//Lista de puntos
  l_dest(&ACM_CurvePoints);
  //Matriz del gradiente
  free(ACM_GradMtx);
   //Matriz del gradiente de Ki
  free(ACM_GradkiMtx);
  //Matriz de la velocidad asociada a cada punto de la imagen
  free(ACM_VelocMtx);
  //MAtriz de Ki
  free(ACM_kiextMtx);
  //Matrices del levelset
  free(ACM_fiMtx);
  //Imagen con la que vamos a tratar que puede ser una subimagen dentro de la original
  free(ACM_Image);
  free(ACM_ImageOri);
  free(ACM_dist);

  
}
/*-----------------------------------------------------------------*/
/*  FUNCIONES PARA DIBUJAR MATRICES ASOCIADAS EL LEVEL SET         */
/*-----------------------------------------------------------------*/
//Almacena el resultado actual en una imagen PGM
void ACM_DrawModel(){
  unsigned char *imagen_reg;
  char nombre[20];
  punto2D pto;
  sprintf(nombre,"ACM_%d",ACM_iter+100);
	strcat(nombre,".pgm");
  imagen_reg = (unsigned char*)calloc(sizeof(unsigned char) ,(ACM_Image->width * ACM_Image->height));
  for(int i=0; i<ACM_Image->width*ACM_Image->height; i++)
      imagen_reg[i] = 0.0;
  l_empieza(ACM_CurvePoints);
  while(l_quedan(ACM_CurvePoints)){
    l_dame(ACM_CurvePoints,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)] = 255;
  }
  guardarImagenPGM(imagen_reg, nombre);
  free(imagen_reg);
}
/*-----------------------------------------------------------------*/
//Devuelve en ls una imagen con el contorno sobreimpresionado
void ACM_PresentModel(unsigned char *ls){
  punto2D pto;
  memcpy(ls,ACM_ImageOri,sizeof(unsigned char)*ACM_ImageOri->height*ACM_ImageOri->width);
  l_empieza(ACM_CurvePoints);
  while(l_quedan(ACM_CurvePoints)){
   l_dame(ACM_CurvePoints,&pto);
   pto.x = pto.x + ACM_MinNBPto.x;
   pto.y = pto.y + ACM_MinNBPto.y;
   ls[calcSingleSubscript(pto.x,pto.y,ACM_ImageOri->width)] = 255.0;
  }  
 }
/*-----------------------------------------------------------------*/
//Filtra las burbujas 
void ACM_FilterBubbles(){
  int i,j;
  punto2D pto;
  lista l=l_nuev(sizeof(punto2D));
  TMatrizD *imagen_reg = CreateMatrizD(ACM_Image->height,ACM_Image->width);
  //binarizo la imagen
  ACM_BinaryImage(imagen_reg,ACM_CurvePoints);
  for(i=0; i<ACM_Image->height-5; i++){
    for(j=0; j<ACM_Image->width-5; j++){
      if((imagen_reg->matrix[calcSingleSubscript(i,j,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i,j+1,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i,j+2,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i,j+3,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i,j+4,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+1,j,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+1,j+1,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+1,j+3,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+1,j+4,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+2,j,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+2,j+2,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+2,j+4,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+3,j,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+3,j+1,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+3,j+3,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+3,j+4,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+4,j,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+4,j+1,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+4,j+2,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+4,j+3,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+4,j+4,ACM_Image->width)]==0.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+1,j+2,ACM_Image->width)]==1.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+2,j+1,ACM_Image->width)]==1.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+2,j+3,ACM_Image->width)]==1.0)&&(
          imagen_reg->matrix[calcSingleSubscript(i+3,j+2,ACM_Image->width)]==1.0)){
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
    
    l_empieza(ACM_CurvePoints);
    punto2D pto2;
    while(l_quedan(ACM_CurvePoints)){
      borrar=false;
      l_dame(ACM_CurvePoints,&pto);
      l_empieza(l);
      while(l_quedan(l)){
        l_dame(l,&pto2);
        if(pto.x == pto2.x && pto.y == pto2.y)
          borrar = true;
      }
      if(!borrar)
        l_meted(l_temp,&pto);
    }
    l_dest(&ACM_CurvePoints);
    ACM_CurvePoints=l_nuev(sizeof(punto2D));
    l_empieza(l_temp);
    while(l_quedan(l_temp)){
      l_dame(l_temp, &pto);
      l_meted(ACM_CurvePoints,&pto);
    }
	l_dest(&l_temp);
    }
	l_dest(&l);
	
    DestroyMatrix(&imagen_reg);
}
/*----------------------------------------------*/
//Determina un pixel interno al level set
bool ACM_InternalCurvePoint(double* imag,punto2D *pto){
   for(int i=1; i<ACM_Image->height-1; i++){
      for(int j=1; j<ACM_Image->width-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,ACM_Image->width)] < 0)&&
          (imag[calcSingleSubscript(i-1,j,ACM_Image->width)] == 0)&&
          (imag[calcSingleSubscript(i,j-1,ACM_Image->width)] == 0)){
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
void ACM_RecoverCurveModel(){	
	unsigned char* pdImag;
	punto2D pto;
	TLlista lGeneral = new (TCab);
	TRegion trReg,trReg2;
	int i;
	
	pdImag=(unsigned char*)calloc(sizeof(unsigned char),ACM_Image->width*ACM_Image->height);
	trReg.lContorno = l_nuev(sizeof(punto2D));
	trReg.lPuntos = l_nuev(sizeof(punto2D));
	llInicializa(lGeneral);

	l_empieza(ACM_CurvePoints);
	while(l_quedan(ACM_CurvePoints)){
		l_dame(ACM_CurvePoints,&pto);
		pdImag[calcSingleSubscript(pto.x,pto.y,ACM_Image->width)] = 255;
	}
	extraerSegmentacion(pdImag,lGeneral,255);
	l_dest(&(trReg.lContorno));
	l_dest(&(trReg.lPuntos));
	delete(lGeneral);
	free(pdImag);	
}
/*----------------------------------------------*/
//Almacena los puntos del contorno en un txt
void ACM_StoreModel(char* nameFile){
	punto2D pto;
	FILE *arch;

	arch = fopen(nameFile,"w");

	l_empieza(ACM_CurvePoints);
	while(l_quedan(ACM_CurvePoints)){
		l_dame(ACM_CurvePoints,&pto);
		pto.x = pto.x + ACM_MinNBPto.x;
		pto.y = pto.y + ACM_MinNBPto.y;
		fprintf(arch,"%d %d\n",pto.x,pto.y);
	}	
	fclose(arch);

}
/*----------------------------------------------*/
//Reposiciona los puntos para poder dibujarlos sobre la imagen completa
void ACM_SetUpCurveModel(lista l){
	punto2D pto;
	
	l_empieza(ACM_CurvePoints);
	while(l_quedan(ACM_CurvePoints)){
		l_dame(ACM_CurvePoints,&pto);
		pto.x = pto.x + ACM_MinNBPto.x;
		pto.y = pto.y + ACM_MinNBPto.y;
		l_meted(l,&pto);
	}	
	
}


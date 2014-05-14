
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
 #include <fstream>
 #include "levelset.h"
 #include "matematics.h" 
 #include "comun.h"
 #include "listalista.h"

using namespace std;

 //Fa: Factor constante en el caso de la velocidad no asociada a la curvatura
 #define fa 1.0
 #define FILE_FUZZY "centros.txt"
double THRESS_KI = 0.4;
 

//Minima diferencia entre dos levelset para el criterio de parada
double MIN_DIF = 0.01;

float dBeta = 0.5;


//Dimensiones de la Imagen
int IMAGE_HEIGHT;
int IMAGE_WIDTH;
int IMAGE_HEIGHT_old;
int IMAGE_WIDTH_old;
int PresentImageWidth;
int PresentImageHeight;
int W_imag_ori;
int H_imag_ori;

//Puntos del Level Set
lista ptosLevelSet;
//Level Set anterior
lista anteriorptosLevelSet;
//Narrow Band
int NARROW;
//Factor de la curvatura
float EPSILON;
float SIGMA_G;
int H_i;
unsigned char *image;
//Matriz con la imagen convolucionada con la gaussiana
double *gradconvgauss;
//Maximo del gradiente conv con gauss
float m1;
//Minimo del gradiente conv con gauss
float m2;
//Imagen Original
unsigned char *image_ori;
//Matrices del levelset
//Matriz de evolucion del level set
double *fi;
double *fi_old;
bool Gauss = false;
bool circul = false;
//Numero de iteraciones

int iter=1;
//Modo de evolucion: 0 = expansion, 1=compresion
int modo;
float inct;

//Criterio de Parada
float kant;
int pto_fmm_actual=0, pto_fmm_viejo=0;
//Numero de repeticiones del mismo valor del level set
int numrep;
//Para localizar un punto interno del level set
punto2D ptoInterno;
punto2D pto_minNB;
punto2D pto_maxNB;
punto2D pto_minNB_old;


//Parte fuzzy
//Matriz de Etiquetas
unsigned char *mat_fuzzy;
//Etiqueta Mayoritaria
int etiqueta_mayor_ls;
int *matrizEtiq_ls;
bool fuzzy = false;
unsigned char *centros;
int num_clus;
int cluster_int;

void pintar_double(double *mat,char *nombre);

puntof2D* GvfMtx;
int *dist;


/*-----------------------------------------------*/
/*    LEVEL SET                                  */
/*-----------------------------------------------*/
//Determina una nueva narrow band
void determina_nb(lista l,int tama_nb,punto2D *pto_min,punto2D *pto_max){
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
  if((*pto_max).x >= H_imag_ori) (*pto_max).x = H_imag_ori-1;
  (*pto_max).y = (y2+tama_nb);
  if((*pto_max).y >= W_imag_ori) (*pto_max).y = W_imag_ori-1;

}


/*-----------------------------------------------*/
//Genera una imagen binaria en la que el fondo el 0.0 y los puntos del level set 1.0
void imagen_bin(int *imagen_reg){
  punto2D pto;
  int i;
  for(i=0; i<IMAGE_WIDTH*IMAGE_HEIGHT;i++)
    imagen_reg[i] = 0.0;
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)] = 1.0;
  }

}
/*-----------------------------------------------*/
//Realiza una determinación de los signos utilizando un etiquetado de regiones
void llenado_region_globos(){
  int valor,i;

  int *imag2;
  int *etiquetas;
  
  imag2 = (int*)calloc(sizeof(int),IMAGE_WIDTH*IMAGE_HEIGHT);
  etiquetas = (int*)calloc(sizeof(int),IMAGE_WIDTH*IMAGE_HEIGHT);

  imagen_bin(imag2);
  
  etiq_region(imag2,etiquetas);
  free(imag2);
  valor = etiquetas[calcSingleSubscript(0, 0,IMAGE_WIDTH)];
  if(modo == 1){
      for(i=0; i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
    //Para los puntos localizados cambio el signo.
      if(etiquetas[i] == valor){
        fi[i] *= -1;
      }
    }
    }else{
      for(i=0; i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
      if(etiquetas[i] != valor){
        fi[i] *= -1;
      }
    }
  }  
  free(etiquetas);
}  
/*-----------------------------------------------*/
//Determina si un punto pertenece o no al level set
bool puntoLevelSet(int i,int j){
  //Si el valor es menos que cero
  if((fi[calcSingleSubscript(i,j,IMAGE_WIDTH)] < 0.0)&&
  //y alguno de sus 4 vecinos es positivo, entonces es un punto del level set.
        ((fi[calcSingleSubscript(i-1,j,IMAGE_WIDTH)] >= 0.0)||  
         (fi[calcSingleSubscript(i+1,j,IMAGE_WIDTH)] >= 0.0)||
         (fi[calcSingleSubscript(i,j-1,IMAGE_WIDTH)] >= 0.0)||
         (fi[calcSingleSubscript(i,j+1,IMAGE_WIDTH)] >= 0.0)))
          return true;
     return false;
}
/*-----------------------------------------------------------------*/
//Determina los puntos del nuevo level set
void determinar_nuevo_ls_opt(){
  int i,j;
  punto2D pto;
  //Reinicializamos la lista de puntos del LevelSet
  l_dest(&ptosLevelSet);
  ptosLevelSet = l_nuev(sizeof(punto2D));
  //para todos los puntos de la imagen determino aquellos que pertenecen al level set
  for(i=2; i<IMAGE_HEIGHT-2; i++){
    for(j =2; j<IMAGE_WIDTH-2; j++){
      //Si pertenece lo inserto en la nueva lista
       if(puntoLevelSet(i,j)){
           pto.x = i;
           pto.y = j;
           l_meted(ptosLevelSet,&pto);

       }
     }
  }

}
/*--------------------------------------------------------------*/
//    PARTE FUZZY
/*--------------------------------------------------------------*/
/*----------------------------------------------*/
//Determina si un punto pertenece al cluster de interes
bool pertenece_cluster_ls(int x, int y){
  if(image[calcSingleSubscript(x,y,IMAGE_WIDTH)]== centros[cluster_int])
    return true;
  return false;
}
/*-----------------------------------*/
//Calcula la matriz de puntos que pertenecen al cluster
void calc_mat_fuzzy(){
  if(fuzzy){
  for(int i=0; i<IMAGE_HEIGHT; i++){
    for(int j=0; j<IMAGE_WIDTH;j++){
      if(pertenece_cluster_ls(i,j))
        mat_fuzzy[calcSingleSubscript(i,j,IMAGE_WIDTH)]=1.0;
      else
        mat_fuzzy[calcSingleSubscript(i,j,IMAGE_WIDTH)]=0.0;
    }
  }
  fuzzy = !fuzzy;
  }
}
/*-----------------------------------------------------------------*/
void inicializar_fuzzylevelset(){
  unsigned char *imag_temp;
  int *imag_etiq;
  imag_temp = (unsigned char*)malloc(sizeof(unsigned char)*IMAGE_HEIGHT*IMAGE_WIDTH);
  imag_etiq = (int*)malloc(sizeof(int)*IMAGE_HEIGHT*IMAGE_WIDTH);
  matrizEtiq_ls = (int*)calloc(sizeof(int),IMAGE_HEIGHT*IMAGE_WIDTH);
  mat_fuzzy = (unsigned char*)malloc(sizeof(unsigned char)*IMAGE_HEIGHT*IMAGE_WIDTH);
  cargar_centros_fuzzy(&centros,FILE_FUZZY,&num_clus);       
  fuzzyCmeansforzado(image,IMAGE_HEIGHT,IMAGE_WIDTH,imag_temp,centros,num_clus,10);
  for(int i=0; i<IMAGE_HEIGHT*IMAGE_WIDTH; i++) imag_etiq[i] = imag_temp[i];
  etiq_region(imag_etiq,matrizEtiq_ls);
  //etiqueta_mayor_ls = etiq_mayoritaria(ptosLevelSet,matrizEtiq_ls);
    /**/memcpy(image,imag_temp,sizeof(unsigned char)*IMAGE_HEIGHT*IMAGE_WIDTH);
  /**/printf("Etiqueta mayoritaria: %d\n",etiqueta_mayor_ls);
  calc_mat_fuzzy();
  free(imag_temp); 
  free(imag_etiq);
  }

/*--------------------------------------------------------------*/
//CALCULO DEL MAPA DE DISTANCIAS
/*--------------------------------------------------------------*/
/*-----------------------------------------------*/
//Calculo el mapa de distancias con signo aplicando un proceso incremental.
void calc_dist_fi(){
  int i,j,a,b;
   //Calculo la Chamfer Distance
  free(dist);
  chamfer_distance4v(fi,fi);
  dist = (int*)malloc(sizeof(int) * IMAGE_HEIGHT*IMAGE_WIDTH);
  for(int n=0; n<IMAGE_HEIGHT*IMAGE_WIDTH; n++)
	  dist[n] = fi[n];
   //Calculo el signo de las distancias a partir del valor anterior
  if (modo == 0){
    for(i=0; i<IMAGE_HEIGHT_old;i++){
        for(j=0; j<IMAGE_WIDTH_old;j++){
          if(fi_old[calcSingleSubscript(i,j,IMAGE_WIDTH_old)]< 0.0){
            a = i+pto_minNB_old.x-pto_minNB.x;
            b = j+pto_minNB_old.y-pto_minNB.y;
            if(a>=0 && b>=0 && a<IMAGE_HEIGHT && b<IMAGE_WIDTH){
              fi[calcSingleSubscript(a,b,IMAGE_WIDTH)]=
                -fi[calcSingleSubscript(a,b,IMAGE_WIDTH)];
            }
          }
        }
    }
    for(i=0;i<IMAGE_HEIGHT;i++){
      if(fi[calcSingleSubscript(i,0,IMAGE_WIDTH)]<0)
        fi[calcSingleSubscript(i,0,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(i,0,IMAGE_WIDTH)];
      if(fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)]<0)
        fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)];

    }
    for(j=0;j<IMAGE_WIDTH;j++){
      if(fi[calcSingleSubscript(0,j,IMAGE_WIDTH)]<0)
        fi[calcSingleSubscript(0,j,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(0,j,IMAGE_WIDTH)];
      if(fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)]<0)
        fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)];
    }
 
  }else{
    for(i=0; i<IMAGE_HEIGHT_old;i++){
      for(j=0; j<IMAGE_WIDTH_old;j++){         
        if(fi_old[calcSingleSubscript(i,j,IMAGE_WIDTH_old)]< 0.0){
         a = i+pto_minNB_old.x-pto_minNB.x;
         b = j+pto_minNB_old.y-pto_minNB.y;
         if(a>=0 && b>=0 && a<IMAGE_HEIGHT && b<IMAGE_WIDTH){
           fi[calcSingleSubscript(a,b,IMAGE_WIDTH)]=
            -fi[calcSingleSubscript(a,b,IMAGE_WIDTH)];
         }
        }
      }
    }
    for(i=0;i<IMAGE_HEIGHT;i++){
      if(fi[calcSingleSubscript(i,0,IMAGE_WIDTH)]>0)
        fi[calcSingleSubscript(i,0,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(i,0,IMAGE_WIDTH)];
      if(fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)]>0)
        fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(i,IMAGE_WIDTH-1,IMAGE_WIDTH)];

    }
    for(j=0;j<IMAGE_WIDTH;j++){
      if(fi[calcSingleSubscript(0,j,IMAGE_WIDTH)]>0)
        fi[calcSingleSubscript(0,j,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(0,j,IMAGE_WIDTH)];
      if(fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)]>0)
        fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)]=
          -fi[calcSingleSubscript(IMAGE_HEIGHT-1,j,IMAGE_WIDTH)];
    } 
  }
  free(fi_old);  
}
/*-----------------------------------------------*/
void calc_GVF(){	
	GradientVF(image,GvfMtx,0.2);	
}
/*-----------------------------------------------*/
//calculo las distancias basadas en funciones recursivas
void calc_dist_recur(){
  //Calculo la Chamfer Distance
 chamfer_distance4v(fi,fi);
  dist = (int*)malloc(sizeof(int) * IMAGE_HEIGHT*IMAGE_WIDTH);
  for(int n=0; n<IMAGE_HEIGHT*IMAGE_WIDTH; n++)
	  dist[n] = fi[n];
  //Aplicamos el algoritmo de relleno para determinar los signos
  llenado_region_globos();

}
/*-----------------------------------------------*/
//Caclculo la velocidad dependiente de la curvatura o Fg
void calc_veloc_curv(double *fg){  
  curvatura(fi,fg,grad_diffin_norm);
  for(int i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
	  fg[i] *= (-EPSILON);
  }
}


/*-----------------------------------*/
//Ca lculo el factor de parada ki
void calc_ki(double *kiext){
  int i;
  float a;
  gradconvgauss = (double*)calloc(sizeof(double),IMAGE_WIDTH*IMAGE_HEIGHT);
  grad_conv_gauss(image,SIGMA_G,&a,&a,gradconvgauss);
  for(i=0; i<IMAGE_WIDTH*IMAGE_HEIGHT; i++){
    kiext[i]= 1.0/(1.0 + gradconvgauss[i]);
  }
  free(gradconvgauss);
}
/*-----------------------------------*/
//Calculo el factor de parada extendido ^ki
void calc_kiext(double *kiext){
  int i,j;
  double *level;
  punto2D *dist_4sed;
  punto2D pto;
  int tama =IMAGE_HEIGHT*IMAGE_WIDTH;
  double v_grad;
  level = (double*)calloc(sizeof(double),tama);
  
 float a;
  gradconvgauss = (double*)calloc(sizeof(double),IMAGE_WIDTH*IMAGE_HEIGHT);
  grad_conv_gauss(image,SIGMA_G,&a,&a,gradconvgauss);
  
  for(i=0;i<tama;i++)
    kiext[i] = 0.0;
  l_empieza(ptosLevelSet);
  //La velocidad para los puntos del Level set
   while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    //El valor es igual a:
    /*
                  1
    ki =  ------------------
          (1 + |G * I(x,y)|)
    */
    v_grad = 1.0/(1.0 + gradconvgauss[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]);
    if(v_grad<THRESS_KI)
      kiext[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 0.0;
    else
      kiext[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= v_grad;
      
    level[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 1.0;
  }
   //Para hallar los puntos de level set mas cercano a cada punto de la imagen
  //aplico un mapa de distancias 4SSED
   free(gradconvgauss);
  dist_4sed = (punto2D*)calloc(sizeof(punto2D),tama);
  vecin4ssed_distance(level,dist_4sed);
   //Los puntos que no pertenecen al level set tienen una valor igual a la velocidad del
  //punto de la level set mas cercana a ese punto.
  int pos;
  for(i=1;i<IMAGE_HEIGHT-1;i++){
    for(j=1;j<IMAGE_WIDTH-1;j++){
        pos= calcSingleSubscript(i,j,IMAGE_WIDTH);
        kiext[pos] =
          kiext[calcSingleSubscript(i - dist_4sed[pos].x,j - dist_4sed[pos].y,IMAGE_WIDTH)];
    }
  }

  free(level);
  free(dist_4sed);
}
/*-----------------------------------------------*/
void reinicializa_ls(){
  punto2D pto;
   
  //Imagen con la que vamos a tratar que puede ser una subimagen dentro de la original
  free(image);
  //free(GvfMtx);
  //Matriz con la imagen convolucionada con la gaussiana
  lista ptosLevelSet_old;
  ptosLevelSet_old = l_nuev(sizeof(punto2D));
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     pto.x += pto_minNB.x;

     pto.y += pto_minNB.y;
     l_meted(ptosLevelSet_old,&pto);
   }
  l_dest(&ptosLevelSet);
  ptosLevelSet = l_nuev(sizeof(punto2D));
  pto_minNB_old=pto_minNB;
  fi_old = (double*) malloc(sizeof(double)*IMAGE_HEIGHT*IMAGE_WIDTH);
  memcpy(fi_old,fi,(sizeof(double)*IMAGE_HEIGHT*IMAGE_WIDTH));
  free(fi);
  determina_nb(ptosLevelSet_old,NARROW,&pto_minNB,&pto_maxNB);
  l_empieza(ptosLevelSet_old);
  while(l_quedan(ptosLevelSet_old)){
     l_dame(ptosLevelSet_old, &pto);
     pto.x -= pto_minNB.x;
     pto.y -= pto_minNB.y;
     l_meted(ptosLevelSet,&pto);
   }
  l_dest(&ptosLevelSet_old);

  IMAGE_HEIGHT_old=IMAGE_HEIGHT;
  IMAGE_WIDTH_old=IMAGE_WIDTH;
  IMAGE_HEIGHT = pto_maxNB.x -pto_minNB.x;
  IMAGE_WIDTH = pto_maxNB.y -pto_minNB.y;
  image = (unsigned char*)calloc(sizeof(unsigned char),(IMAGE_HEIGHT*IMAGE_WIDTH));
  int i,j,k=0;
  for(i=pto_minNB.x; i<pto_maxNB.x; i++){
    for(j=pto_minNB.y; j<pto_maxNB.y; j++){
     image[k]=image_ori[calcSingleSubscript(i,j,W_imag_ori)];
     k++;
    }
  }

  
  Set_X_Y(IMAGE_WIDTH,IMAGE_HEIGHT);
  //Matriz con la imagen convolucionada con la gaussiana
  gradconvgauss=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
   fi=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
   
  l_empieza(ptosLevelSet);
   while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     fi[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 1.0;
   }
   if(fuzzy){
	   free(mat_fuzzy);
		free(matrizEtiq_ls);
		inicializar_fuzzylevelset();
   }
   
  
    grad_conv_gauss(image,SIGMA_G,&m1,&m2,gradconvgauss);
  //Calculamos la matriz de distancias con signo aplicando recursion
   calc_dist_fi();  
  
}
/*-------------------------------------------------------*/
void reinicializa_lsGVF(){
  punto2D pto;
  //Almacenamos la matriz actual de fi
  fi_old = (double*) malloc(sizeof(double)*IMAGE_HEIGHT*IMAGE_WIDTH);
  memcpy(fi_old,fi,(sizeof(double)*IMAGE_HEIGHT*IMAGE_WIDTH));
  free(gradconvgauss);
  //Matriz del gradiente
  free(fi);
  //Imagen con la que vamos a tratar que puede ser una subimagen dentro de la original
  free(image);
  free(GvfMtx);
  //Matriz con la imagen convolucionada con la gaussiana
  lista ptosLevelSet_old;
  ptosLevelSet_old = l_nuev(sizeof(punto2D));
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     pto.x += pto_minNB.x;

     pto.y += pto_minNB.y;
     l_meted(ptosLevelSet_old,&pto);
   }
  l_dest(&ptosLevelSet);
  ptosLevelSet = l_nuev(sizeof(punto2D));
  pto_minNB_old=pto_minNB;
  determina_nb(ptosLevelSet_old,NARROW,&pto_minNB,&pto_maxNB);
  l_empieza(ptosLevelSet_old);
  while(l_quedan(ptosLevelSet_old)){
     l_dame(ptosLevelSet_old, &pto);
     pto.x -= pto_minNB.x;
     pto.y -= pto_minNB.y;
     l_meted(ptosLevelSet,&pto);
   }
  l_dest(&ptosLevelSet_old);

  IMAGE_HEIGHT_old=IMAGE_HEIGHT;
  IMAGE_WIDTH_old=IMAGE_WIDTH;
  IMAGE_HEIGHT = pto_maxNB.x -pto_minNB.x;
  IMAGE_WIDTH = pto_maxNB.y -pto_minNB.y;
  image = (unsigned char*)calloc(sizeof(unsigned char),(IMAGE_HEIGHT*IMAGE_WIDTH));
  int i,j,k=0;
  for(i=pto_minNB.x; i<pto_maxNB.x; i++){
    for(j=pto_minNB.y; j<pto_maxNB.y; j++){
     image[k]=image_ori[calcSingleSubscript(i,j,W_imag_ori)];
     k++;
    }
  }
  
  Set_X_Y(IMAGE_WIDTH,IMAGE_HEIGHT);
  //Matriz con la imagen convolucionada con la gaussiana
  gradconvgauss=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  fi=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  GvfMtx = (puntof2D*)calloc(sizeof(puntof2D),(IMAGE_HEIGHT*IMAGE_WIDTH));
 
  l_empieza(ptosLevelSet);
   while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     fi[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 1.0;
   }  
   //Calculo del GFV
   calc_GVF();
  //Convolucionamos la imagen con una gausiana
    grad_conv_gauss(image,SIGMA_G,&m1,&m2,gradconvgauss);
  //Calculamos la matriz de distancias con signo aplicando recursion
  calc_dist_fi();  
  
}
/*-----------------------------------------------------------------*/
void filtroLevelSet(){
  punto2D pto,ptoant;
  l_empieza(ptosLevelSet);
  lista filtrols;
  filtrols = l_nuev(sizeof(punto2D));
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    l_meted(filtrols,&pto);
  }

  l_dest(&ptosLevelSet);
  ptosLevelSet = l_nuev(sizeof(punto2D));
  l_empieza(filtrols);
  ptoant.x = ptoant.y = 0;
  while(l_quedan(filtrols)){
    l_dame(filtrols,&pto);
    if((pto.x>0)&&(pto.y>0)&&((pto.x != ptoant.x)||(pto.y != ptoant.y)))
      l_meted(ptosLevelSet,&pto);
    ptoant.x = pto.x;
    ptoant.y = pto.y;
  }
  l_dest(&filtrols);
}
/*-----------------------------------------------------------------*/
void almacena_anteriorls(){
  punto2D pto;
  l_empieza(ptosLevelSet);
  l_dest(&anteriorptosLevelSet);
  anteriorptosLevelSet = l_nuev(sizeof(punto2D));
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    l_meted(anteriorptosLevelSet,&pto);
  }
}


/*-----------------------------------------------------------------*/
void siguiente_LevelSet(double **MatrizEnerg){
  float valor,nuevovalor,kiv,gradv,fgv;
  int i;
  double *fg  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_veloc_curv(fg);
  //Calculamos la matriz asociada al criterio de parada
  if((*MatrizEnerg) != NULL) free(*MatrizEnerg);
  (*MatrizEnerg)  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_kiext((*MatrizEnerg));
  //Calculo del gradiente
  //Matriz del gradiente
	double *grad;

  grad=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  grad_diffin(fi,grad);
  for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){

    valor = fi[i];
    if(fuzzy)
      kiv =mat_fuzzy[i];
    else
      kiv =(*MatrizEnerg)[i];
    fgv = fg[i];
    gradv =grad[i];
    nuevovalor = valor - (inct *(kiv*(fa + fgv) * gradv));
    fi[i]= nuevovalor;
  }
 
  determinar_nuevo_ls_opt();
 
  iter++;
  free(grad);  
  free(fg);
}

/*-----------------------------------------------------------------*/

/*-----------------------------------------------------------------*/
void siguiente_Geodesica(double **MatrizEnerg){
  float valor,nuevovalor,kiv,gradv,gradkiv,fgv;
  int i;
  
  double *fg  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_veloc_curv(fg);
  if((*MatrizEnerg) != NULL) free(*MatrizEnerg);
  (*MatrizEnerg)  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_kiext((*MatrizEnerg));

  //Calculo del gradiente del level set
  	double *grad;

  grad=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  grad_diffin(fi,grad);
   //Calculo el gradiente de kiv
  //ki = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  double *ki = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_ki(ki);
  //grad_Prewitt(ki,gradki);
  double *gradki = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  grad_diffin((*MatrizEnerg),gradki);
  //free(ki);
  for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
    valor = fi[i];
    if(fuzzy)
      kiv =mat_fuzzy[i];
    else
      kiv =(*MatrizEnerg)[i];
    fgv = fg[i];
    gradv =grad[i];
    gradkiv = gradki[i];
    nuevovalor = valor - (inct *((kiv*(fa + fgv)*gradv)+(gradkiv*gradv)));
    fi[i]= nuevovalor;
  }
	determinar_nuevo_ls_opt();
  iter++;
   free(ki);
  free(gradki);
  free(fg);
  free(grad);
}
/*-----------------------------------------------------------------*/
void siguiente_Geodesica_GVF(double **MatrizEnerg){
	
	float valor,nuevovalor,kiv,gradv,gradkiv,fgv,GVFAndGradFi,gvfv;
	int i;
	//Matrices para los gradientes del level set y criterio de parada
	double *gradX = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
	double *gradY = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
	double *GVFGrad = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
	double *geodes= (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
	double *gvf= (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
	double *grad;
	
	grad=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));

	//Velocidad asociada a la curvatura
	double *fg  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_veloc_curv(fg);
  if((*MatrizEnerg) != NULL) free(*MatrizEnerg);
  (*MatrizEnerg)  = (double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  calc_kiext((*MatrizEnerg));

	//Calculo del gradiente del level set
	grad_diffin(fi,gradX,gradY,grad);
	for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
		GVFGrad[i] = ((GvfMtx[i].x * gradX[i]) + (GvfMtx[i].y * gradY[i]));
	}

	VectorNormalized(GVFGrad,GVFGrad);
	
	StoreVectImg(GvfMtx,"GVF.txt");
	guardarImagenTXT(GVFGrad,"GvfGRad.txt");
	unsigned char *iTemp = (unsigned char*) calloc(sizeof(unsigned char),(IMAGE_HEIGHT*IMAGE_WIDTH));
	doubleTounsignedchar(GVFGrad,iTemp);
	guardarImagenPGM(iTemp, "GradGVF.pgm");
	for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
		gradX[i] = sqrt((GvfMtx[i].x * GvfMtx[i].x) + (GvfMtx[i].y*GvfMtx[i].y));
	}
	doubleTounsignedchar(gradX,iTemp);
	guardarImagenPGM(iTemp, "GVF.pgm");
	doubleTounsignedchar(grad,iTemp);
	guardarImagenPGM(iTemp, "Grad.pgm");
	free(iTemp);
	
	for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
		geodes[i] = ((*MatrizEnerg)[i])*(fa+fg[i])*grad[i];
		gvf[i] = GVFGrad[i] * ((*MatrizEnerg)[i]);
	}
	VectorNormalized(geodes,geodes);
	VectorNormalized(gvf,gvf);

	for(i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
		fi[i] = fi[i] - (inct * ((dBeta * geodes[i]) - ((1.0 - dBeta) * gvf[i])));
	}
  free(gradX);
  free(gradY);
  free(GVFGrad);
  free(geodes);
  free(gvf);
  free(grad);
  free(fg);
  //Determina un nuevo level set
 
  determinar_nuevo_ls_opt();
 
  iter++;
}
/*-----------------------------------------------------------------*/
void inicializar_ls_ffm(unsigned char* imag,int mod,float sigma,int narrow,float epsilon,
	lista selec,int width_imag,int height_imag, float umbral,double th_ki){


}
/*-----------------------------------------------------------------*/
void inicializar_level_setGVF(unsigned char* imagen_reg,float increm,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,bool fuz,
  int num_cl,int cl_int,bool interp,float aporte){

	  punto2D pto;
  int i,j,k;
  modo = mod;
  NARROW = narrow;
  EPSILON = epsilon;
  MIN_DIF = umbral;
  SIGMA_G = sigma;
  THRESS_KI = th_ki;
  W_imag_ori = width_imag;
  H_imag_ori = height_imag;
  iter = 1;
  fuzzy = fuz;
  inct = increm;
  num_clus=num_cl;
  cluster_int = cl_int;
  dBeta = aporte;

  determina_nb(selec,narrow,&pto_minNB,&pto_maxNB);


  IMAGE_HEIGHT = pto_maxNB.x -pto_minNB.x;
  IMAGE_WIDTH = pto_maxNB.y -pto_minNB.y;
  image = (unsigned char*)calloc(sizeof(unsigned char),(IMAGE_HEIGHT*IMAGE_WIDTH));
  image_ori = (unsigned char*)calloc(sizeof(unsigned char),(W_imag_ori*H_imag_ori));
  k=0;
  memcpy(image_ori,imagen_reg,sizeof(unsigned char)*W_imag_ori*H_imag_ori);
  for(i=pto_minNB.x; i<pto_maxNB.x; i++){
    for(j=pto_minNB.y; j<pto_maxNB.y; j++){
     image[k]=imagen_reg[calcSingleSubscript(i,j,W_imag_ori)];
     k++;
    }
  }
  

  Set_X_Y(IMAGE_WIDTH,IMAGE_HEIGHT);

  //Matriz con la imagen convolucionada con la gaussiana
  gradconvgauss=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  //Matriz del gradiente

  fi=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));

  GvfMtx = (puntof2D*)calloc(sizeof(puntof2D),(IMAGE_HEIGHT*IMAGE_WIDTH));
 
  ptosLevelSet = l_nuev(sizeof(punto2D));
  lista pto_temp= l_nuev(sizeof(punto2D));
  if(!circ){
	  if(interp){
		  cierre_lineal(selec,l_elementos(selec),pto_temp);
		  l_empieza(pto_temp);
		  while(l_quedan(pto_temp)){
			  l_dame(pto_temp,&pto);
			  pto.x -= pto_minNB.x;
			  pto.y -= pto_minNB.y;
			  l_meted(ptosLevelSet,&pto);
		  }
		  l_dest(&pto_temp);
	  }else{
		  l_dest(&pto_temp);
		  l_empieza(selec);
		  while(l_quedan(selec)){
			  l_dame(selec,&pto);
			  pto.x -= pto_minNB.x;
			  pto.y -= pto_minNB.y;
			  l_meted(ptosLevelSet,&pto);
		  }
	  }
	 
   //Filtramos si hay puntos repetidos
   filtroLevelSet();
   }else{
     circul = true;
     l_empieza(selec);
	 while(l_quedan(selec)){
		 l_dame(selec, &pto);
		 pto.x -= pto_minNB.x;
		 pto.y -= pto_minNB.y;
		 l_meted(ptosLevelSet,&pto);
	 }
   }
   l_empieza(ptosLevelSet);
   while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     fi[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 1.0;
   }
   if(fuzzy)
      inicializar_fuzzylevelset();
  //Calculo del GFV
  calc_GVF();
  //Convolucionamos la imagen con una gausiana
  grad_conv_gauss(image,SIGMA_G,&m1,&m2,gradconvgauss);
  //Calculamos la matriz de distancias con signo aplicando recursion
  calc_dist_recur(); 
}
/*-----------------------------------------------------------------*/
 void inicializar_level_set(unsigned char* imagen_reg,float increm,int mod,float sigma,int narrow,float epsilon,
  lista selec,int width_imag,int height_imag, float umbral, bool circ, double th_ki,bool fuz,
  int num_cl,int cl_int,bool interp){
  punto2D pto;
  int i,j,k;
  modo = mod;
  NARROW = narrow;
  EPSILON = epsilon;
  MIN_DIF = umbral;
  SIGMA_G = sigma;
  THRESS_KI = th_ki;
  W_imag_ori = width_imag;
  H_imag_ori = height_imag;
  iter = 1;
  inct = increm;
  fuzzy = fuz;
  num_clus=num_cl;
  cluster_int = cl_int;

  determina_nb(selec,narrow,&pto_minNB,&pto_maxNB);


  IMAGE_HEIGHT = pto_maxNB.x -pto_minNB.x;
  IMAGE_WIDTH = pto_maxNB.y -pto_minNB.y;
  image = (unsigned char*)calloc(sizeof(unsigned char),(IMAGE_HEIGHT*IMAGE_WIDTH));
  image_ori = (unsigned char*)calloc(sizeof(unsigned char),(W_imag_ori*H_imag_ori));
  k=0;
  memcpy(image_ori,imagen_reg,sizeof(unsigned char)*W_imag_ori*H_imag_ori);
  for(i=pto_minNB.x; i<pto_maxNB.x; i++){
    for(j=pto_minNB.y; j<pto_maxNB.y; j++){
     image[k]=imagen_reg[calcSingleSubscript(i,j,W_imag_ori)];
     k++;
    }
  }
  

  Set_X_Y(IMAGE_WIDTH,IMAGE_HEIGHT);

  //Matriz con la imagen convolucionada con la gaussiana
  gradconvgauss=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
   //Matriz de la velocidad asociada a cada punto de la imagen
  //Matrices del levelset
  fi=(double*)calloc(sizeof(double),(IMAGE_HEIGHT*IMAGE_WIDTH));
  ptosLevelSet = l_nuev(sizeof(punto2D));
  lista pto_temp= l_nuev(sizeof(punto2D));
  if(!circ){
	  if(interp){
		  cierre_lineal(selec,l_elementos(selec),pto_temp);
		  l_empieza(pto_temp);
		  while(l_quedan(pto_temp)){
			  l_dame(pto_temp,&pto);
			  pto.x -= pto_minNB.x;
			  pto.y -= pto_minNB.y;
			  l_meted(ptosLevelSet,&pto);
		  }
		  l_dest(&pto_temp);
	  }else{
		  l_dest(&pto_temp);
		  l_empieza(selec);
		  while(l_quedan(selec)){
			  l_dame(selec,&pto);
			  pto.x -= pto_minNB.x;
			  pto.y -= pto_minNB.y;
			  l_meted(ptosLevelSet,&pto);
		  }
	  }
   //Filtramos si hay puntos repetidos
   filtroLevelSet();
   }else{
     circul = true;
     l_empieza(selec);
	 while(l_quedan(selec)){
		 l_dame(selec, &pto);
		 pto.x -= pto_minNB.x;
		 pto.y -= pto_minNB.y;
		 l_meted(ptosLevelSet,&pto);
	 }
   }
   l_empieza(ptosLevelSet);
   while(l_quedan(ptosLevelSet)){
     l_dame(ptosLevelSet, &pto);
     fi[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]= 1.0;
   }
   if(fuzzy)
      inicializar_fuzzylevelset();
  //Convolucionamos la imagen con una gausiana
  grad_conv_gauss(image,SIGMA_G,&m1,&m2,gradconvgauss);
  //Calculamos la matriz de distancias con signo aplicando recursion
  calc_dist_recur(); 
}
/*-----------------------------------------------------------------*/
bool parada(double *kiext){
  double kexttot = 0;
  float valor = MIN_DIF;
  l_empieza(ptosLevelSet);
  punto2D pto;
	while(l_quedan(ptosLevelSet)){
		l_dame(ptosLevelSet,&pto);

    kexttot+= kiext[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)];
  }
	kexttot /= l_elementos(ptosLevelSet);
  /**/printf("Kant= %f Kextt= %f\n",kant,kexttot);
  if ((kexttot < 0.9)&&((kant - kexttot) < valor)&& ((kant - kexttot)>(-valor)))
  		numrep++;
	else
		numrep = 0;
	kant = kexttot;
  if (numrep >= (NARROW-3)){
    return(true);
	}else{
		return(false);
  }
}
/*-----------------------------------------------------------------*/
bool fin_narrow(){
  punto2D pto;
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    if((pto.x == 1)||(pto.y == 1)||
       (pto.x == IMAGE_HEIGHT-1)||(pto.y == IMAGE_WIDTH-1))
        return true;
  }
  return false;
}

bool fin_narrow_dist(){
  punto2D pto;
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    if(dist[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)] >= NARROW)
        return true;
  }
  return false;
}


/*-----------------------------------------------------------------*/
void liberar_ls(){
  l_dest(&ptosLevelSet);
  //Matriz con la imagen convolucionada con la gaussiana
//  free(gradconvgauss);
  //Matriz del gradiente

  free(fi);
  //Imagen con la que vamos a tratar que puede ser una subimagen dentro de la original
  free(image);
  free(GvfMtx);
  if(fuzzy){
	   free(mat_fuzzy);
		free(matrizEtiq_ls);
		inicializar_fuzzylevelset();
   }
  free(image_ori);
  free(dist);
}
 

/*-----------------------------------------------------------------*/
/*  FUNCIONES PARA DIBUJAR MATRICES ASOCIADAS EL LEVEL SET         */

/*-----------------------------------------------------------------*/
void pintarLevelSet(){
  unsigned char *imagen_reg;
  char nombre[20];
  punto2D pto;
  sprintf(nombre,"LSEXT_%d",iter+100);
	strcat(nombre,".pgm");
  imagen_reg = (unsigned char*)calloc(sizeof(unsigned char) ,(IMAGE_WIDTH * IMAGE_HEIGHT));
  for(int i=0; i<IMAGE_WIDTH*IMAGE_HEIGHT; i++)
      imagen_reg[i] = 0.0;
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)] = 255;
  }
  guardarImagenPGM(imagen_reg, nombre);
  free(imagen_reg);
}
/*-----------------------------------------------------------------*/
void pintar_double(double *mat,char *nombre){
  unsigned char *imagen_reg;
  char cadena[50];

  imagen_reg = (unsigned char*) calloc (sizeof(unsigned char),(IMAGE_WIDTH * IMAGE_HEIGHT));
  doubleTounsignedchar(mat, imagen_reg);
  sprintf(cadena,"%s_%d.pgm",nombre,iter+100);
  guardarImagenPGM(imagen_reg, cadena);

  free(imagen_reg);
}

/*-----------------------------------------------------------------*/
void almacenar_double(double *mat, char* nombre){
  FILE *arch;
  char ext[40];
  sprintf(ext,"%s_%d.txt",nombre,iter+100);
  arch = fopen(ext,"w");
  for(int i=0; i< IMAGE_HEIGHT; i++){
    for(int j=0; j< IMAGE_WIDTH; j++){
    	fprintf(arch,"%.2f\t",mat[calcSingleSubscript(i,j,IMAGE_WIDTH)]);
    }
    fprintf(arch,"\n");
   }
   fclose(arch);
}
/*-----------------------------------------------------------------*/
void almacenar_int(int *mat, char* nombre){
  FILE *arch;
  char cadena[40];
  char ext[10];
  sprintf(ext,"_%d",iter+100);
  sprintf(cadena,nombre);
	strcat(cadena,ext);
	strcat(cadena,".txt");
  arch = fopen(cadena,"w");
  for(int i=0; i< IMAGE_HEIGHT; i++){
    for(int j=0; j< IMAGE_WIDTH; j++){
    	fprintf(arch,"%d\t",mat[calcSingleSubscript(i,j,IMAGE_WIDTH)]);
    }

    fprintf(arch,"\n");
   }
   fclose(arch);
}
void almacenar_uchar(unsigned char *mat, char* nombre){
  FILE *arch;
  char cadena[40];
  char ext[10];
  sprintf(ext,"_%d",iter+100);
  sprintf(cadena,nombre);
	strcat(cadena,ext);
	strcat(cadena,".txt");
  arch = fopen(cadena,"w");
  for(int i=0; i< IMAGE_HEIGHT; i++){
    for(int j=0; j< IMAGE_WIDTH; j++){
    	fprintf(arch,"%d\t",mat[calcSingleSubscript(i,j,IMAGE_WIDTH)]);
    }
    fprintf(arch,"\n");
   }
   fclose(arch);
}

/*-----------------------------------------------------------------*/
void pintar_ls(){
  unsigned char *imagen_reg;
  punto2D pto;
  char nombre[20];
  imagen_reg = (unsigned char*) calloc (sizeof(unsigned char) , (IMAGE_WIDTH * IMAGE_HEIGHT));
  for(int i=0; i<IMAGE_HEIGHT; i++)
    for(int j=0; j<IMAGE_WIDTH; j++)
	    imagen_reg[calcSingleSubscript(i,j,IMAGE_WIDTH)] = image[calcSingleSubscript(i,j,IMAGE_WIDTH)];
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
  	l_dame(ptosLevelSet,&pto);
  	imagen_reg[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)]=0.0;
  }
  sprintf(nombre,"LS_%d",iter+100);
	strcat(nombre,".pgm");
  guardarImagenPGM(imagen_reg, nombre);
  free(imagen_reg);
}
/*-----------------------------------------------------------------*/
void almacena_ls(){
  punto2D  pto;
  FILE *arch;
  arch = fopen("LevelSet.txt","w");

  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
    l_dame(ptosLevelSet,&pto);
    fprintf(arch,"( %d,%d)\n",pto.x, pto.y);
  }
  fclose(arch);
}
/*-----------------------------------------------------------------*/
//Dibuja la curva en una imagen.
void pintar_curva(unsigned char *imagen_reg){
//Inicializa la imagen
  for(int i=0; i<IMAGE_HEIGHT; i++)
    for(int j=0; j< IMAGE_WIDTH; j++)
      imagen_reg[calcSingleSubscript(i,j,IMAGE_WIDTH)] = 0;
  punto2D pto;
  lista l;
  l = l_nuev(sizeof(punto2D));
  crear_circunferencia_cerrada(128,128,60,l);
  l_empieza(l);
  for(int k=0; k<l_elementos(l); k++){
    l_dame(l,&pto);
    imagen_reg[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)] = 255;
  }
  l_dest(&l);
}
/*-----------------------------------------------------------------*/
void actual_level_set_FFM(unsigned char *ls){  
  punto2D pto;
  memcpy(ls,image_ori,sizeof(unsigned char)*H_imag_ori*W_imag_ori);
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
   l_dame(ptosLevelSet,&pto);
   ls[calcSingleSubscript(pto.x,pto.y,W_imag_ori)] = 255.0;
   //   (-1)*(imagen_reg[calcSingleSubscript(pto.x,pto.y,TAMA)]-255.0);
  }  
 }
/*-----------------------------------------------------------------*/
 void actual_level_set(unsigned char *ls){
  punto2D pto;
  memcpy(ls,image_ori,sizeof(unsigned char)*H_imag_ori*W_imag_ori);
  l_empieza(ptosLevelSet);
  while(l_quedan(ptosLevelSet)){
   l_dame(ptosLevelSet,&pto);
   pto.x = pto.x + pto_minNB.x;
   pto.y = pto.y + pto_minNB.y;
   ls[calcSingleSubscript(pto.x,pto.y,W_imag_ori)] = 255.0;
   //   (-1)*(imagen_reg[calcSingleSubscript(pto.x,pto.y,TAMA)]-255.0);
  }  
 }
/*-----------------------------------------------------------------*/
void filtrado_burb(){
  int i,j;
  punto2D pto;
  lista l=l_nuev(sizeof(punto2D));
  int *imagen_reg;
  imagen_reg = (int*)malloc(sizeof(int)*IMAGE_WIDTH*IMAGE_HEIGHT);
  //binarizo la imagen
  imagen_bin(imagen_reg);
  for(i=0; i<IMAGE_HEIGHT-5; i++){
    for(j=0; j<IMAGE_WIDTH-5; j++){
      if((imagen_reg[calcSingleSubscript(i,j,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+1,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+2,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+3,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i,j+4,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+1,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+3,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+4,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+2,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+4,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+1,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+3,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+4,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+1,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+2,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+3,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+4,j+4,IMAGE_WIDTH)]==0.0)&&(
          imagen_reg[calcSingleSubscript(i+1,j+2,IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+1,IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+2,j+3,IMAGE_WIDTH)]==1.0)&&(
          imagen_reg[calcSingleSubscript(i+3,j+2,IMAGE_WIDTH)]==1.0)){
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
    
    l_empieza(ptosLevelSet);
    punto2D pto2;
    while(l_quedan(ptosLevelSet)){
      borrar=false;
      l_dame(ptosLevelSet,&pto);
      l_empieza(l);
      while(l_quedan(l)){
        l_dame(l,&pto2);
        if(pto.x == pto2.x && pto.y == pto2.y)
          borrar = true;
      }
      if(!borrar)
        l_meted(l_temp,&pto);
    }
    l_dest(&ptosLevelSet);
    ptosLevelSet=l_nuev(sizeof(punto2D));
    l_empieza(l_temp);
    while(l_quedan(l_temp)){
      l_dame(l_temp, &pto);
      l_meted(ptosLevelSet,&pto);
    }
	l_dest(&l_temp);
    }
	l_dest(&l);
	
    free(imagen_reg);
}
/*----------------------------------------------*/
bool puntoInternoLS(double* imag,punto2D *pto){
   for(int i=1; i<IMAGE_HEIGHT-1; i++){
      for(int j=1; j<IMAGE_WIDTH-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,IMAGE_WIDTH)] < 0)&&
          (imag[calcSingleSubscript(i-1,j,IMAGE_WIDTH)] == 0)&&
          (imag[calcSingleSubscript(i,j-1,IMAGE_WIDTH)] == 0)){
            (*pto).x = i;
            (*pto).y = j;
           return true;
        }
      }
    }
   return false;
}

/*----------------------------------------------*/
void extraerLeveSet(){
	
	unsigned char* pdImag;
	punto2D pto;
	TLlista lGeneral = new (TCab);
	TRegion trReg,trReg2;
	int i;
	
	pdImag=(unsigned char*)calloc(sizeof(unsigned char),IMAGE_WIDTH*IMAGE_HEIGHT);
	trReg.lContorno = l_nuev(sizeof(punto2D));
	trReg.lPuntos = l_nuev(sizeof(punto2D));
	llInicializa(lGeneral);

	l_empieza(ptosLevelSet);
	while(l_quedan(ptosLevelSet)){
		l_dame(ptosLevelSet,&pto);
		pdImag[calcSingleSubscript(pto.x,pto.y,IMAGE_WIDTH)] = 255;
	}
	extraerSegmentacion(pdImag,lGeneral,255);
	l_dest(&(trReg.lContorno));
	l_dest(&(trReg.lPuntos));
	delete(lGeneral);
}

bool parada_fmm(){
	pto_fmm_actual=l_elementos(ptosLevelSet);
	printf("%d  Ptos LS actual= %d\t Ptos. LS anterior= %d\tDif= %d\n",iter,pto_fmm_actual,pto_fmm_viejo,(pto_fmm_actual-pto_fmm_viejo));
	
	if((abs(pto_fmm_actual-pto_fmm_viejo))<1)
		numrep++;
	else 
		numrep=0;
	pto_fmm_viejo=pto_fmm_actual;
	if(numrep>2){
		numrep=0;
		return true;
	}
	return false;
}

void almacenarLevelSet(char* nameFile){
	punto2D pto;
	FILE *arch;

	arch = fopen(nameFile,"w");

	l_empieza(ptosLevelSet);
	while(l_quedan(ptosLevelSet)){
		l_dame(ptosLevelSet,&pto);
		pto.x = pto.x + pto_minNB.x;
		pto.y = pto.y + pto_minNB.y;
		fprintf(arch,"%d %d\n",pto.x,pto.y);
	}	
	fclose(arch);

}

void actualizar_selec_ls(lista l){
	punto2D pto;
	
	l_empieza(ptosLevelSet);
	while(l_quedan(ptosLevelSet)){
		l_dame(ptosLevelSet,&pto);
		pto.x = pto.x + pto_minNB.x;
		pto.y = pto.y + pto_minNB.y;
		l_meted(l,&pto);
	}	
	
}


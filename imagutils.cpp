


#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <malloc.h>
#include "listagen.h"
#include "matematics.h"
#include "imagutils.h"

using namespace std;

#define SIZE1 100
#define SIZE2 50


enum {DERECHA,IZQUIERDA,ARRIBA,ABAJO,DERABAJO,DERARRIBA,IZQABAJO,IZQARRIBA};

int W_imagen=  512;
int H_imagen = 512;




#define MUESTREO_CIRC 8


//Variable para el calculo de la grad_conv_gauss
float top;
double *cuadrado = NULL;
int num_c;
//Inicializa las dimensiones de la imagen
void Set_X_Y ( int w, int h){
	W_imagen = w;
	H_imagen = h;
}
/*----------------------------------------------*/
/*int calcSingleSubscript(int i,int j,int w){
   return((i * w) + j);
}
*/  
/*----------------------------------------------*/
//Devuelve las dimensiones de la imagen
void Get_X_Y ( int *w, int *h){
	*w = W_imagen;
	*h = H_imagen;
}
/*----------------------------------------------*/
//Carga en la matriz una imagen desde un archivo PGM

int cargarImagenPGM(unsigned char **imag, char *arch){
  ifstream archivo;
  char cabecera[100];   
  archivo.open(arch,fstream::binary);
  if(!archivo) return 1;
  //Quito la cabecera
  for(int i=0; i<2; i++){
	archivo.getline(cabecera,100);
  }
	if(archivo.bad() || archivo.fail()) return 2;
	//Almaceno las dimensiones de la imagen
	archivo>>W_imagen>>H_imagen;
	//descarto la cabecera restante
	archivo.getline(cabecera,100);
	if(archivo.bad() || archivo.fail()) return 2;
	archivo.getline(cabecera,100);
	if(archivo.bad() || archivo.fail()) return 2;
	//Inicializo el vector que contendra la imagen
	*imag = (unsigned char *)calloc(W_imagen*H_imagen,sizeof(unsigned char));
	//Cargo la imagen
	archivo.read((char*)*imag,(W_imagen*H_imagen));
	if(archivo.bad() || archivo.fail()) return 2;
  //cierro el archivo
	archivo.close();
	return 0;

}
/*
int cargarImagenPGM(unsigned char **imag, char *arch){

    FILE* archivo;
    int h,w;
    char cadena[100];
    

    archivo = fopen(arch,"r");
    if(!archivo) return 1;
	
	//Quito la cabecera
	fgets(cadena,100,archivo);
	fgets(cadena,100,archivo);
	printf("Inicio\n");
	//Almaceno las dimensiones de la imagen
	fscanf(archivo,"%d %d\n",&W_imagen,&H_imagen);
	printf("Dimensiones: %d, %d\n",W_imagen,H_imagen);
	//descarto la cabecera restante
	fgets(cadena,100,archivo);

	//Inicializo el vector que contendrá la imagen
	*imag = (unsigned char*)malloc(W_imagen*H_imagen*sizeof(unsigned char));

	//Cargo la imagen
	fread(*imag,sizeof(unsigned char),(W_imagen*H_imagen),archivo);

	//cierro el archivo
	fclose(archivo);
	return 0;
}
*/
/*----------------------------------------------*/
//Carga en la matriz una imagen desde un archivo RAW
int cargarImagenRAW(unsigned char **imag, char *arch,int width, int height){
	ifstream archivo;
  //Almaceno las dimensiones de la imagen
  archivo.open(arch,fstream::binary);
  if(!archivo) return 1;
	W_imagen= width;
    H_imagen= height;
	//Inicializo el vector que contendra la imagen
	*imag = (unsigned char*)calloc(W_imagen*H_imagen,sizeof(unsigned char));
	//Cargo la imagen
	archivo.read((char*)*imag,(W_imagen*H_imagen));
	if(archivo.bad() || archivo.fail()) return 2;
	//cierro el archivo
	archivo.close();
	return 0;
}

/*----------------------------------------------*/
//Guarda la imagen en un archivo PGM
void guardarImagenPGM(unsigned char *imag, char arch[20]){
	ofstream desc;
	desc.open(arch,fstream::binary);	
	desc<<"P5\n"<<"# Created by Antonio Martinez\n"<<W_imagen<<" "<<H_imagen<<"\n"<<"255"<<"\n";
	desc.write((char*)imag,W_imagen*H_imagen);
	desc.close();	
}
/*----------------------------------------------*/
//Guarda la imagen en un archivo RAW
void guardarImagenRAW(unsigned char *imag, char arch[20]){
	ofstream desc;
	desc.open(arch,fstream::binary);
	desc.write((char*)imag,W_imagen*H_imagen);
	desc.close();
}
/*----------------------------------------------*/
//Guarda la imagen en un archivo TXT
void guardarImagenTXT(unsigned char *imag, char arch[20]){
  FILE *desc;
  desc = fopen(arch,"w");
  for(int i=0; i< H_imagen; i++){
    for(int j=0; j< W_imagen; j++){
    	fprintf(desc,"%d\t",(int)imag[calcSingleSubscript(i,j,W_imagen)]);
    }
    fprintf(desc,"\n");
   }
   fclose(desc);
}
/*----------------------------------------------*/
//Guarda una matriz de double en un archivo PGM
void guardarImagenTXT(double *imag, char arch[20]){
  FILE *desc;
  desc = fopen(arch,"w");
  for(int i=0; i< H_imagen; i++){
    for(int j=0; j< W_imagen; j++){
    	fprintf(desc,"%.2f\t",imag[calcSingleSubscript(i,j,W_imagen)]);
    }
    fprintf(desc,"\n");
   }
   fclose(desc);
}
/*----------------------------------------------*/
//Guarda en un archivo TXT las dimensiones de la imagen y los puntos del contorno
void guardarContornoCNT(lista l, char arch[20]){
  FILE *desc;
  desc = fopen(arch,"w");
  punto2D pto;
  fprintf(desc,"C2\n#Creado por Antonio Martinez\n%d %d\n 1\nv %d z 0\n{",W_imagen,H_imagen,l_elementos(l));
  l_empieza(l);
  while(l_quedan(l)){
    l_dame(l,&pto);
    fprintf(desc,"%d %d\n",pto.x, pto.y);
  }
  fprintf(desc,"}");
  fclose(desc);
}
/*----------------------------------------------*/
//Gradiente por diferencias finitas con las dos direcciones del gradiente
void grad_diffin(double *imag,double *gradx, double *grady,double *grad){
    int i,j;
  double valorx, valory;
  int tamaH = H_imagen-1;
  int tamaW = W_imagen-1;
  for(i=1; i<tamaH; i++){
	for(j=1; j<tamaW; j++){
        valorx = (imag[calcSingleSubscript(i+1,j,W_imagen)]-
			imag[calcSingleSubscript(i-1,j,W_imagen)])/2.0;
		valory = (imag[calcSingleSubscript(i,j+1,W_imagen)]-
			imag[calcSingleSubscript(i,j-1,W_imagen)])/2.0;
		 gradx[calcSingleSubscript(i,j,W_imagen)] = valorx;
		grady[calcSingleSubscript(i,j,W_imagen)] = valory;
        grad[calcSingleSubscript(i,j,W_imagen)] = sqrt((valorx*valorx)+(valory*valory));
		}
	}
}

  
/*----------------------------------------------*/
//Gradiente por diferencias finitas
void grad_diffin(double *imag,double *grad){
    int i,j;
  double valorx, valory;
  int tamaH = H_imagen-1;
  int tamaW = W_imagen-1;
  for(i=1; i<tamaH; i++){
	for(j=1; j<tamaW; j++){
        	valorx = (imag[calcSingleSubscript(i+1,j,W_imagen)]-imag[calcSingleSubscript(i-1,j,W_imagen)])/2.0;
		valory = (imag[calcSingleSubscript(i,j+1,W_imagen)]-imag[calcSingleSubscript(i,j-1,W_imagen)])/2.0;
		grad[calcSingleSubscript(i,j,W_imagen)] = sqrt((valorx*valorx)+(valory*valory));
	}
  }
}
/*----------------------------------------------*/
//Transforma una matriz de double en una imagen de unsigned char
void doubleTounsignedchar(double *imag, unsigned char *imag2){
  double rmax, rmin, rng;
  rmax = rmin = imag[0];
  int i;
  int tama = (W_imagen*H_imagen);
  for(i=1; i<tama; i++){
      if (imag[i]>rmax)
        rmax=imag[i];
      if (imag[i]<rmin)
        rmin=imag[i];
  }
    rng=(rmax-rmin);

  for(i=0;i<tama;i++)
       imag2[i] = (unsigned char)(((float)(imag[i]-rmin)/rng)*255.0);
}
/*----------------------------------------------*/
//El negativo de una imagen
void negativo(unsigned char *imag, unsigned char *imagres){
	int i;
  int tama = (W_imagen*H_imagen);
	for(i=1; i<tama; i++)
		imagres[i] = (-1)*(imag[i]-255.0);
}
/*-----------------------------------------------*/
/* CALCULO DEL GRADIENTE CONVOLUCIONADO CON LA   */
/*                 GAUSSIANA                     */
/*-----------------------------------------------*/
 //convolucion en la direccion X
 void convo_vectorx(double *m, float *mask,int len,float *max,float *min,double **resul){
   float *conv;
   int i,j,k,first;
   conv = (float *)calloc(sizeof(float),W_imagen*SIZE1);
   first = 1;
   int tamaW = W_imagen+len;
   int prevW = W_imagen-1;

   *max=-2147483647;
	 *min=2147483647;
   
   
   for(i=0;i<H_imagen;i++){
     for(k=0;k<tamaW;k++){
       conv[k]=0;
       for(j= (int)maximo(0,k-len);j<=(int)minimo(k,prevW);j++)
         conv[k]+= m[calcSingleSubscript(i,j,W_imagen)] * mask[k-j];
    	if(conv[k]<*min) *min = conv[k];
      if(conv[k]>*max) *max = conv[k];

    }
    int len2 = len/2;
    for(j=0;j<W_imagen;j++)
       (*resul)[calcSingleSubscript(i,j,W_imagen)]=conv[j+(len2)];
  }
   free(conv);
}

/*-----------------------------------------------*/
//Convolucion en la direccion Y
void convo_vectory(double *m,float *mask,int len,float *max,float *min, double **resul){
  float *conv;
  int i,j,k,first;
  first = 1;                     
  conv = (float *)calloc(sizeof(float),W_imagen*SIZE1);

   *max=-2147483647;
	 *min=2147483647;
    int tamaH = H_imagen+len;
    int prevH = H_imagen-1;
  
  for(j=0;j<W_imagen;j++){
    for(k=0;k<tamaH;k++){
    	conv[k]=0;                     
	    for(i=(int)maximo(0,k-len);i<=(int)minimo(k,prevH);i++)
	      conv[k]+= m[calcSingleSubscript(i,j,W_imagen)]* mask[k-i];
    	if(conv[k]<*min) *min = conv[k];
      if(conv[k]>*max) *max = conv[k];
    }
    int len2 = len/2;
    for(i=0;i<H_imagen;i++)
      (*resul)[calcSingleSubscript(i,j,W_imagen)]=conv[i+(len2)];
  }
  free(conv);
}
/*-----------------------------------------------*/
 //Calculo una mascara vectorial gaussiana con una desviacion estandar sigma
 //Termina cuando alcanza el 1% del valor maximo
void get_gaussian(float s, float **y, int *len){
  int r,i;

  float gaussian_r;


  r=1;

  gaussian_r=1;
  *len=0;
  (*y)[*len]=gaussian_r;
  while(gaussian_r >= 0.01){
    gaussian_r=exp(-0.5*pow(r,2)/pow(s,2));
    if (gaussian_r>=0.01){
      for (i=*len;i>=0;i--)
        (*y)[i+1]=(*y)[i];
      i=*len;
      (*y)[i+2]=gaussian_r;
      (*y)[0]=gaussian_r;
      r=r+1;
      *len=i+2;
    }/*end if*/
  }/*end while*/
}/*end get_gaussian*/

/*-----------------------------------------------*/
 //Calculo la derivada vectorial de una mascara gaussiana con una desviacion estandar sigma
 //Termina cuando alcanza el 1% del valor maximo
void get_derigaussian(float s, float **y, int *len){
  int r,i;
  float deriv_gaussian_r;
  float max_val;

  r=1;
  max_val=0;
  *len=0;
  deriv_gaussian_r=max_val;
  (*y)[*len]=deriv_gaussian_r;
  while (deriv_gaussian_r>=0.01*max_val){
    deriv_gaussian_r=r/(pow(s,2))*exp(-0.5*pow(r,2)/pow(s,2));
    for(i=*len;i>=0;i--)
      (*y)[i+1]=(*y)[i];
    i=*len;
    (*y)[i+2]= (- deriv_gaussian_r);
    (*y)[0]=deriv_gaussian_r;

    r+=1;
    *len=i+2;
    if (deriv_gaussian_r > max_val)
      max_val=deriv_gaussian_r;
  }/*end while*/
  top = max_val;

}/*end get_derigaussian*/
/*-----------------------------------------------*/
 //Calculo la segunda derivada vectorial de una mascara gaussiana con una desviacion estandar sigma
 //Termina cuando alcanza el 1% del valor maximo
void get_2derigaussian(float s, float **y, int *len){

  int r,i,negative;
  float max_val, d2_gaussian_r; 
  r = 1;

  max_val = - 1/pow(s,2);
  *len = 0;
  d2_gaussian_r = max_val;
  (*y)[*len] = d2_gaussian_r;
  if (max_val<0)
    negative=1;
  while ((d2_gaussian_r>=0.01*max_val)||(negative==1)){
    d2_gaussian_r=((pow(r,2)-pow(s,2))/pow(s,4))*exp(-0.5*pow(r,2)/pow(s,2));
    for (i=*len;i>=0;i--)
      (*y)[i+1]=(*y)[i];
    i=*len;
    (*y)[i+2] = d2_gaussian_r;
    (*y)[0] = d2_gaussian_r;
    r+=1;
    *len=i+2;
    if (d2_gaussian_r > max_val)
      max_val = d2_gaussian_r;
    if (max_val>=0)
      negative=2;
  }/* end while */
  top=max_val;
}/* end get_2derigaussian */
/*-----------------------------------------------*/
//Calculo el gradiente convolucionado con la gaussiana
 void grad_conv_gauss(unsigned char* imag, float sigma, float *max,float *min,double* total_gradient){
   int i,length_gaussian,length_deriv;
   float *gaussian;
   float *deriv_gaussian;
   double* dgaux;
   double* dgauy;
   double* dgaux_gauy;
   double* dgauy_gaux;
   double* m;
   float rmin,rmax;
   int tama = (W_imagen*H_imagen);

   rmax=-2147483647;
	 rmin=2147483647;
   gaussian = (float*) calloc (sizeof(float),SIZE1);
   deriv_gaussian = (float*) calloc (sizeof(float),SIZE2);

   get_gaussian(sigma,&gaussian,&length_gaussian);
   get_derigaussian(sigma,&deriv_gaussian,&length_deriv);

   m = (double*)calloc( sizeof(double) , (W_imagen * H_imagen));
   dgaux = (double*)calloc( sizeof(double) , (W_imagen * H_imagen));
   dgauy = (double*)calloc( sizeof(double) , (W_imagen * H_imagen));
   dgaux_gauy = (double*)calloc( sizeof(double), (W_imagen * H_imagen));

   dgauy_gaux = (double*)calloc( sizeof(double) , (W_imagen * H_imagen));

   for(i=0;i<tama;i++)
          m[i]=(double) imag[i];
  //Calculo la convolucion de la imagen con la derivada de la gaussiana por filas
  convo_vectorx(m,deriv_gaussian,length_deriv,&rmax,&rmin,&dgaux);

  //Calcula el suavizado de la gaussiana en la direccion y
  convo_vectory(dgaux,gaussian,length_gaussian,&rmax,&rmin,&dgaux_gauy);

  //Calcula la convolucion de la imagen con la derivada de la gaussiana por columnas
  convo_vectory(m,deriv_gaussian,length_deriv,&rmax,&rmin,&dgauy);

  //Calcula el suavizado de la gaussiana(y) en la direccion x
  convo_vectorx(dgauy, gaussian, length_gaussian, &rmax, &rmin,&dgauy_gaux);

  for(i=0;i<tama;i++){
      total_gradient[i]=sqrt((dgaux_gauy[i]*dgaux_gauy[i])+(dgauy_gaux[i]*dgauy_gaux[i]));
      if(total_gradient[i]>rmax)
        rmax=total_gradient[i];
      if (total_gradient[i]<rmin)
        rmin=total_gradient[i];
  }  
  *min = rmin;
  *max = rmax;
  free(gaussian);
  free(deriv_gaussian);
  free(m);
  free(dgaux);
  free(dgauy);
  free(dgaux_gauy);
  free(dgauy_gaux);
}
/*-----------------------------------------------*/
//Calculo de una imagen suavizada por la gaussiana
void gaussiana(float sigma, unsigned char *imagen, double *imagsuav){
	double valor,mascara_gausiana[SIZE2][SIZE2];
	float  *gausiana;
	double  *convo;
	int long_gausiana;
    int i,j;

  gausiana =(float *) calloc(sizeof (float),SIZE1);
  convo = (double*) calloc (sizeof(double),W_imagen*H_imagen);
	get_gaussian(sigma,&gausiana,&long_gausiana);

	//solo necesitamos los valores mayores de 0.01
	int start = 0;
	valor = pow(gausiana[0],2);
	while(valor<0.01){
		start++;
		valor = pow(gausiana[start],2);
	}

	//Elimino los valores < 0.01 y construyo la mascara de la gausiana
	for(i=start; i<= long_gausiana - start;i++)
		for(j= start; j<=long_gausiana - start; j++)
			mascara_gausiana[i - start][j-start] = gausiana[i]*gausiana[j];

	//Actualizo el valor del tamaño de la gausiana despues
	//de eliminar los valores < 0.01
	long_gausiana = long_gausiana - (2 * start);

	//Calculo el máximo y el minimo para el reescalado
	int lon = long_gausiana/2;
	double rmax = imagen[0];
	double rmin = rmax;
	double k;

	for(i=lon; i<H_imagen-lon; i++)
		for(j=lon; j<W_imagen-lon; j++){
			k = 0;
			for(int n=0; n<long_gausiana; n++)
				for(int m=0; m<long_gausiana; m++)
					k += mascara_gausiana[n][m] * imagen[calcSingleSubscript((i+n-lon),(j+m-lon),W_imagen)];
			if(k>rmax)
				rmax = k;
			if(k<rmin)
				rmin = k;
			convo[calcSingleSubscript(i,j,W_imagen)] = k;
		}
	double rng = (float)(rmax - rmin);

	for(i=lon; i<H_imagen-lon; i++)
		for(int j=lon; j<W_imagen-lon; j++){
			imagsuav [calcSingleSubscript(i,j,W_imagen)] =
        (((convo[calcSingleSubscript(i,j,W_imagen)] - rmin)/rng)*255.0);
		}
  free(gausiana);
  free(convo);
}
/*-----------------------------------------------*/
/*-----------------------------------------------*/
double minmod(double a, double b){
  if(a > 0.0){
    if(b >0.0){
      return(minimo(fabs(a),fabs(b)));
    }
  }else{
    if(b < 0.0){
      return(-minimo(fabs(a),fabs(b)));
    }
  }
  return 0;
}
/*-----------------------------------------------*/
//Curvatura basada en distancias entre puntos
  void curv_div(double *matriz,double *curv){
   int i,j;
   int pos, posxmas, posxmenos, posymas, posymenos, posxymas, posxymenos;
   int posxmasymenos, posxmenosymas;
   double factor1, factor2, factor3, factor4;

  for(i=1; i<H_imagen-1; i++){
		for(j=1; j<W_imagen-1; j++){
      pos = calcSingleSubscript(i,j,W_imagen);
      posymas = calcSingleSubscript(i,j+1,W_imagen);
      posymenos = calcSingleSubscript(i,j-1,W_imagen);
      posxmas = calcSingleSubscript(i+1,j,W_imagen);
      posxmenos = calcSingleSubscript(i-1,j,W_imagen);
      posxymenos =calcSingleSubscript(i-1,j-1,W_imagen);
      posxymas =calcSingleSubscript(i+1,j+1,W_imagen);
      posxmasymenos =calcSingleSubscript(i+1,j-1,W_imagen);
      posxmenosymas =calcSingleSubscript(i-1,j+1,W_imagen);

      factor1=pow(matriz[posymas]-matriz[pos],2)+
          pow(minmod(matriz[pos]-matriz[posxmenos],matriz[posxmas]-matriz[pos]),2);
      if(factor1 != 0)
        factor1 = (matriz[posymas]-matriz[pos])/sqrt(factor1);
      else
        factor1 = 0.0;

      factor2 = pow(matriz[posymenos]-matriz[pos],2)+

          pow(minmod(matriz[posymenos]-matriz[posxymenos],matriz[posxmasymenos]-matriz[posymenos]),2);
      if(factor2 != 0)
        factor2 =(matriz[posymenos]-matriz[pos])/ sqrt(factor2);
      else
        factor2 = 0.0;

      factor3 = pow(matriz[posxmas]-matriz[pos],2)+
          pow(minmod(matriz[pos]-matriz[posymenos],matriz[posymas]-matriz[pos]),2);
      if(factor3 != 0)
        factor3 =  (matriz[posxmas]-matriz[pos])/sqrt(factor3);
      else
        factor3 = 0.0;

      factor4 =pow(matriz[posxmenos]-matriz[pos],2)+
          pow(minmod(matriz[posxmenos]-matriz[posxymenos],matriz[posxmenosymas]-matriz[posxmenos]),2);
      if(factor4 != 0)
        factor4 = (matriz[posxmenos]-matriz[pos])/sqrt(factor4);
      else
        factor4 = 0.0;
      
        
      curv[pos] = factor1 + factor2 + factor3 + factor4;
    }
  }
}
     
  
/*-----------------------------------------------*/
//Calculo de la curvatura de la matriz Fi
void curvatura(double *matriz,double *curv,void gradiente(double * m, double * gx, double* gy, double* g)){
	double *derivx, *derivy, *derivxx, *derivyy, *derivxy, *derivyx, *grad;
	int pos;

	derivx = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	derivxx = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	derivy = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	derivyy = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	derivxy = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	derivyx = (double*) calloc (sizeof(double) , (W_imagen * H_imagen));
	grad = (double*) malloc (sizeof(double)*(W_imagen * H_imagen));
  //Calculo el gradiente en X, en Y y gradiente de la imagen con mascaras de Prewitt
	gradiente(matriz,derivx,derivy,grad);
  //Calculo el gradiente en X y en Y y gradiente del gradiente en X con mascaras de Prewitt
	gradiente(derivx,derivxx,derivxy,grad);
  //Calculo el gradiente en X y en Y y el gradiente del gradiente en Y con mascaras de Prewitt

	gradiente(derivy,derivyx,derivyy,grad);
  //Calculo la curvatura para los puntos de Fi
  free(grad);
  /*Basada en la siguiente ecuacion:
               2                      2
         Wxx Wy - 2 Wx Wy Wxy + Wyy Wx
  K = - --------------------------------
                       2     2  3/2
                    (Wx  + Wy )
   */

	for(int i=0; i<H_imagen; i++){
    for(int j=0; j<W_imagen; j++){
      pos = calcSingleSubscript(i,j,W_imagen);
      if((derivx[pos]== 0)&&(derivy[pos]==0))
	      curv[pos] = 0.0;
	     else
	      curv[pos] =-((derivxx[pos] * (derivy[pos]*derivy[pos])) - (2 * derivx[pos] * derivy[pos] * derivxy[pos]) +
				(derivyy[pos] * (derivx[pos]*derivx[pos])))/pow(((derivx[pos]*derivx[pos]) + (derivy[pos]*derivy[pos])),1.5);
   }
	}
  //Libero las matrices temporales
  free(derivx);
  free(derivxx);
	free(derivy);
	free(derivyy);
	free(derivxy);
	free(derivyx); 
}
/*----------------------------------------------*/
/*----------------------------------------------*/
double curvatura(punto2D pto, punto2D pto_ant, punto2D pto_sig){
   double result;
   result = sqrt(pow(pto_ant.x - (2.0 * pto.x) +pto_sig.x,2)+
					      pow(pto_ant.y - (2.0 * pto.y) +pto_sig.y,2));

   return result;
}      
/*------------------------------------------------------------*/
//Función accesoria del la función Special Signed Euqlidean Distance (SSED)
void test(int i,int j,int i2,int j2, punto2D *im_dist){

  double dist_x1,dist_y1,dist_x2,dist_y2;
  
//  if((H_imagen>(i+i2))&&((i+i2)>=0)&&(W_imagen>(j+j2))&&((j+j2)>=0)){
    dist_x1 = cuadrado[abs(im_dist[calcSingleSubscript(i+i2,j+j2,W_imagen)].x - i2)];
    dist_y1 = cuadrado[abs(im_dist[calcSingleSubscript(i+i2,j+j2,W_imagen)].y - j2)];
    dist_x2 = cuadrado[abs(im_dist[calcSingleSubscript(i,j,W_imagen)].x)];
    dist_y2 = cuadrado[abs(im_dist[calcSingleSubscript(i,j,W_imagen)].y)];

    if((dist_x1+dist_y1)<(dist_x2+dist_y2)){
      im_dist[calcSingleSubscript(i,j,W_imagen)].x =
        im_dist[calcSingleSubscript(i+i2,j+j2,W_imagen)].x - i2;
      im_dist[calcSingleSubscript(i,j,W_imagen)].y =
        im_dist[calcSingleSubscript(i+i2,j+j2,W_imagen)].y - j2;
    }
//  }
}
/*----------------------------------------------------------*/
//Función que calculo un mapa de distancia con signo
void vecin4ssed_distance(double *im_input,punto2D *im_dist){

  int i,j;
  int tama = W_imagen*H_imagen;
  int tamaH = H_imagen-1;
  int tamaH2 = H_imagen-2;
  int tamaW = W_imagen-1;
  int tamaW2 = W_imagen-2;
  
  if(cuadrado==NULL){  
    cuadrado=(double*)malloc(sizeof(double)*1000);
    for(i=0; i<1000; i++)
      cuadrado[i]=i*i;
  }

  for(i=0; i<tama; i++){
    if(im_input[i]==0){
      im_dist[i].x = 900;
      im_dist[i].y = 900;
    }
    else{
      im_dist[i].x = 0;
      im_dist[i].y = 0;
    }
  }
  //Barrido hacia delante
  for(i=1;i<tamaH;i++){
    for(j=1; j<tamaW; j++){
      test(i,j,-1,0,im_dist);
      test(i,j,0,-1,im_dist);
     } 
    for(j=tamaW2; j>=1; j--)
      test(i,j,1,0,im_dist);      
       
  }
  
  //Barrido hacia atras

  for(i=tamaH2;i>=1; i--){
    for(j=tamaW2; j>=1; j--){
      test(i,j,1,0,im_dist);
      test(i,j,0,1,im_dist);     
    }   
    for(j=1; j<tamaW; j++)
      test(i,j,-1,0,im_dist);        
  } 
}
 
  
/*------------------------------------------------------------*/
//Calculo del mapa de distancias 8 vecinos
void chamfer_distance8v(double *im_input,double *im_output)
 {
  int i,j,k;
  int dist[5],min;
  int d1=1,d2=1;
  for(i=0;i<H_imagen;i++)
    for(j=0;j<W_imagen;j++){
      if(im_input[calcSingleSubscript(i,j,W_imagen)]!=0)
  	    im_output[calcSingleSubscript(i,j,W_imagen)]=0;  /* || max=1000 si se hace la superficie inversa*/
  	 else
	      im_output[calcSingleSubscript(i,j,W_imagen)]=1000; /* || min=0 si se hace la superficie inversa*/
    }

//Iteracion hacia delante
  for(i=1;i<H_imagen-1;i++)
    for(j=1;j<W_imagen-1;j++){
      dist[0]=im_output[calcSingleSubscript(i,j,W_imagen)];
      dist[1]=im_output[calcSingleSubscript(i,j-1,W_imagen)]+d1;

      dist[2]=im_output[calcSingleSubscript(i-1,j-1,W_imagen)]+d2;
      dist[3]=im_output[calcSingleSubscript(i-1,j,W_imagen)]+d1;

      dist[4]=im_output[calcSingleSubscript(i-1,j+1,W_imagen)]+d2;

      min=dist[0];
      im_output[calcSingleSubscript(i,j,W_imagen)]=dist[0];
      for(k=1;k<5;k++)
	      if(dist[k]<min){
          min=dist[k];
          im_output[calcSingleSubscript(i,j,W_imagen)]=dist[k];
	      }


    }
//Iteracion hacia atras
  for(i=H_imagen-2;i>0;i--)
    for(j=W_imagen-2;j>0;j--){
      dist[0]=im_output[calcSingleSubscript(i,j,W_imagen)];
      dist[1]=im_output[calcSingleSubscript(i,j+1,W_imagen)]+d1;
      dist[2]=im_output[calcSingleSubscript(i+1,j+1,W_imagen)]+d2;
      dist[3]=im_output[calcSingleSubscript(i+1,j,W_imagen)]+d1;
      dist[4]=im_output[calcSingleSubscript(i+1,j-1,W_imagen)]+d2;
      min=dist[0];
      im_output[calcSingleSubscript(i,j,W_imagen)]=dist[0];
      for(k=1;k<5;k++)

	      if(dist[k]<min){
          min=dist[k];
          im_output[calcSingleSubscript(i,j,W_imagen)]=dist[k];
	      }
    }
    for(i=0;i<H_imagen;i++)
       im_output[calcSingleSubscript(i,0,W_imagen)]=im_output[calcSingleSubscript(i,1,W_imagen)];
    for(j=0;j<W_imagen;j++)
       im_output[calcSingleSubscript(0,j,W_imagen)]=im_output[calcSingleSubscript(1,j,W_imagen)];
    for(i=0;i<H_imagen;i++)
       im_output[calcSingleSubscript(i,W_imagen-1,W_imagen)]=im_output[calcSingleSubscript(i,W_imagen-2,W_imagen)];
    for(j=0;j<W_imagen;j++)
       im_output[calcSingleSubscript(H_imagen-1,j,W_imagen)]=im_output[calcSingleSubscript(H_imagen-2,j,W_imagen)];
}
/*------------------------------------------------------------*/
//Calculo del mapa de distancias 4 vecinos
void chamfer_distance4v(double *im_input,double *im_output)
 {
  int i,j,k;
  int dist[3],min;
  int d1=1;
  for(i=0;i<H_imagen;i++)
    for(j=0;j<W_imagen;j++){
      if(im_input[calcSingleSubscript(i,j,W_imagen)]!=0)
  	    im_output[calcSingleSubscript(i,j,W_imagen)]=0;  /* || max=1000 si se hace la superficie inversa*/
  	 else
	      im_output[calcSingleSubscript(i,j,W_imagen)]=1000; /* || min=0 si se hace la superficie inversa*/

    }

//Iteracion hacia delante
  for(i=1;i<H_imagen-1;i++)
    for(j=1;j<W_imagen-1;j++){
      dist[0]=im_output[calcSingleSubscript(i,j,W_imagen)];
      dist[1]=im_output[calcSingleSubscript(i,j-1,W_imagen)]+d1;
      dist[2]=im_output[calcSingleSubscript(i-1,j,W_imagen)]+d1;
    

      min=dist[0];
      im_output[calcSingleSubscript(i,j,W_imagen)]=dist[0];
      for(k=1;k<3;k++)
	      if(dist[k]<min){
          min=dist[k];
          im_output[calcSingleSubscript(i,j,W_imagen)]=dist[k];
	      }


    }
//Iteracion hacia atras
  for(i=H_imagen-2;i>0;i--)
    for(j=W_imagen-2;j>0;j--){
      dist[0]=im_output[calcSingleSubscript(i,j,W_imagen)];
      dist[1]=im_output[calcSingleSubscript(i,j+1,W_imagen)]+d1;
      dist[2]=im_output[calcSingleSubscript(i+1,j,W_imagen)]+d1;
  
      min=dist[0];
      im_output[calcSingleSubscript(i,j,W_imagen)]=dist[0];
      for(k=1;k<3;k++)

	      if(dist[k]<min){
          min=dist[k];
          im_output[calcSingleSubscript(i,j,W_imagen)]=dist[k];

	      }
    }
    for(i=0;i<H_imagen;i++)
       im_output[calcSingleSubscript(i,0,W_imagen)]=im_output[calcSingleSubscript(i,1,W_imagen)];
    for(j=0;j<W_imagen;j++)
       im_output[calcSingleSubscript(0,j,W_imagen)]=im_output[calcSingleSubscript(1,j,W_imagen)];
    for(i=0;i<H_imagen;i++)
       im_output[calcSingleSubscript(i,W_imagen-1,W_imagen)]=im_output[calcSingleSubscript(i,W_imagen-2,W_imagen)];
    for(j=0;j<W_imagen;j++)
       im_output[calcSingleSubscript(H_imagen-1,j,W_imagen)]=im_output[calcSingleSubscript(H_imagen-2,j,W_imagen)];
}


/*------------------------------------------------------------*/
//Devuelve el primer punto del interior del level set
bool determinar_pto_interno(int *imag, punto2D *pto){
   for(int i=1; i<H_imagen-1; i++){
      for(int j=1; j<W_imagen-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,W_imagen)] == 0)&&
          (imag[calcSingleSubscript(i-1,j,W_imagen)] == 1)&&
          (imag[calcSingleSubscript(i,j-1,W_imagen)] == 1)){
            (*pto).x = i;
            (*pto).y = j;
           return true;
        }
      }
    }
   return false;
}
/*------------------------------------------------------------*/
//Determinar por recursion 4 conexa los puntos pertenecientes a una region
void evol_punto(int *imagen_reg,int x, int y){
    //Que no se salga de los limites de la imagen
    if((x>=0)&&(x <H_imagen) && (y>=0)&&(y<W_imagen)){
      //Si el punto pertenece a la region
      if(imagen_reg[calcSingleSubscript(x,y,W_imagen)] == 0.0){
          imagen_reg[calcSingleSubscript(x,y,W_imagen)] = 2;
          //Llamadas recursivas por la cola
          evol_punto(imagen_reg,x+1,y);
          evol_punto(imagen_reg,x,y+1);
          evol_punto(imagen_reg,x-1,y);
          evol_punto(imagen_reg,x,y-1);
    
      }
    }
}       
/*-----------------------------------------------*/
//Funcion para interpolar los puntos de una curva de forma lineal
 void interpol_lineal(lista ptos_interp, int x0,int y0,int x1,int y1)
 {
  double x,y,dx,dy,m;
  dy=(double)(y1-y0);
  dx=(double)(x1-x0);
  x=(double)x0;
  y=(double)y0;

  punto2D pto;
  pto.x = x0;

  pto.y = y0;
  //l_meted(ptos_interp, &pto);
  //Si no tienen ninguna coordenada igual
 if(dx!=0.0 && dy!=0.0)

   {
    m=dy/dx;

       if( dx>0 && dy>0 && m<=1.0)
         for(x=x0;x<=x1;x++)
	    {
       pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
	     y=y+m;
	    }


     else if( dx>0 && dy>0 && m>=1.0 )
         for(y=y0;y<=y1;y++)
	    {
	     pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
	     x=x+(1.0/m);
	    }

     else if( dx>0 && dy<0 && m<=-1.0 )
         for(y=y0;y>=y1;y--)
	    {
	     pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);

	     x=x+(1.0/fabs(m));
	    }


     else if( dx>0 && dy<0 && m>=-1.0 )
         for(x=x0;x<=x1;x++)
	    {
	     pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);

	     y=y+m;
	    }

     else if( dx<0 && dy>0 && m>=-1.0 )
         for(x=x0;x>=x1;x--)
	    {

	     pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);

	     y=y+fabs(m);
	    }
     else if( dx<0 && dy>0 && m<=-1.0 )
         for(y=y0;y<=y1;y++)
	    {
	     pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
	     x=x+(1.0/m);
	    }
     else if( dx<0 && dy<0 && m<=1.0 )
         for(x=x0;x>=x1;x--)
	    {
	     pto.x =(int)(0.5+x);

       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
	     y=y-m;
	    }
     else /*if( dx<0 && dy<0 && fabs(dx)<=fabs(dy))*/
         for(y=y0;y>=y1;y--)
	    {
       pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
	     x=x-(1.0/m);
	    }
    }
//Si tienen coordenadas iguales.
  else if( dx==0.0 &&  dy>0.0 )

         for(y=y0;y<=y1;y++)
            {
	        pto.x =(int)(0.5+x);

         pto.y = (int)(0.5+y);
         l_meted(ptos_interp, &pto);

            }

  else if( dx==0.0 &&  dy<0.0 )
         for(y=y0;y>=y1;y--)
            {
	        pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);

       l_meted(ptos_interp, &pto);

            }

  else if( dx>0.0 &&  dy==0.0 )
         for(x=x0;x<=x1;x++)
            {
	        pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
            }
  else if( dx<0.0 &&  dy==0.0 )
         for(x=x0;x>=x1;x--)
            {
	        pto.x =(int)(0.5+x);
       pto.y = (int)(0.5+y);
       l_meted(ptos_interp, &pto);
    }
  }
/*-----------------------------------------------------------------*/
//Funcion que cierra los puntos de una curva usando interpolacion lineal
void cierre_lineal(lista ptos,int numptos, lista curva_cerrada)
 {
  int p;
  punto2D pto, pto_sig;
  punto2D pto_int;
  lista ptos_interp;
  ptos_interp = l_nuev(sizeof(punto2D));
  l_empieza(ptos);
  l_dame(ptos, &pto);
  for(p=0;p<numptos-1;p++){
    l_dame(ptos, &pto_sig);
    interpol_lineal(ptos_interp,pto.x,pto.y,pto_sig.x,pto_sig.y);
	if(!l_vacia(ptos_interp)){
		l_empieza(ptos_interp);
		for(int i =0; i<l_elementos(ptos_interp); i++){
			l_dame(ptos_interp,&pto_int);
			l_meted(curva_cerrada,&pto_int);
		}
		pto.x = pto_sig.x;
		pto.y = pto_sig.y;
	}

  }
  pto.x = pto_sig.x;
  pto.y = pto_sig.y;
  l_empieza(ptos);
  l_dame(ptos, &pto_sig);
  interpol_lineal(ptos_interp,pto.x,pto.y,pto_sig.x,pto_sig.y);
  l_empieza(ptos_interp);

  for(int i =0; i<l_elementos(ptos_interp)-1; i++){

      l_dame(ptos_interp,&pto_int);
      l_meted(curva_cerrada,&pto_int);
    }
  l_dest(&ptos_interp);
 }
/*--------------------------------------------------*/
//Funcion que genera ruido aleatorio
void ruido_aleatorio(unsigned char* imag, unsigned char *imag_res){
	int i,n;


	for (i= 0; i<W_imagen*H_imagen; i++){
		if ((float)rand()/(float)RAND_MAX < 0.2){ 
  		n = (int)(90.0-((float)rand()/(float)RAND_MAX)*180.0);
	  	n = n + imag[i];
		if (n>255)
      n = 255;
		else
      if (n<0)
        n = 0;
    }else{
      n =imag[i];
    }
		imag_res[i] = n;
	}
}


/*-----------------------------------------------------------------*/
//Determina los puntos de un circulo dentro de la imagen
void crear_circunferencia_centrada(int orig_x, int orig_y, float radio,int muestreo,lista puntos){
  punto2D pto;
  for(int i=0; i<(muestreo*2); i++){
    pto.x = orig_x + (radio * cos(i*(Pi/muestreo)));
    pto.y = orig_y + (radio * sin(i*(Pi/muestreo)));
    l_meted(puntos,&pto);
	}
}
/*-----------------------------------------------------------------*/
//crea una circunferencia cerrada, con todos los puntos definidos
void crear_circunferencia_cerrada(int orig_x, int orig_y, float radio,lista curva){

  lista puntos_curva;
  puntos_curva = l_nuev(sizeof(punto2D));
  crear_circunferencia_centrada(orig_x, orig_y, radio,MUESTREO_CIRC, puntos_curva);
  cierre_lineal(puntos_curva,l_elementos(puntos_curva),curva);
  l_dest(&puntos_curva);
}

/*-----------------------------------------------*/
//Crea una lista de puntos formada por circulos en toda la imagen
void multiples_circulos(int num_circ, lista l){
  lista temp_ls;
  punto2D pto;
  num_c = num_circ;
  for(int i=1; i<num_circ+1; i++){
     for(int j=1; j<num_circ+1; j++){
          int res = H_imagen/(num_circ+1);
          int resto =H_imagen%(num_circ+1);
          int radio = res/2;
          int separa = (resto + res)/(num_circ+1);
          temp_ls= l_nuev(sizeof(punto2D));
          crear_circunferencia_cerrada((i*(separa+(radio*2)))-radio,(j*(separa+(radio*2)))-radio,radio,temp_ls);
          l_empieza(temp_ls);
          while(l_quedan(temp_ls)){
            l_dame(temp_ls,&pto);
            l_meted(l,&pto);
          }
          l_dest(&temp_ls);
        }
      }


}
/*-----------------------------------------------*/
//Determina el centro de todos los círculos de la inicialización con círculos
void centro_circulos(lista l){
	punto2D pto;

	for(int i=1; i<num_c+1; i++){
     for(int j=1; j<num_c+1; j++){
		 int res = H_imagen/(num_c+1);
         int resto =H_imagen%(num_c+1);
         int radio = res/2;
         int separa = (resto + res)/(num_c+1);
		 pto.x = (i*(separa+(radio*2)))-radio;
		 pto.y = (j*(separa+(radio*2)))-radio,radio;
		 l_meted(l,&pto);
	 }
	}
}
/*----------------------------------------------*/
/*----------------------------------------------*/
/* rutinas de Fuzzy C means                     */
/*----------------------------------------------*/
int fuzzycmeans_pertenencias(unsigned char*imag, int nDimX, int nDimY, int nNumFases,
				 unsigned char*cnt, double dExp, float **perten);
int fuzzycmeans_centros(float **perten, unsigned char* imag, int nDimX, int nDimY,
				   int nNumFases, double dExp, unsigned char*cnt);
int fuzzycmeans_error(unsigned char*cnt_new, unsigned char*cnt_old,
								int nNumFases, double *pdError);


//Función principal del fuzzy
void fuzzyCmeans(unsigned char *imag,int dim_x, int dim_y,
				 unsigned char *cluster_result, int numclusters,int epsilon){


	int	i,j;
	bool seguir;
	unsigned char *centroides_old;
	unsigned char *centroides_new;

	centroides_old = (unsigned char*)calloc(numclusters, sizeof(unsigned char));
	centroides_new = (unsigned char*)calloc(numclusters, sizeof(unsigned char));
	seguir = true;
	int dwDim = dim_x * dim_y;
  double error;
	//Inicializamos la matriz de pertenencias
	float **m_pert;
	m_pert = (float**)calloc(numclusters, sizeof(float*));
	for(i = 0; i<numclusters; i++)
		m_pert[i]= (float*)calloc(dwDim, sizeof(float));

	//Inicializamos la matriz de centros.
	for(i=0; i<numclusters; i++)
		centroides_old[i]=(unsigned char)((255*i)/(numclusters-1));

	// Calculo de pertenencias
	fuzzycmeans_pertenencias(imag,dim_x,dim_y,numclusters,centroides_old,2.0,m_pert);
	do{
		//Cálculo de centros
		fuzzycmeans_centros(m_pert,imag,dim_x,dim_y,numclusters,2.0,centroides_new);
		//Cálculo de pertenencias
		fuzzycmeans_pertenencias(imag,dim_x,dim_y,numclusters,centroides_new, 2.0,m_pert);
		//Cálculo del error
		fuzzycmeans_error(centroides_new, centroides_old, numclusters, &error);
		printf("Error: %f\n",error);
		if(error<=epsilon){			
			seguir = false;
			break;
		}else{
			for(i=0; i<numclusters;i++){
				centroides_old[i] = centroides_new[i];			
			}
		}
	}while(seguir == true);
  /**/centroides_new[0] = 0.0;
	for(i=0; i<numclusters;i++)
		centroides_old[i] = centroides_new[i];

		unsigned char centro;
		float maxx;
		maxx=0.0;
		for(j=0; j<dwDim;j++){
			for(i=0; i<numclusters; i++){
				if(m_pert[i][j]> maxx){
					maxx =m_pert[i][j];
					centro = (unsigned char)i;
				}
			}
			imag[j] = centroides_new[centro];
			maxx=0.0;
		}
		for(j=0; j<dwDim; j++)
			cluster_result[j] = (int)imag[j];
   /**/
   for(j=0;j<numclusters;j++)
    printf("Centroide %d: %d\n",j,centroides_old[j]);
    /**/
}
/*----------------------------------------------*/
//Función fuzzy pasándole los centroides
void fuzzyCmeansforzado(unsigned char *imag,int dim_x, int dim_y,
				 unsigned char *cluster_result, unsigned char *centroides,int numclusters,
				 int epsilon){

		int	i,j;
	bool seguir;
	unsigned char *centroides_old;
	unsigned char *centroides_new;

	centroides_old = (unsigned char*)calloc(numclusters, sizeof(unsigned char));
	centroides_new = (unsigned char*)calloc(numclusters, sizeof(unsigned char));
	seguir = true;
	int dwDim = dim_x * dim_y;
  double error;
	//Inicializamos la matriz de pertenencias
	float **m_pert;
	m_pert = (float**)calloc(numclusters, sizeof(float*));
	for(i = 0; i<numclusters; i++)
		m_pert[i]= (float*)calloc(dwDim, sizeof(float));

	//Inicializamos la matriz de centros.
	for(i=0; i<numclusters; i++)
		centroides_old[i]=centroides[i];
	// Calculo de pertenencias
	fuzzycmeans_pertenencias(imag,dim_x,dim_y,numclusters,centroides_old,2.0,m_pert);
	do{
		//Cálculo de centros
		fuzzycmeans_centros(m_pert,imag,dim_x,dim_y,numclusters,2.0,centroides_new);
		//Cálculo de pertenencias
		fuzzycmeans_pertenencias(imag,dim_x,dim_y,numclusters,centroides_new, 2.0,m_pert);
		//Cálculo del error
		fuzzycmeans_error(centroides_new, centroides_old, numclusters, &error);
		printf("Error: %f\n",error);
		if(error<=epsilon){
			seguir = false;
			break;
		}else{
			for(i=0; i<numclusters;i++){
				centroides_old[i] = centroides_new[i];				
			}
		}
	}while(seguir == true);
	for(i=0; i<numclusters;i++)
		centroides_old[i] = centroides_new[i];
		unsigned char centro;
		float maxx;
		maxx=0.0;
		for(j=0; j<dwDim;j++){
			for(i=0; i<numclusters; i++){
				if(m_pert[i][j]> maxx){
					maxx =m_pert[i][j];
					centro = (unsigned char)i;
				}
			}
			imag[j] = centroides_new[centro];
			maxx=0.0;
		}
		for(j=0; j<dwDim; j++)
			cluster_result[j] = (int)imag[j];
   /**/
   for(j=0;j<numclusters;j++)
    printf("Centroide %d: %d\n",j,centroides_old[j]);
   /**/
}

/*----------------------------------------------*/
//Calculo de pertenencias fuzzy de cada punto a los centroides
int fuzzycmeans_pertenencias(unsigned char*imag, int nDimX, int nDimY, int nNumFases,
				 unsigned char*cnt, double dExp, float **perten){

	double	dAux1, dAux2, dAux3;
	double	dExponente;
	int	dwIndex;

	int	i, j;
	int	dwDim;

	dwDim = nDimX * nDimY;

	dExponente = 2 / (dExp - 1);
	for(i = 0; i < nNumFases; i++){
		for( dwIndex = 0; dwIndex<dwDim ;  dwIndex++){
			//if(imag[dwIndex] !=(unsigned char)255){
				dAux2 = (imag[dwIndex] > cnt[i])?
					imag[dwIndex] - cnt[i]:cnt[i] - imag[dwIndex];
				dAux1 = 0.0;
				for(j = 0; j < nNumFases; j++){
					dAux3 = (imag[dwIndex] > cnt[j])? imag[dwIndex] - cnt[j]:
															cnt[j] - imag[dwIndex];
					if(dAux3 == 0)
						dAux3 = 0.5; // Evitar división por 0
					if (dExp != 1)
						dAux1 += pow(dAux2 / dAux3, dExponente);
					else
						dAux1 = 0;
				}
			//}
			if(dAux1 != 0)
				perten[i][dwIndex] = (float)(1 / dAux1);
			else
				perten[i][dwIndex] = 1;

		}
	}
	return 0;
}

/*----------------------------------------------*/
//Calculo de los nuevos centroides
int fuzzycmeans_centros(float **perten, unsigned char* imag, int nDimX, int nDimY,
				   int nNumFases, double dExp, unsigned char*cnt){
	double	sdAuxColor;
	double	dAux1, dAux2;
	int	dwIndex;
	int		i;
	int		dwDim;

	dwDim = nDimX * nDimY;
	for(i = 0; i < nNumFases; i++)
	{
		sdAuxColor = 0;
		dAux2 = 0;
		for( dwIndex = 0;dwIndex<dwDim ;  dwIndex++){
			//if(imag[dwIndex] !=(unsigned char)255){
				dAux1 = pow((double)perten[i][dwIndex], dExp);
				sdAuxColor += dAux1*imag[dwIndex];
				dAux2 += dAux1;
			//}

		}
		if (dAux2 == 0)
		{
			cnt[i] = 0;
		}
		else
		{
			cnt[i] = (unsigned char)(sdAuxColor / dAux2);
		}
	}
	return(0);

}
/*----------------------------------------------*/
//Calculo del error
int fuzzycmeans_error(unsigned char*cnt_new, unsigned char*cnt_old,
								int nNumFases, double *pdError){
	double	dAux1;

	int		i;


    (*pdError)=0;
	for(i = 0; i < nNumFases; i++)
	{
		dAux1 = (cnt_new[i] > cnt_old[i])? cnt_new[i] - cnt_old[i]:
												 cnt_old[i] - cnt_new[i];

		(*pdError) += dAux1;
	}

	return(0);
}

/*----------------------------------------------*/
//Devuelve el valor de la etiqueta mayoritaria después de un etiquetado de regiones
int etiq_mayoritaria (lista contorno,int *matriz_etiq){

  punto2D pto;
  int i,max_et=0,max = 0;
  int *vect;
  int *imagen_reg;
  imagen_reg = (int*) calloc (sizeof(double),W_imagen*H_imagen);
  l_empieza(contorno);
  while(l_quedan(contorno)){
    l_dame(contorno,&pto);

    imagen_reg[calcSingleSubscript(pto.x,pto.y,W_imagen)] = 1.0;

  }
  determinar_pto_interno(imagen_reg,&pto);
  evol_punto(imagen_reg,pto.x, pto.y);
  for(i=0; i<W_imagen*H_imagen;i++){
    if(imagen_reg[i]==2){
      if(matriz_etiq[i]>max) max=matriz_etiq[i];
    }
  }
  vect=(int*)calloc(sizeof(int),max);
  for(i=0; i<W_imagen*H_imagen;i++){
    if(imagen_reg[i]==2){
      vect[matriz_etiq[i]]++;
    }
  }

  int pos;
  for(i=0; i<max;i++){
    if(vect[i]>max_et){
       max_et = vect[i];
       pos = i;
    }
  }
  return pos;
}
/*----------------------------------------------*/
//Carga los centroides desde un archivo 
void cargar_centros_fuzzy(unsigned char** vec, char *nombre, int *num_c){
  FILE *f = fopen(nombre,"r");
  int num;
  fscanf(f,"%d",&num);
  *num_c = num;
  (*vec) = (unsigned char*)calloc(sizeof(unsigned char),num);
  for(int i = 0; i<num; i++)
    fscanf(f,"%d\n",&(*vec)[i]);
  fclose(f);
}
/*----------------------------------------------*/
//Determina el Umbral por Otsu
unsigned char umbralOtsu(int histo[256]){
	unsigned char umbral = 0;
	double valor = 0.0;
 	double p1, m1,m2,aux;
	double pi_prima=0,num_ptos=0;
	double histo_acum[256];
	int i;

	//Calculo el total del histograma
	histo_acum[0] = histo[0];

	for(i=1; i<256; i++){
		histo_acum[i] = histo_acum[i-1] + histo[i];
		pi_prima += i * histo[i];
	}

	num_ptos = H_imagen * W_imagen;

	//Inicializo
	p1 = histo_acum[0]/num_ptos;
	m1 = 0.0;
	aux= 0.0;
	//Voy calculando los valores de las variables
	for(i=1; i<255; i++){
		p1 += ((double)histo[i]/(double)num_ptos);
		m1= 0.0;
		m2 = 0.0;
		aux+= i * histo[i];
		m1 = aux/histo_acum[i];
		pi_prima -= i * histo[i];
		m2 = pi_prima /(histo_acum[255] - histo_acum[i]);
		if((p1*(1 - p1)*(pow((m1 - m2),2)))>valor){
			valor = (p1*(1 - p1)*(pow((m1 - m2),2)));
			umbral = i;
		}
	}
	return umbral;
}
/*-----------------------------------------------*/
//Calculo el laplaciano de la Imagen
void laplaciano(unsigned char *imag,double *laplac){
  for(int i=1; i<H_imagen-1; i++){
		for(int j=1; j<W_imagen-1; j++){
		  laplac[calcSingleSubscript(i,j,W_imagen)] =
        (4 * imag[calcSingleSubscript(i,j,W_imagen)]) -
            (imag[calcSingleSubscript(i-1,j,W_imagen)] +
            imag[calcSingleSubscript(i,j-1,W_imagen)]+
            imag[calcSingleSubscript(i+1,j,W_imagen)] +
            imag[calcSingleSubscript(i,j+1,W_imagen)]);
		}
	}
}
/*----------------------------------------------*/
//Binariza una imagen segun el umbral que le pasemos
void binarizaImag(unsigned char *imag, unsigned char umb)
{
	int y;

	for(y = 0; y < H_imagen*W_imagen; y++){
			if(imag[y] < umb)
				imag[y] = 255;
			else
				imag[y] = 0;
		
	}
}
/*----------------------------------------------*/
//Calcula el histograma de una imagen
void histograma(unsigned char* imagen, int histo[256]){
	unsigned char *ori_imag;
  int i;

	ori_imag = imagen;
	for(i=0; i<256; i++)
		histo[i] = 0;
	for(i=0; i<H_imagen*W_imagen; i++){
		histo[(*imagen)]++;
		imagen++;
	}
	imagen = ori_imag;
}    

/*----------------------------------------------*/
int maxHisto(int histo[256],int limite){
	//Devuelve el valor de gris más numeroso menor que el límite
	int max=0,i=0,pos=0;
	for (i=0;i<256 && i<limite;i++){
		if (histo[i]>max){
			max=histo[i];
			pos=i;
		}
	}
	return pos;
}
/*----------------------------------------------*/
int area(unsigned char* imagen, int val){
	long i;
	int areaG=-1;
	int hist[256];
	if (areaG<0){
		areaG=0;
		histograma(imagen,hist);
		for(i=val+1;i<256;i++)
			areaG+=hist[i];
	}
	return areaG;
}
/*----------------------------------------------*/
float entropia(unsigned char* imagen){
	float entro = 0;
	int i;
	float aux=0.0F;
	int hist[256];
	int x,y;
	int LIM = 50;
	int iArea;
	int umbral=maxHisto(hist,LIM);
	iArea = area(imagen,umbral);
	if(entro==0){
		Get_X_Y(&x,&y);
		histograma(imagen,hist);
		aux=(float)iArea;
		for(i=umbral+1;i<256;i++){
			if(((float)(hist[i]/aux))>(float)1e-5)
				entro-=(hist[i]/aux) * (float)log(hist[i]/aux);
		}
	}
	return entro;
}
/*----------------------------------------------*/
bool puntoInterno(double* imag,punto2D *pto){
   for(int i=1; i<H_imagen-1; i++){
      for(int j=1; j<W_imagen-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,W_imagen)] < 0)&&
          (imag[calcSingleSubscript(i-1,j,W_imagen)] == 0)&&
          (imag[calcSingleSubscript(i,j-1,W_imagen)] == 0)){
            (*pto).x = i;
            (*pto).y = j;
           return true;
        }
      }
    }
   return false;
}
/*----------------------------------------------*/
int clase(int y,int *vect,int tama){
	int result = vect[y];
	int pos=y;

	while(result!=0){
		pos=result;
		result = vect[pos];
	}
	return pos;
}
/*----------------------------------------------*/
int reorganizar_etiquetas(int max,lista l,int *m_et){
	int *equivalencias;
   int i,vmax=0;
	punto2D pto;

  equivalencias = (int*)calloc(sizeof(int),max);

	l_empieza(l);
	while(l_quedan(l)){
		l_dame(l,&pto);
    if((equivalencias[pto.x]!= 0)&&(equivalencias[pto.x]<pto.y)){
      equivalencias[pto.y] = equivalencias[pto.x];
    }else{
      equivalencias[pto.x]=clase(pto.y,equivalencias,max);
    }

	}
 	for(i=0; i<H_imagen*W_imagen; i++){
    if(m_et[i]!= 0)
      if(equivalencias[m_et[i]]!= 0)
        m_et[i] = equivalencias[m_et[i]];
	if(m_et[i]>vmax)
		vmax=m_et[i];

	}
	free(equivalencias);
	return(vmax);
}

/*----------------------------------------------*/
int vecinos(int i, int j, int k,lista l, int *m,int *imag){
  int min=k;
  int pos_izq, pos_arr,pos;

  pos= calcSingleSubscript(i,j,W_imagen);
  if(j>0)
    pos_izq = calcSingleSubscript(i,j-1,W_imagen);
  else
    pos_izq =-1;
  if (i>0)
    pos_arr = calcSingleSubscript(i-1,j,W_imagen);
  else
    pos_arr = -1;
  if(pos_izq != -1)
    if((imag[pos_izq]!=1)&&(imag[pos_izq]==imag[pos])){
      if(m[pos_izq]<min)
        min = m[pos_izq];
      l_meted(l,&m[pos_izq]);

    }

  if(pos_arr != -1)

    if((imag[pos_arr]!=1)&&(imag[pos_arr]==imag[pos])){
      if(m[pos_arr]<min)
        min = m[pos_arr];
      l_meted(l,&m[pos_arr]);
    }
  return min;
}
/*----------------------------------------------*/
int vecinos(int i, int j, int k,lista l, int *m,unsigned char *imag){
  int min=k;
  int pos_izq, pos_arr,pos;

  pos= calcSingleSubscript(i,j,W_imagen);
  if(j>0)
    pos_izq = calcSingleSubscript(i,j-1,W_imagen);
  else
    pos_izq =-1;
  if (i>0)
    pos_arr = calcSingleSubscript(i-1,j,W_imagen);
  else
    pos_arr = -1;
  if(pos_izq != -1)
    if((imag[pos_izq]!=1)&&(imag[pos_izq]==imag[pos])){
      if(m[pos_izq]<min)
        min = m[pos_izq];
      l_meted(l,&m[pos_izq]);

    }

  if(pos_arr != -1)

    if((imag[pos_arr]!=1)&&(imag[pos_arr]==imag[pos])){
      if(m[pos_arr]<min)
        min = m[pos_arr];
      l_meted(l,&m[pos_arr]);
    }
  return min;
}
/*----------------------------------------------*/
int etiq_region(int* imagen,int *m_etiq){
	int i,j,k,vmax;
	int pos;
	lista l_equiv;
	k=1;
	l_equiv= l_nuev(sizeof(punto2D));
	for(i=0; i<H_imagen; i++){
		for(j=0; j<W_imagen; j++){
			pos = calcSingleSubscript(i,j,W_imagen);
			if(imagen[pos]!= 1){
				int etiq;
				lista l;
				l = l_nuev(sizeof(int));
				etiq=vecinos(i,j,k,l,m_etiq,imagen);
				if(l_vacia(l)){
					m_etiq[pos]=k;
					k++;
				}else{
					m_etiq[pos] = etiq;
					int temp;
					l_empieza(l);
					while(l_quedan(l)){
						l_dame(l,&temp);
						if(temp != etiq){
							punto2D pto;
							pto.x=temp;
							pto.y=etiq;
							l_meted(l_equiv,&pto);
						}
					}
				}
				l_dest(&l);
			}
		}
	}
	if(!l_vacia(l_equiv))
		vmax=reorganizar_etiquetas(k,l_equiv,m_etiq);
	l_dest(&l_equiv);
	return(vmax);
}

/*----------------------------------------------*/
//crecimiento de regiones 
void crecRegion(int* pdIm,punto2D p2DPto,lista lListaC, lista lListaP){
	punto2D p2DPtoActual;
	lista lPuntosVisitar;
	lPuntosVisitar = l_nuev(sizeof(punto2D));
	l_meted(lPuntosVisitar,&p2DPto);
	while (!l_vacia(lPuntosVisitar)){
		l_empieza(lPuntosVisitar);
		l_sacad(lPuntosVisitar,&p2DPtoActual);
		//Si es un punto del level set
		if(pdIm[calcSingleSubscript(p2DPtoActual.x,p2DPtoActual.y,W_imagen)]==0||
			pdIm[calcSingleSubscript(p2DPtoActual.x,p2DPtoActual.y,W_imagen)]==3){
			l_meted(lListaC,&p2DPtoActual);
			pdIm[calcSingleSubscript(p2DPtoActual.x,p2DPtoActual.y,W_imagen)]=3;
		}
		//Si es parte de la region contenida
		if(pdIm[calcSingleSubscript(p2DPtoActual.x,p2DPtoActual.y,W_imagen)]==-1){
			l_meted(lListaP,&p2DPtoActual);
			pdIm[calcSingleSubscript(p2DPtoActual.x,p2DPtoActual.y,W_imagen)]=2;
			p2DPtoActual.x++;
			l_meted(lPuntosVisitar,&p2DPtoActual);
			p2DPtoActual.x -=2;
			l_meted(lPuntosVisitar,&p2DPtoActual);
			p2DPtoActual.y++;
			p2DPtoActual.x++;
			l_meted(lPuntosVisitar,&p2DPtoActual);
			p2DPtoActual.y -=2;
			l_meted(lPuntosVisitar,&p2DPtoActual);
		}
	}
	l_dest(&lPuntosVisitar);
}
/*----------------------------------------------*/
bool puntoInterno(int* imag,punto2D *pto){
   for(int i=1; i<H_imagen-1; i++){
      for(int j=1; j<W_imagen-1; j++){
        //Busca que sea interno (igual a 0) y que los puntos anteriores sean del
        //level set (iguales a 1)
        if((imag[calcSingleSubscript(i,j,W_imagen)] < 0)&&
          (imag[calcSingleSubscript(i-1,j,W_imagen)] == 0)&&
          (imag[calcSingleSubscript(i,j-1,W_imagen)] == 0)){
            (*pto).x = i;
            (*pto).y = j;
           return true;
        }
      }
    }
   return false;
}
/*----------------------------------------------*/
int extraerSegmentacion(unsigned char *imag, TLlista lGeneral, int ValContour){

	int *imag2;
	punto2D pto;
	TRegion trReg;
	int i,valor;
	int *etiquetas;
	
	
	imag2 = (int*)calloc(sizeof(int),W_imagen*H_imagen);
	etiquetas = (int*)calloc(sizeof(int),W_imagen*H_imagen);
  	
	trReg.lContorno = l_nuev(sizeof(punto2D));
	trReg.lPuntos = l_nuev(sizeof(punto2D));
	llInicializa(lGeneral);
	//Pongo los pixeles del contorno a 1
	for(i = 0;i<W_imagen*H_imagen; i++)
		if(imag[i]==ValContour) imag2[i] = 1;
	//Aplico el etiquetado de regiones
	etiq_region(imag2,etiquetas);
	//Localizo el valor de etiqueta del fondo
	valor = etiquetas[calcSingleSubscript(0,0,W_imagen)];
	for(i=0; i<W_imagen*H_imagen;i++){
	  //Pongo los contornos a 0 y lo demas a 1
	  imag2[i] = !imag2[i];
	  //Para los puntos que no son el fondo cambio el signo.
	  if(etiquetas[i] != valor) imag2[i] *= -1;      
    }
//Construyo la estructura
	while(puntoInterno(imag2,&pto)){
		l_dest(&trReg.lContorno);
		l_dest(&trReg.lPuntos);
		trReg.lContorno=l_nuev(sizeof(punto2D));
		trReg.lPuntos=l_nuev(sizeof(punto2D));
		//Voy haciendo un crecimiento de regiones de las zonas en que los puntos son negativos.
		crecRegion(imag2,pto,trReg.lContorno,trReg.lPuntos);
		llMete(lGeneral,trReg);		
	}
	free(imag2);
	free(etiquetas);
			
	llEmpieza(lGeneral);
	int orden = 0;
	FILE *arch;
	arch = fopen("ListaRegiones.txt","w");
	while(llQuedan(lGeneral)){
		llDame(lGeneral,&trReg);
		fprintf(arch,"\nLevelSet %d\n",orden);
		orden++;
		fprintf(arch,"\tLista de Puntos Region %d: ",l_elementos(trReg.lPuntos));
		l_empieza(trReg.lPuntos);
		while(l_quedan(trReg.lPuntos)){		
			l_dame(trReg.lPuntos,&pto);
			fprintf(arch,"(%d, %d, %d) ",pto.x, pto.y, imag[calcSingleSubscript(pto.x,pto.y,W_imagen)]);
		}
		fprintf(arch,"\n\tLista del Contorno %d: ",l_elementos(trReg.lContorno));
		l_empieza(trReg.lContorno);
		while(l_quedan(trReg.lContorno)){		
			l_dame(trReg.lContorno,&pto);
			fprintf(arch,"(%d, %d) ",pto.x, pto.y);
		}
		l_dest(&(trReg.lContorno));
		l_dest(&(trReg.lPuntos));
	}
	fclose(arch);
	return(orden);
			
}
/*----------------------------------------------*/
void CartesianasToPolares2D(lista lCartesianas,lista lPolares, punto2D p2DCentroGravedad){

	punto2D pto;
	polar2D pol;

	l_empieza(lCartesianas);
	while(!l_vacia(lCartesianas)){
		l_dame(lCartesianas,&pto);
		pol.r = sqrt((pto.x - p2DCentroGravedad.x)*(pto.x - p2DCentroGravedad.x) + 
					 (pto.y - p2DCentroGravedad.y)*(pto.y - p2DCentroGravedad.y));
		pol.a = atan((pto.y - p2DCentroGravedad.y)/(pto.x - p2DCentroGravedad.x));
		l_meted(lPolares,&pol);
	}
}
/*----------------------------------------------*/
void CodigoCadena(unsigned char *imag, int iIntensidad, lista l){
	bool encontrado = false;
	int tama = H_imagen*W_imagen;
	int pos = 0;

	while((encontrado == false)&&(tama>pos)){
		if(imag[pos] == iIntensidad) encontrado = true;
		pos++;
	}
	punto2D p,q;
	p.x = pos/W_imagen;
	p.y = pos%W_imagen;
	int direccion = DERECHA;
	q.x = p.x;
	q.y = p.y;

	do{
		switch(direccion){
		case DERECHA:
			p.y++;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = DERABAJO;
			}
			break;
		case DERABAJO:
			p.y++;
			p.x++;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = ABAJO;
			}
			break;

		case ABAJO:
			p.x++;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = IZQABAJO;
			}
			break;
		case IZQABAJO:
			p.y--;
			p.x++;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = IZQUIERDA;
			}
			break;
		case IZQUIERDA:
			p.y--;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = IZQARRIBA;
			}
			break;
		case IZQARRIBA:
			p.y--;
			p.x--;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = ARRIBA;
			}
			break;
		case ARRIBA:
			p.x--;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = DERARRIBA;
			}
			break;
		case DERARRIBA:
			p.x--;
			p.y++;
			if(imag[calcSingleSubscript(p.x,p.y,W_imagen)] == iIntensidad){				
				l_meted(l,&p);
			}else{
				direccion = DERECHA;
			}
			break;
		}
	}while((p.x!=q.x)||(p.y!=q.y));
	
}
/*----------------------------------------------*/
void SelectorPuntos(lista l1, lista l2, int numPuntos){
	int intervalo = (l_elementos(l1))/numPuntos;
	punto2D pto;

	for(int i=0; i<l_elementos(l1); i++){
		l_dame(l1,&pto);
		if((i%intervalo) == 0){
			l_meted(l2,&pto);
		}
	}	
}
/*----------------------------------------------*/
void DoubleMtxToTXT(char *FileName,double* DoubleMtx){
	FILE *arch = fopen(FileName,"w");
	for(int i=0; i<H_imagen; i++){
		for(int j=0; j<W_imagen;j++){
			fprintf(arch,"%.2f ",DoubleMtx[calcSingleSubscript(i,j,W_imagen)]);
		}
		fprintf(arch,"\n");
	}
	fclose(arch);
}
/*----------------------------------------------*/
bool equalMatrix(puntof2D* A, puntof2D* B){
	for(int i=0; i< H_imagen*W_imagen; i++){
		if((fabs(A[i].x - B[i].x)>0.0)||(fabs(A[i].y - B[i].y)>0.0)){
			return false;
		}
	}
	return(true);
}

/*----------------------------------------------*/
void GradientVF(unsigned char* iImage, puntof2D* mGVF,double alpha){
	double *BMtx;
	double *C1Mtx;
	double *C2Mtx;
	double *GradXMtx;
	double *GradYMtx;
	puntof2D *GVFActual;	
	puntof2D *GVFPrev;
	float max,min;
	float gmax, gmin;

	int i;

	double *dEdgeMap = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	BMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	C2Mtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	C1Mtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GradXMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GradYMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GVFActual = (puntof2D*)calloc(sizeof(puntof2D),H_imagen*W_imagen);	
	GVFPrev = (puntof2D*)calloc(sizeof(puntof2D),H_imagen*W_imagen);	
	
	grad_conv_gauss(iImage, 0.8, &gmax,&gmin,dEdgeMap);


	
	max = min = dEdgeMap[0];
	for(i=1; i< H_imagen*W_imagen; i++){
		if(dEdgeMap[i] > max) max = dEdgeMap[i];
		if(dEdgeMap[i] < min) min = dEdgeMap[i];
	}

	for(i=0; i< H_imagen*W_imagen; i++){
		dEdgeMap[i] = (dEdgeMap[i] - min)/(max-min); 
	}
	
	for(i=0; i< H_imagen*W_imagen; i++){
		dEdgeMap[i] = 1.0 - dEdgeMap[i];
	}
	
	grad_diffin(dEdgeMap,GradXMtx,GradYMtx,BMtx);

	for(i=0; i< H_imagen*W_imagen; i++){
		BMtx[i] = (GradXMtx[i]*GradXMtx[i]) + (GradYMtx[i]*GradYMtx[i]);
		C1Mtx[i] = BMtx[i]*GradXMtx[i];
		C2Mtx[i] = BMtx[i]*GradYMtx[i];
		GVFPrev[i].x = GradXMtx[i];
		GVFPrev[i].y = GradYMtx[i];		
	}

	int vueltas=0;	
	do{
		/**/
			vueltas++;
			printf("Vuetal del GVF: %d\n",vueltas);		
		/**/
		max = 0;
		for(i=W_imagen; i<(H_imagen*(W_imagen-1)); i++){
			GVFActual[i].x =  ((1 - BMtx[i])*GVFPrev[i].x) +
				(alpha*4.0*
				(((GVFPrev[i+1].x+GVFPrev[i-1].x+GVFPrev[i-W_imagen].x+GVFPrev[i+W_imagen].x)*0.25) 
				- GVFPrev[i].x)
				) 
				+ C1Mtx[i];

			GVFActual[i].y = ((1.0 - BMtx[i])*GVFPrev[i].y) + 
				(alpha*4.0*
				(((GVFPrev[i+1].y + GVFPrev[i-1].y + GVFPrev[i-W_imagen].y + GVFPrev[i+W_imagen].y)*0.25) 
				- GVFPrev[i].y)) + C2Mtx[i];		
		}
		for(i=0; i<H_imagen*W_imagen; i++){
			GVFPrev[i].x = GVFActual[i].x;
			GVFPrev[i].y = GVFActual[i].y;
		}
		
	}while(/*(!(equalMatrix(GVFActual,GVFPrev)))&&*/vueltas<(H_imagen));
	for(i=0; i<H_imagen*W_imagen; i++){
		GradXMtx[i] = mGVF[i].x = GVFPrev[i].x;
		GradYMtx[i] = mGVF[i].y = GVFPrev[i].y;
	}
	/*
		DoubleMtxToTXT("GVF_X.txt",GradXMtx);
		DoubleMtxToTXT("GVF_Y.txt",GradYMtx);
	*/
	free(dEdgeMap);
	free(BMtx);
	free(C2Mtx);
	free(C1Mtx);
	free(GradXMtx);
	free(GradYMtx);
	free(GVFActual);
	free(GVFPrev);
}
/*----------------------------------------------*/
void GradientVF_Paragios(unsigned char* iImage, puntof2D* mGVF,double alpha){
	double *BMtx;
	double *C1Mtx;
	double *C2Mtx;
	double *GradXMtx;
	double *GradYMtx;
	puntof2D *GVFActual;	
	puntof2D *GVFPrev;
	float max,min;
	float gmax, gmin;

	int i;

	double *dEdgeMap = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	BMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	C2Mtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	C1Mtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GradXMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GradYMtx = (double*)calloc(sizeof(double),H_imagen*W_imagen);
	GVFActual = (puntof2D*)calloc(sizeof(puntof2D),H_imagen*W_imagen);	
	GVFPrev = (puntof2D*)calloc(sizeof(puntof2D),H_imagen*W_imagen);	
	
	grad_conv_gauss(iImage, 1.0, &gmax,&gmin,dEdgeMap);
	
	max = min = dEdgeMap[0];
	for(i=1; i< H_imagen*W_imagen; i++){
		if(dEdgeMap[i] > max) max = dEdgeMap[i];
		if(dEdgeMap[i] < min) min = dEdgeMap[i];
	}

	for(i=0; i< H_imagen*W_imagen; i++){
		dEdgeMap[i] = (dEdgeMap[i] - min)/(max-min); 
	}
	grad_diffin(dEdgeMap,GradXMtx,GradYMtx,BMtx);
	for(i=0; i< H_imagen*W_imagen; i++){
		BMtx[i] = ((GradXMtx[i]*GradXMtx[i]) + (GradYMtx[i]*GradYMtx[i])) * dEdgeMap[i];
		C1Mtx[i] = BMtx[i]*GradXMtx[i];
		C2Mtx[i] = BMtx[i]*GradYMtx[i];
		GVFPrev[i].x = GradXMtx[i];
		GVFPrev[i].y = GradYMtx[i];	
	}
	int vueltas=0;	
	do{
		/**/
			vueltas++;
			printf("Vuetal del GVF: %d\n",vueltas);		
		/**/
		max = 0;
		for(i=W_imagen; i<(H_imagen*(W_imagen-1)); i++){
			GVFActual[i].x =  ((1 - BMtx[i])*GVFPrev[i].x) +
				(alpha*4.0*
				(((GVFPrev[i+1].x+GVFPrev[i-1].x+GVFPrev[i-W_imagen].x+GVFPrev[i+W_imagen].x)*0.25) 
				- GVFPrev[i].x)
				) 
				+ C1Mtx[i];

			GVFActual[i].y = ((1.0 - BMtx[i])*GVFPrev[i].y) + 
				(alpha*4.0*
				(((GVFPrev[i+1].y + GVFPrev[i-1].y + GVFPrev[i-W_imagen].y + GVFPrev[i+W_imagen].y)*0.25) 
				- GVFPrev[i].y)) + C2Mtx[i];		
		}
		for(i=0; i<H_imagen*W_imagen; i++){
			GVFPrev[i].x = GVFActual[i].x;
			GVFPrev[i].y = GVFActual[i].y;
		}
		
	}while(/*(!(equalMatrix(GVFActual,GVFPrev)))&&*/(vueltas<H_imagen));
	for(i=0; i<H_imagen*W_imagen; i++){
		GradXMtx[i] = mGVF[i].x = GVFPrev[i].x;
		GradYMtx[i] = mGVF[i].y = GVFPrev[i].y;
	}
	/*
		DoubleMtxToTXT("GVF_X.txt",GradXMtx);
		DoubleMtxToTXT("GVF_Y.txt",GradYMtx);
	*/
	free(dEdgeMap);
	free(BMtx);
	free(C2Mtx);
	free(C1Mtx);
	free(GradXMtx);
	free(GradYMtx);
	free(GVFActual);
	free(GVFPrev);
}
/*----------------------------------------------*/
void DrawVectImag(puntof2D* ImgVectMtx, unsigned char* Img){
	int i,j,m,n,k;
	lista lInterpPoints = l_nuev(sizeof(punto2D));
	punto2D pto;
	for(i=0; i<H_imagen*W_imagen; i++)
		Img[i] = 0;

	for(i=0; i<H_imagen; i=i+8){
		for(j=0; j<W_imagen; j=j+8){
			m = i + (ImgVectMtx[calcSingleSubscript(i,j,W_imagen)].x * 10.0);
			n = j + (ImgVectMtx[calcSingleSubscript(i,j,W_imagen)].y * 10.0);
			if((m>0)&&(m<H_imagen)&&(n>0)&&(n<W_imagen)){
				if((i==m)&&(j==n)){
					pto.x = m;
					pto.y = n;
					l_meted(lInterpPoints,&pto);
				}else{
					interpol_lineal(lInterpPoints,i,j,m,n);	
				}
				l_empieza(lInterpPoints);
				while(l_quedan(lInterpPoints)){
					l_dame(lInterpPoints,&pto);
					Img[calcSingleSubscript(pto.x,pto.y,W_imagen)] = 128;
				}
				for(k=-1; k<2;k++){
					Img[calcSingleSubscript(m+k,n+k,W_imagen)] = 255;
					Img[calcSingleSubscript(m+k,n,W_imagen)] = 255;
					Img[calcSingleSubscript(m,n+k,W_imagen)] = 255;
					Img[calcSingleSubscript(m+k,n-k,W_imagen)] = 255;
				}

			}
		}		
	}
	l_dest(&lInterpPoints);
}
/*----------------------------------------------*/
void DrawVectImag(punto2D* ImgVectMtx, unsigned char* Img){
	int i,j,m,n,k;
	lista lInterpPoints = l_nuev(sizeof(punto2D));
	punto2D pto;
	for(i=0; i<H_imagen*W_imagen; i++)
		Img[i] = 0;

	for(i=0; i<H_imagen; i=i+8){
		for(j=0; j<W_imagen; j=j+8){
			m = i - (ImgVectMtx[calcSingleSubscript(i,j,W_imagen)].x * 10.0);
			n = j - (ImgVectMtx[calcSingleSubscript(i,j,W_imagen)].y * 10.0);
			if((m>0)&&(m<H_imagen)&&(n>0)&&(n<W_imagen)){
				if((i==m)&&(j==n)){
					pto.x = m;
					pto.y = n;
					l_meted(lInterpPoints,&pto);
				}else{
					interpol_lineal(lInterpPoints,i,j,m,n);	
				}
				l_empieza(lInterpPoints);
				while(l_quedan(lInterpPoints)){
					l_dame(lInterpPoints,&pto);
					Img[calcSingleSubscript(pto.x,pto.y,W_imagen)] = 128;
				}
				for(k=-1; k<2;k++){
					Img[calcSingleSubscript(m+k,n+k,W_imagen)] = 255;
					Img[calcSingleSubscript(m+k,n,W_imagen)] = 255;
					Img[calcSingleSubscript(m,n+k,W_imagen)] = 255;
					Img[calcSingleSubscript(m+k,n-k,W_imagen)] = 255;
				}

			}
		}		
	}
	l_dest(&lInterpPoints);
}
/*----------------------------------------------*/
void StoreVectImg(puntof2D *P2dfMtx,char *name){
	FILE* arch;
	int i,j;
	arch=fopen(name,"w");
	for(i=0; i<H_imagen; i++){
		for(j=0; j<W_imagen; j++){
			fprintf(arch,"(%.2f,%.2f)\t",P2dfMtx[calcSingleSubscript(i,j,W_imagen)].x,
				P2dfMtx[calcSingleSubscript(i,j,W_imagen)].y);
		}
		fprintf(arch,"\n");
	}
	fclose(arch);
}
/*----------------------------------------------*/
void StoreVectImg(punto2D *P2dfMtx,char *name){
	FILE* arch;
	int i,j;
	arch=fopen(name,"w");
	for(i=0; i<H_imagen; i++){
		for(j=0; j<W_imagen; j++){
			fprintf(arch,"(%d,%d)\t",P2dfMtx[calcSingleSubscript(i,j,W_imagen)].x,
				P2dfMtx[calcSingleSubscript(i,j,W_imagen)].y);
		}
		fprintf(arch,"\n");
	}
	fclose(arch);
}
/*----------------------------------------------*/
void grad_diffin_norm(double *imag,double *gradx, double *grady,double *grad){

	int i;
	double max,min;

	grad_diffin(imag,gradx,grady,grad);
	
	max = min = gradx[0];
	for(i=1; i<H_imagen*W_imagen; i++){
		if(gradx[i]> max) max = gradx[i];
		if(grady[i]> max) max = grady[i];		
		if(gradx[i]< min) min = gradx[i];
		if(grady[i]< min) min = grady[i];		
	}
	for(i=0; i<H_imagen*W_imagen; i++){
		gradx[i] = gradx[i] - min/(max-min) ;
		grady[i] = grady[i] - min/(max-min);		
	}
	VectorNormalized(grad,grad);

}
/*----------------------------------------------*/
void VectorModule(puntof2D *P2dfMtx, double *ModuleMtx){
	for(int i=0; i<H_imagen*W_imagen; i++){
		ModuleMtx[i] = sqrt((P2dfMtx[i].x*P2dfMtx[i].x)+(P2dfMtx[i].y*P2dfMtx[i].y));
	}
}
/*----------------------------------------------*/
void VectorModule(punto2D *P2dfMtx, double *ModuleMtx){
	for(int i=0; i<H_imagen*W_imagen; i++){
		ModuleMtx[i] = sqrt((P2dfMtx[i].x*P2dfMtx[i].x)+(P2dfMtx[i].y*P2dfMtx[i].y));
	}
}
/*----------------------------------------------*/
void ScaleTransform(double* iImageOri,double intervA,double intervB,double* iImageFin){
	
	double max, min; 
	int i;
	max = min = iImageOri[0];
	for(i=1; i<H_imagen*W_imagen; i++){
		if(iImageOri[i]> max) max = iImageOri[i];
		if(iImageOri[i]< min) min = iImageOri[i];		
	}
	double EscalaA = max - min;
	double EscalaB = intervB - intervA;
	for(i=0; i<H_imagen*W_imagen; i++){
		iImageFin[i] = (((iImageOri[i] - min)*EscalaB)+(intervA * EscalaA))/EscalaA;
	}
}
/*----------------------------------------------*/
void VectorNormalized(puntof2D* P2dfMtxIni,puntof2D* P2dfMtxFin){

	double max, min; 
	int i;
	max = min = P2dfMtxIni[0].x;
	for(i=1; i<H_imagen*W_imagen; i++){
		if(P2dfMtxIni[i].x> max) max = P2dfMtxIni[i].x;
		if(P2dfMtxIni[i].y> max) max = P2dfMtxIni[i].y;
		if(P2dfMtxIni[i].x< min) min = P2dfMtxIni[i].x;		
		if(P2dfMtxIni[i].y< min) min = P2dfMtxIni[i].y;		
	}
	double EscalaA = max - min;
	for(i=0; i<H_imagen*W_imagen; i++){
		P2dfMtxFin[i].x = (P2dfMtxIni[i].x - min)/EscalaA;
		P2dfMtxFin[i].y = (P2dfMtxIni[i].y - min)/EscalaA;

	}
}
/*----------------------------------------------*/
void VectorNormalized(double* MtxIni,double* MtxFin){
	double max, min; 
	int i;
	max = min = MtxIni[0];
	for(i=1; i<H_imagen*W_imagen; i++){
		if(MtxIni[i]> max) max = MtxIni[i];
		if(MtxIni[i]< min) min = MtxIni[i];		
	}
	double EscalaA = max - min;
	for(i=0; i<H_imagen*W_imagen; i++){
		MtxFin[i] = (MtxIni[i] - min)/EscalaA;
	}
}
/*----------------------------------------------*/
TMatrizI* CreateMatrizI(TMatrizI *Matrix){
	TMatrizI *MatrizFin = (TMatrizI*) malloc(sizeof(TMatrizI));
	MatrizFin->height = Matrix->height;
	MatrizFin->width = Matrix->width;
	MatrizFin->matrix = (int*)malloc(sizeof(int)*Matrix->height * Matrix->width);
	memcpy(MatrizFin->matrix,Matrix->matrix,(sizeof(int)*Matrix->height * Matrix->width));
	return MatrizFin;
}
/*----------------------------------------------*/
TMatrizI* CreateMatrizI(TMatrizD* Matrix){
	TMatrizI *MatrizFin = (TMatrizI*) malloc(sizeof(TMatrizI));
	MatrizFin->height = Matrix->height;
	MatrizFin->width = Matrix->width;
	MatrizFin->matrix = (int*)malloc(sizeof(int)*Matrix->height * Matrix->width);
	for(int i=0; i<Matrix->height * Matrix->width; i++) MatrizFin->matrix[i] = Matrix->matrix[i];
	return MatrizFin;
}
/*----------------------------------------------*/
TMatrizI* CreateMatrizI(int height, int width){
	TMatrizI *MatrizFin = (TMatrizI*) malloc(sizeof(TMatrizI));
	MatrizFin->height = height;
	MatrizFin->width = width;
	MatrizFin->matrix = (int*)calloc(sizeof(int),height * width);
	return MatrizFin;
}
/*----------------------------------------------*/
TMatrizD* CreateMatrizD(TMatrizD* Matrix){
	TMatrizD *MatrizFin = (TMatrizD*) malloc(sizeof(TMatrizD));
	MatrizFin->height = Matrix->height;
	MatrizFin->width = Matrix->width;
	MatrizFin->matrix = (double*)malloc(sizeof(double)*Matrix->height * Matrix->width);
	memcpy(MatrizFin->matrix,Matrix->matrix,(sizeof(double)*Matrix->height * Matrix->width));
	return MatrizFin;
}
/*----------------------------------------------*/
TMatrizD* CreateMatrizD(int height, int width){
	TMatrizD *MatrizFin = (TMatrizD*) malloc(sizeof(TMatrizD));
	MatrizFin->height = height;
	MatrizFin->width = width;
	MatrizFin->matrix = (double*)calloc(sizeof(double),height * width);
	return MatrizFin;
}
/*----------------------------------------------*/
void DestroyMatrix(TMatrizI **Matrix){
	free((*Matrix)->matrix);
	free(*Matrix);
	(*Matrix)=NULL;
}
/*----------------------------------------------*/
void DestroyMatrix(TMatrizD **Matrix){
	free((*Matrix)->matrix);
	free(*Matrix);
	(*Matrix)=NULL;
}
/*----------------------------------------------*/
TImagen* CreateImagen(int height, int width){
	TImagen *MatrizFin = (TImagen*) malloc(sizeof(TImagen));
	MatrizFin->height = height;
	MatrizFin->width = width;
	MatrizFin->matrix = (unsigned char*)calloc(sizeof(unsigned char),height * width);
	return MatrizFin;
}
/*----------------------------------------------*/
TImagen* CreateImagen(TImagen *Imag){
	TImagen *MatrizFin = (TImagen*) malloc(sizeof(TImagen));
	MatrizFin->height = Imag->height;
	MatrizFin->width = Imag->width;
	MatrizFin->matrix = (unsigned char*)calloc(sizeof(unsigned char),Imag->height * Imag->width);
	memcpy(MatrizFin->matrix,Imag->matrix,(sizeof(unsigned char)*Imag->height * Imag->width));
	return MatrizFin;
}
/*----------------------------------------------*/
void DestroyImagen(TImagen **Imag){
	free((*Imag)->matrix);
	free(*Imag);
	(*Imag)=NULL;
}
/*----------------------------------------------*/
int maximov(int *vector,int dim){

	int max = vector[0];
	for(int i=0; i<dim; i++)
		if(vector[i]>max) max=vector[i];
	return(max);
}
/*----------------------------------------------*/
int numPicos(lista lPicos, int vecindad, double umbral, unsigned char* imag){
	int histogram[256];
	double histNorm[256];
	bool pico = false;
	int i,n,nPicos=0;
	histograma(imag,histogram);
	double maxim = maximov(histogram,256);
	for(n=0; n<256; n++)
		histNorm[n] = ((double)histogram[n])/maxim;
	/**/
	FILE *f;
	f = fopen("Histograma.txt","w");
	for(n=0; n<256; n++)
		fprintf(f,"%d \t %.3f\n",n,histNorm[n]);

	fclose(f);
	/**/
	for(i=vecindad; i<(255-vecindad); i++){
		if(histNorm[i] > umbral){
			pico = true;
			for(int n=vecindad; n>0; n--){				
				if(histogram[i - n]>histogram[i])
					pico=false;
				if(histogram[i + n]>histogram[i])
					pico=false;
			}
			if(pico == true){
				nPicos++;
				l_meted(lPicos,&i);
				i += vecindad+1;
			}
		}
	}
	return(nPicos);
}
/*----------------------------------------------*/
void RegionesInteres(lista lPicos, int vecindad, double  umbral, unsigned char* imag){
	
	int t;
	int nP = numPicos(lPicos, vecindad, umbral,imag);
	if(nP>2){
		l_empieza(lPicos);
		l_sacai(lPicos,&t);
	}
	
	/**/
	l_empieza(lPicos);
	while(l_quedan(lPicos)){
		l_dame(lPicos,&t);
		printf("Pico: %d\n",t);
	}
	/**/
}

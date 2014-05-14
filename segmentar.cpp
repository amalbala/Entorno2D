/*--------------------------------------------------*/
/*					UCLM - ISA						*/
/*				E.T.S.I. INDUSTRIALES				*/
/*													*/
/*	AvnD. Camilo Jose Cela s/n - 13071 Ciudad Real	*/
/*													*/
/*	Proyecto: Funciones Separables					*/
/*	Investigador Principal: Gloria Bueno Garcia		*/
/*  Programadores: Antonio Martinez					*/
/*													*/
/*	Fecha: Junio - 2005								*/
/*--------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "comun.h"
#include "imagutils.h"
#include "snakes.h"
#include "levelset.h"





enum {PDM,LS,GEOD,EXPAN,CONTR};

double *MatrizEnerg = NULL;

/*---------------- SUB MENU -------------------------------*/
void cargarPuntos(char arch[40],lista list){
  FILE *archivo;
  archivo = fopen(arch,"r");
  punto2D pto;
  int a,b;
  
  Get_X_Y(&a,&b);
  while(fscanf(archivo,"%d %d\n",&pto.x,&pto.y)!= EOF){
	pto.x = b - pto.x;
	l_meted(list,&pto);
  }
  fclose(archivo);
}
/*---------------- SUB MENU -------------------------------*/
int main(int argc, char *argv[]){
  //Tiempo Computacional
  clock_t tiempo_total;
  clock_t tiempo1, tiempo2;
  int opcion;
  unsigned char *imag;
  int a,b;
  float increm;
  int nb;
  lista puntosselec;
  puntosselec=l_nuev(sizeof(punto2D));
  tiempo_total = tiempo1 = tiempo2 = 0;
  if(argc == 1){
	  printf(CABECERA_UCLM);
    printf("segmentar -P <imagen_orginal> <imagen_resultado> <contorno_inicial>\n\t<alfa> <beta> <gamma> <w_lineas> <w_bordes> <w_term> <w_presion>\n");
    printf("segmentar -L|G -C|-E  <imagen_orginal> <imagen_resultado> <contorno_inicial>\n\t<sigma> <Narrow_Band> <Increm. t> <Epsilon> <Thresshold>\n");
  }else{
    switch(argv[1][1]){
      case 'P':
        if (argc != 12){
          printf("Numero de parametros incorrecto. El correcto es 11.\n");
          exit(0);
        }else{
           if(cargarImagenPGM(&imagen,argv[2])==0){
              opcion = PDM;
                cargarPuntos(argv[4],puntosselec);
                Get_X_Y(&a,&b);
                float alpha_sn,beta_sn,teta_sn,w_term,w_line,w_edge,w_press;
                alpha_sn = atof(argv[5]);
                beta_sn = atof(argv[6]);
                teta_sn = atof(argv[7]);
                w_line = atof(argv[8]);
                w_edge = atof(argv[9]);
                w_term = atof(argv[10]);
                w_press = atof(argv[11]);
                inicializarSnake(imagen,puntosselec,a,b,
                alpha_sn,beta_sn,teta_sn,w_term,w_line,w_edge,w_press,false,0,0);
           }else{
             printf("Error al abrir archivo de imagen.\n");
             exit(0);
           }		 
        }
      break;
      case 'L' :
        if (argc != 11){
          printf("Numero de parametros incorrecto. El correcto es 10.\n");
          exit(0);
        }else{
          if(cargarImagenPGM(&imagen,argv[3])==0){
              opcion = LS;
              cargarPuntos(argv[5],puntosselec);
                Get_X_Y(&a,&b);
                float sigma,epsil,thr;
                sigma = atof(argv[6]);
                nb = atoi(argv[7]);
                increm = atof(argv[8]);
                epsil = atof(argv[9]);
                thr = atof(argv[10]);
                if(strcmp(argv[2],"-C")==0){
                  inicializar_level_set(imagen,increm,1,sigma,nb,epsil,puntosselec,a,b,thr,false,0,false,0,0,true);				  
                }else{
                  if(strcmp(argv[2],"-E")==0){
                    inicializar_level_set(imagen,increm,0,sigma,nb,epsil,puntosselec,a,b,thr,false,0,false,0,0,true);
                  }else{
                    printf("Defina compresion -C o extension -E.\n");
                    exit(0);
                  }
                }
				
              siguiente_LevelSet(&MatrizEnerg);            
           }else{
             printf("Error al abrir archivo de imagen.\n");
             exit(0);
           }
        }
       break;
       case 'G' :
        if (argc != 11){
          printf("Numero de parametros incorrecto. El correcto es 10.\n");
          exit(0);
        }else{
           if(cargarImagenPGM(&imagen,argv[3])==0){
              opcion = GEOD;
                cargarPuntos(argv[5],puntosselec);
                Get_X_Y(&a,&b);
                float sigma,epsil,thr;
                sigma = atof(argv[6]);
                nb = atoi(argv[7]);
                increm = atof(argv[8]);
                epsil = atof(argv[9]);
                thr = atof(argv[10]);
                if(strcmp(argv[2],"-C")==0){
                  
				   inicializar_level_set(imagen,increm,1,sigma,nb,epsil,puntosselec,a,b,thr,false,0,false,0,0,true);
                }else{
                  if(strcmp(argv[2],"-E")==0){
                     inicializar_level_set(imagen,increm,0,sigma,nb,epsil,puntosselec,a,b,thr,false,0,false,0,0,true);
                  }else{
                    printf("Defina compresion -C o extension -E.\n");
                    exit(0);
                  }
                }
				siguiente_Geodesica(&MatrizEnerg);
           }else{
             printf("Error al abrir archivo de imagen.\n");
             exit(0);
           }
        }
       break;
      }
      int iterac = 1;
      imag = (unsigned char*)calloc(sizeof(unsigned char),a*b);
      switch(opcion){
        case LS: case GEOD:
          while(!parada(MatrizEnerg)){
            tiempo1 = clock();
            if((iterac%nb==0) || fin_narrow()){
                reinicializa_ls();
            }
            if(opcion == LS)
              siguiente_LevelSet(&MatrizEnerg);
            else
              siguiente_Geodesica(&MatrizEnerg);
            tiempo2 = clock();
            tiempo_total += tiempo2-tiempo1;
            iterac++;
            printf(".");
          }
          printf("\nFin de la segmentacion.\n");
          printf("Numero de Iteraciones: %d\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
          filtrado_burb();
          actual_level_set(imag);
          Set_X_Y(a,b);
          if(opcion == LS)
            guardarImagenPGM(imag,argv[4]);
          else
            guardarImagenPGM(imag,argv[4]);
          liberar_ls();
		  free(imag);
		  free(imagen);
        break;
        case PDM:
         while(!parada_snk()||(iterac==1)){
           tiempo1 = clock();
           siguienteSnake();
           iterac++;
           printf("%d\n",iterac);
           tiempo2 = clock();
           tiempo_total += tiempo2-tiempo1;
           printf(".");
         }
         printf("\nFin de la segmentacion.\n");
         printf("Numero de Iteraciones: %d\nTiempo Computacional= %f\n",iterac,(double)tiempo_total/CLOCKS_PER_SEC);
         pintaguardaSnakeInterp(imag);
         guardarImagenPGM(imag,argv[3]);
         libera_snk();
		 free(imag);
		 free(imagen);
        break;
      }
  }
  l_dest(&puntosselec);
  return EXIT_SUCCESS;
}

           

  
  
  

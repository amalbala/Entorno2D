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



#ifndef _SNAKES_H_
#define _SNAKES_H_

 void inicializarSnake(unsigned char *imag,lista puntos,int x_imag, int y_imag,
 float alpha_s,float beta_s,float teta_s,float wt_s,float wl_s,float we_s,float wp_s,
 	int fuzzy,int num_c,int cl_int); 
 void siguienteSnake();
 void siguienteSnakeFuzzy();
 void pintaguardaSnake(unsigned char *imag);
 void pintaguardaSnakeInterp(unsigned char *imag);
 int parada_snk();
   
 void libera_snk();
 void guardaSnakeInterp();

#endif

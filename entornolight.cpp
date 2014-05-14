
//Librerias utilizadas

#include <GL/glut.h>
#include <GL/glui.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <fstream.h>

#include "imagutils.h"
#include "listagen.h"




#define PATH_IMAG ""


//Lista de constantes para los menus
enum {CARGAR,SUAVIZAR,BINAR,BINAR_UMB,GRADIEN,LAPLAC,NEGAT,GAUSS_GRAD,GUARDAR,
  ENSUCIAR,NUEVA_SELEC};


int extension;
int direccion;

//angulo de la camara
GLint angulo = 0;

//Dimensiones de la ventana de vision
GLsizei weight = 800;
GLsizei height = 800;
//Deplazamiento de centrado de la imagen
GLfloat desplaz_x;
GLfloat desplaz_y;
//Seleccion del menu
GLint opcion = 10000;

//Posiciones relativas del raton
GLfloat x_m = 0.0;
GLfloat y_m = 0.0;

int winIdMain;
GLUI* glui;
char label[100];
char name[sizeof(GLUI_String)];
char arch[40];
int x_imag,y_imag;
static unsigned char *imag;
unsigned char* imagen;

lista puntosselec;

float SIGMA_GAUSS = 0.6;

GLUI_Rollout *roFuncGenerales;
GLUI_Rollout *roFuncAvanz;
void visualizacion_avanzada();
void visualizacion_simple();

/*--------------------------------------------------*/
/* FUNCIONES DE DIBUJADO DE OPENGL					*/
/*--------------------------------------------------*/
/*--------Tralacion de pixeles ---------------------*/
void swappixelOGL(unsigned char *image){
  unsigned char *imag_tras;
  imag_tras = (unsigned char*)calloc(sizeof(unsigned char),x_imag*y_imag);

  int k=0;
  int i,j;

  for(j=y_imag; j>0; j--){
    for(i=0;i<x_imag; i++){
      imag_tras[k] = image[calcSingleSubscript(i,j-1,x_imag)];
      k++;
    }
  }

  for(i=0; i<x_imag*y_imag;i++)
    image[i] = imag_tras[i];

}

/*------------------- DRAW STRING -----------------------*/
//Dibuja una cadena en el lugar actual
void drawString (char *s)
{
  for(int i=0; i< (int)strlen(s); i++)
    glutBitmapCharacter (GLUT_BITMAP_HELVETICA_10, s[i]);

}
/*---------------- SUB MENU -------------------------------*/
void guardarPuntos(lista list,char arch[40]){
  ofstream desc;
  desc.open(arch);
  punto2D pto;
  l_empieza(list);
  while(l_quedan(list)){
    l_dame(list,&pto);
    desc<<pto.x<<" "<<pto.y<<"\n";
  }
  desc.close();
}
/*---------------- SUB MENU -------------------------------*/
void cargarPuntos(char arch[40],lista list){
  FILE *archivo;
  archivo = fopen(arch,"r");
  punto2D pto;
  while(fscanf(archivo,"%d %d\n",&pto.x,&pto.y)!= EOF){
    l_meted(list,&pto);
  }
  fclose(archivo);
}

/*------------------- CONTROL MENU -----------------------*/
//Para controlar el menu
//AQUI HAY QUE INSERTAR NUEVAS FUNCIONES.
void controlmenu(int value){//Casos del menu de opciones.
   int i;
	switch(value){
    case CARGAR:
      sprintf(arch,PATH_IMAG);
      strcat(arch,name);

      switch(extension){
        case 0:
          strcat(arch,".pgm");
          cargarImagenPGM(&imagen,arch);
        break;
        case 1:
          strcat(arch,".raw");
          cargarImagenRAW(&imagen,arch,512,512);
        break;
      }
      Get_X_Y(&x_imag,&y_imag);
      swappixelOGL(imagen);
      imag = (unsigned char *)calloc(sizeof(unsigned char),x_imag * y_imag);
      memcpy(imag,imagen,sizeof(unsigned char)*x_imag*y_imag);
      for(i=0; i<(int)strlen(name);i++)
        name[i] = '\0';
    break;
      case SUAVIZAR:
       double *imagen_suav;
        imagen_suav = (double*)calloc(sizeof(double),x_imag*y_imag);
       gaussiana(SIGMA_GAUSS,imag,imagen_suav);
       doubleTounsignedchar(imagen_suav,imag);
       free(imagen_suav);
      break;
      case BINAR:
      
      break;
      case BINAR_UMB:
      break;
      case NEGAT:
       negativo(imag,imag);
      break;
      case ENSUCIAR:
       ruido_aleatorio(imagen,imagen);
       memcpy(imag,imagen,(sizeof(unsigned char) * x_imag * y_imag));
      break;
	  case GAUSS_GRAD:
		  double *image_gradgauss;
		  float max,min;
		  image_gradgauss = (double*)calloc(sizeof(double),x_imag*y_imag);
		  grad_conv_gauss(imag,0.5,&max,&min,image_gradgauss );
		  doubleTounsignedchar(image_gradgauss,imag);
		  free(image_gradgauss);
	  break;
      case GRADIEN:
        double *image;
		    double *image_prw;
        image = (double*)calloc(sizeof(double),x_imag*y_imag);
		    image_prw = (double*)calloc(sizeof(double),x_imag*y_imag);
        for(i=0; i<(x_imag*y_imag); i++)
          image[i] = imag[i];
        grad_diffin(image,image_prw);
        doubleTounsignedchar(image_prw,imag);
        free(image);
      break;
      case LAPLAC:
        double *imag_dou;
		imag_dou = (double*)calloc(sizeof(double),x_imag*y_imag);
        laplaciano(imag,imag_dou);
        doubleTounsignedchar(imag_dou,imag);
        free(imag_dou);
      break;
      case NUEVA_SELEC:
        l_dest(&puntosselec);
        puntosselec = l_nuev(sizeof(punto2D));
        memcpy(imag,imagen,(sizeof(unsigned char) * x_imag * y_imag));
      break;
      case GUARDAR:
        swappixelOGL(imag);
        guardarImagenPGM(imag,"pantalla.pgm");
        swappixelOGL(imag);
      break;
	}
 glutPostRedisplay();
}


/*---------------- VENTANA PRINCIPAL -------------------------------*/
void mainDisplay (void)
{
  /* Clean drawing board */
  glutSetWindow (winIdMain);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  if(!l_vacia(puntosselec)){
    l_empieza(puntosselec);
    while(l_quedan(puntosselec)){
      punto2D pto;
      l_dame(puntosselec,&pto);
      imag[calcSingleSubscript(pto.x,pto.y,y_imag)] = 255.0;
    }
  }
    
  desplaz_x = -1.0 + (((weight - x_imag)/2.0)/(weight/2.0));
  desplaz_y = -1.0 + (((height - x_imag)/2.0)/(height/2.0));;
  glRasterPos2f(desplaz_x,desplaz_y);
  glDrawPixels(x_imag,y_imag,GL_LUMINANCE,GL_UNSIGNED_BYTE,imag);
  glutSwapBuffers();
}
/*---------------- INIT -------------------------------*/
//Inicializacion.
void init() {
	//Definicion de la luz
	GLfloat light_ambient[] = { 0.75, 0.75, 0.75, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.1, 0.25, 1.0, 0.0 };

	glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv (GL_LIGHT0, GL_POSITION, light_position);

	//Activamos las luces
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT0);

  //Tipo de modelos
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	//Color de fondo
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
//  glPixelStorei(GL_UNPACK_ALIGNMENT,2);
  puntosselec = l_nuev(sizeof(punto2D));
}

/*------------- RESHAPE -------------------------------------*/
void mainReshape (int w, int h)
{
  int tx,ty,tw,th;
  GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);

  glViewport (tx, ty, tw, th);
  weight = w;
  height = h;
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D (-1.0F, 1.0F, -1.0F, 1.0F);
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();

}

/*------------- idle -------------------------------------*/
void idle (void){
  /* Update main and sub window */
  glutSetWindow (winIdMain);
  glutPostRedisplay ();

}

/*------------- MOTION -------------------------------------*/

//Rutina de control del movimiento del raton
void motion(int x, int y)
 {
   punto2D pto;
   pto.x = x - (height - x_imag) + ((weight - x_imag)/2.0);
   pto.y = (y_imag - y) + (weight - y_imag) - ((height  - y_imag)/2.0);
   l_meted(puntosselec,&pto);
 }

/*--------------- MOUSE -----------------------------------*/

//Rutina de control del raton
void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		motion(x,y);
	}

}

/*--------------- MAIN ---------------------------------*/



int main(int argc,char** argv){
	for(int i=0; i<(int)strlen(name); i++)
    name[i] = '\0';
	glutInit(&argc,argv);
	//inicializa la ventana en el modo Doble Buffer
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	//Tamao de la ventana
	glutInitWindowSize (600, 600);
	//Posicion de la ventana
	glutInitWindowPosition (0, 0);
	//Crea la Ventana
	winIdMain = glutCreateWindow ("Visualizador Imagenes");
	//sentencias de Inicio
	init ();
	//Funcion de dibujo
	glutDisplayFunc(mainDisplay);
	//Funcion cuando cambiemos el tamao de la ventana
  //GLUI_Master.set_glutReshapeFunc(mainReshape);
  glutReshapeFunc(mainReshape);
	//Control del teclado
 // GLUI_Master.set_glutKeyboardFunc(keyboard);
  //Control del Raton
  GLUI_Master.set_glutMouseFunc(mouse);
	glutMotionFunc(motion);
  //Control de las ventanas
  GLUI_Master.set_glutIdleFunc(idle);
  //Creacion del menu de level set

  //Creacion de la subventana
  glui = GLUI_Master.create_glui("Controles...");
  glui->set_main_gfx_window(winIdMain);
 //Panel de Apertura
  GLUI_Panel *pnAbrir = glui->add_panel("Abir MRI");
  glui->add_statictext_to_panel(pnAbrir,"Archivo:");
  glui->add_edittext_to_panel(pnAbrir,"",GLUI_EDITTEXT_TEXT,
    &name);
  GLUI_RadioGroup *rgExtension = glui->add_radiogroup_to_panel(pnAbrir,&extension);
  glui->add_radiobutton_to_group(rgExtension,"PGM");
  glui->add_radiobutton_to_group(rgExtension,"RAW");
  glui->add_separator_to_panel(pnAbrir);
  glui->add_button_to_panel(pnAbrir,"Abrir",CARGAR,(GLUI_Update_CB)controlmenu);
  glui->add_button("Guardar Pantalla",GUARDAR,controlmenu);
  glui->add_button("Salir",0,exit);
  glui->add_column(true);
  //Funciones Generales
  roFuncGenerales = glui->add_rollout("Funciones Generales",true);
  visualizacion_simple();
  glui->add_column(true);
  //Panel de Gaussianas
  roFuncAvanz = glui->add_rollout("Funciones Generales",true);
  visualizacion_avanzada();
  //inicia el proceso de eventos
	glutMainLoop();
	return 0;
}


void add_func_simple(char* text, int constante){
  glui->add_button_to_panel(roFuncGenerales,text,
     constante,(GLUI_Update_CB)controlmenu);
}



void visualizacion_simple(){
  glui->add_button_to_panel(roFuncGenerales,"Nueva Selec.",
     NUEVA_SELEC,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Binarizar OTSU",
     BINAR,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Binarizar 125",
    BINAR_UMB,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Grad Prewitt",
    GRADIEN,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Laplaciano",
    LAPLAC,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Negativo",
    NEGAT,(GLUI_Update_CB)controlmenu);
  glui->add_button_to_panel(roFuncGenerales,"Ensuciar",
    ENSUCIAR,(GLUI_Update_CB)controlmenu);    
}

void visualizacion_avanzada(){
  glui->add_button_to_panel(roFuncAvanz,"Gausiana",
    SUAVIZAR,(GLUI_Update_CB)controlmenu);
  GLUI_Spinner *spGauss = glui->add_spinner_to_panel(roFuncAvanz,"Sigma: ",
    GLUI_SPINNER_FLOAT,&SIGMA_GAUSS);
  spGauss->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spGauss->set_speed(1.0);
  glui->add_separator_to_panel(roFuncAvanz);
  glui->add_button_to_panel(roFuncAvanz,"GradGauss",
    GAUSS_GRAD,(GLUI_Update_CB)controlmenu);
  GLUI_Spinner *spGradGauss = glui->add_spinner_to_panel(roFuncAvanz,"Sigma: ",
    GLUI_SPINNER_FLOAT,&SIGMA_GAUSS);
  spGradGauss->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  spGradGauss->set_speed(1.0);
}
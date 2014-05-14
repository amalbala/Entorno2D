#include "Imagen2D.h"
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//Constructor por defecto
Imagen2D::Imagen2D(){}


//Constructor pasandole dimensiones
Imagen2D::Imagen2D(int NumRow, int NumCol){
	this->nc = NumCol;
	this->nr = NumRow;
	this->data = (unsigned char*)calloc(sizeof(unsigned char),nr*nc);
}

//Constructor desde archivo
Imagen2D::Imagen2D(char *sNombre, char* sExtension){
	if((strcmp(sExtension,"pgm") == 0)||(strcmp(sExtension,"PGM") == 0)){
		char sNombreCompleto[40];
		sNombreCompleto[0]='\0';
		strcat(sNombreCompleto,sNombre);
		strcat(sNombreCompleto,".");
		strcat(sNombreCompleto,sExtension);
		this->ReadImagePGM(sNombreCompleto);
	}
}

//Leer de una imagen PGM
int Imagen2D::ReadImagePGM(char *arch){
  ifstream archivo;
  char cabecera[100];   
  archivo.open(arch,fstream::binary);   
  int i;
  if(!archivo) return 1;
  //Quito la cabecera
  for(i=0; i<2; i++){
	archivo.getline(cabecera,100);
  }
	if(archivo.bad() || archivo.fail()) return 2;
	//Almaceno las dimensiones de la imagen
	archivo>>nc>>nr;
	//descarto la cabecera restante
	archivo.getline(cabecera,100);
	if(archivo.bad() || archivo.fail()) return 2;
	archivo.getline(cabecera,100);
	if(archivo.bad() || archivo.fail()) return 2;
	//Inicializo el vector que contendra la imagen
	
	//Cargo la imagen
	archivo.read((char*)this->data,(nc*nr));
	if(archivo.bad() || archivo.fail()) return 2;
  //cierro el archivo
	archivo.close();
	return 0;
}

//Escribir en Fichero
Imagen2D::WriteImage(char *sNombre, char *sExtension){
	ofstream desc;
	desc.open(nombre,fstream::binary);	
	desc<<"P5\n"<<"# Created by Antonio Martinez\n"<<nc<<" "<<nr<<"\n"<<"255"<<"\n";
	desc.write((char*)(this->data),nc*nr);
	desc.close();
}


//Destructor
Imagen2D::~Imagen2D(){
	free(this->data);
}
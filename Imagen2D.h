


//Clase Imagen2D
class Imagen2D {

//Miembros protegidos, solo accesible por la clase y sus derivados
protected:

	int nc;
	int nr;
	int color;
	
	unsigned char *data;

//Miembros privados, solo accesible por la clase
private:
	ReadImagePGM(char *arch);

//Miembros publicos
public:
	Imagen2D();
	~Imagen2D();
	Imagen2D(int nr, int nc);
	Imagen2D(char *sNombre, char* sExtension);

	WriteImage(char *sNombre, char *sExtension);

};
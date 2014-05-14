MAIN=Entorno2D
CC=gcc
CFLAGS= -Wall
INCLUDES=
LFLAGS=
SRCS =  error.cpp listagen.cpp imagutils.cpp listalista.cpp snakes.cpp GARModel.cpp levelset.cpp ImageVisor.cpp 

CFLAGS += 
LFLAGS += -L. -L/usr/X11R6/lib -lm -lGL -lGLU -lglut -lglui -lpthread 

OBJS = $(SRCS:.cpp=.o)

$(MAIN):  $(OBJS)
	$(CC) $(OBJS) -o $(MAIN) $(LFLAGS)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *.o *~ $(MAIN)



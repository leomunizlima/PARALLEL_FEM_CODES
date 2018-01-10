CC=gcc
IDIR= 
CFLAGS= -lm -lrt -c -Wall -Ofast -march=native
GFLAGS= -lm -lrt -c -Wall -Ofast -march=native
LDFLAGS= 
SOURCES= main.c \
reordering_rcm.c \
reading_mesh.c \
writing_mesh.c \
checking_mesh.c \
mycalloc.c \
myfopen.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=Mesh_Reordering

all: $(SOURCES) $(EXECUTABLE) 
all: 
	gfortran -o triangulation_rcm triangulation_rcm.f90 -O2 -march=native

debug: CFLAGS= -lm -lrt -c -g -Wall
debug: GCFLAGS= 
debug: all

 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) triangulation_rcm






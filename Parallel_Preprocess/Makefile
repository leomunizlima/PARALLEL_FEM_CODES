CC         = gcc
CFLAGS     = -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS    = ./hsl_mc73/hsl_mc73d.o ./hsl_mc73/libhsl_mc73.a ./hsl_mc64/mc64d.o ./hsl_mc64/libmc64.a -lrt -lm -lgfortran -lblas -llapack -Wall -Ofast -march=native
#GCFLAGS    = -lrt -lm -Wall -Ofast -march=native
SOURCES    = Preprocess.c\
partitioning.c\
symrcm.c \
spectral.c
		
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = Preprocess 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -c -lm -lrt -Wall 
debug: GCFLAGS=./hsl_mc73/hsl_mc73d.o ./hsl_mc73/libhsl_mc73.a ./hsl_mc64/mc64d.o ./hsl_mc64/libmc64.a -lrt -lm -lgfortran -lblas -llapack -Wall 
#debug: GCFLAGS= -lrt -lm -Wall -g
debug: all

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clear_ranks:
	gcc -o clear_ranks clear_ranks.c



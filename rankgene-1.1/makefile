CC      = gcc
CXX     = g++
#CFLAGS  = -O1
CFLAGS  = -g
LIBS    = -lm

SRC     = util.c load_data.c load_genes.C impurity_measures.c rankgene.C
OBJ     = util.o load_data.o load_genes.o impurity_measures.o rankgene.o
INC     = oc1.h load_genes.h
MAKEFILE= makefile

default: rankgene

clean:
	rm -rf $(OBJ) rankgene

.c.o: $(INC) $(MAKEFILE)
	$(CC) $(CFLAGS) -c $*.c

.C.o: $(INC) $(MAKEFILE)
	$(CXX) $(CFLAGS) -c $*.C

rankgene: $(MAKEFILE) $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) -o rankgene $(LIBS)

$(OBJ):$(INC) $(MAKEFILE)


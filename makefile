CC=gcc
CFLAGS+= -m64 -g -Wall
LDFLAGS= -L$$GUROBI_HOME/lib -lgurobi91
LIBS=
INC= $$GUROBI_HOME/include/
ALLDEP:= util.o ArrayList.o HashTable.o #AVLTree.o

.PHONY: test

test: $(ALLDEP) test.out partitionByLayers.out

product: CFLAGS = -O3
product: $(ALLDEP) partitionByLayers.out partitionByLayersCheckByCenters.out partitionByLayersCheckByCentersWithP1M1.out partitionByLayersCheckByNeighborsWithP1M1.out partitionByLayersCheckByNeighbors.out partitionByLayersCheckByNeighborsOnlyP1M1.out assignAll.out

%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.c makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

%.MaxID: %.txt
	grep -oE '[0-9]+' $^ | tail -1 > $@; \
	sed -n '7 s/\r$$//p' $^ | sed 's/ /\n/g' >> $@

#gurobi make
%_ILP: %_ILP.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ -I$(INC) $(LDFLAGS) -lm

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

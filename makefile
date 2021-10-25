CC=gcc
CFLAGS+= -g -Wall
LDFLAGS=
LIBS=
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

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

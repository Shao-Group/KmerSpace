CC=gcc
CFLAGS+= -g -Wall
LDFLAGS=
LIBS=
ALLDEP:= util.o ArrayList.o #AVLTree.o

.PHONY: test

test: $(ALLDEP) test.out partitionByLayers.out

product: CFLAGS = -O3
product: $(ALLDEP) partitionByLayers.out

%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.c makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

%.MaxID: %.txt makefile
	grep -oE '[0-9]+' $< | tail -1 > $@; \
	#sed -n '7 s/\r$$//p' $< | sed 's/ /\n/g' >> $@
	sed -n '5p' $< | sed 's/ /\n/g' >> $@

%.MaxID.win: %.txt makefile
	grep -oE '[0-9]+' $< | tail -1 > $@; \
	sed -n '7 s/\r$$//p' $< | sed 's/ /\n/g' >> $@; \
	mv $@ $(basename $@)

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

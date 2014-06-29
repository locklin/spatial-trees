#
# ranger Makefile mark 2
#
CC = cc

# Where to find libraries and include files
#  export LD_LIBRARY_PATH=~/libranger.so.1:$LD_LIBRARY_PATH
COMPILER = gcc
LIBS =  -lm
CFLAGS = #-I
COPTS = -fPIC -g -c -Wall -o
LIBPATH = 

ranger.so: distrib.o pqueue.o eigens.o metric.o naive.o naivekd.o optkd.o \
	vptree.o sproullkd.o
	$(COMPILER) -shared -Wl,-soname,libranger.so.1 -o libranger.so.1 \
	distrib.o pqueue.o eigens.o metric.o naive.o naivekd.o optkd.o \
	vptree.o sproullkd.o -lc	

clean:
	rm *.o *.so.1

dist:
	tar -cf ../spatial-trees.tar *.c *.h Makefile README LICENSE sample_data

fileoper.o: fileoper.c
	$(COMPILER) $(COPTS) fileoper.o fileoper.c  $(LIBS) $(CFLAGS)

eigens.o: eigens.c 
	$(COMPILER) $(COPTS) eigens.o eigens.c  $(LIBS) $(CFLAGS)

distrib.o: distrib.c
	$(COMPILER) $(COPTS) distrib.o distrib.c -lm  $(LIBS) $(CFLAGS)

pqueue.o: pqueue.c 
	$(COMPILER) $(COPTS) pqueue.o pqueue.c  $(LIBS) $(CFLAGS)

metric.o: metric.c
	$(COMPILER) $(COPTS) metric.o metric.c  $(LIBS) $(CFLAGS)

naive.o: naive.c
	$(COMPILER) $(COPTS) naive.o naive.c  $(LIBS) $(CFLAGS)

naivekd.o: naivekd.c naivekd.h
	$(COMPILER) $(COPTS) naivekd.o naivekd.c  $(LIBS) $(CFLAGS)

optkd.o: optkd.c optkd.h 
	$(COMPILER) $(COPTS) optkd.o optkd.c  $(LIBS) $(CFLAGS)

vptree.o: vptree.c vptree.h
	$(COMPILER) $(COPTS) vptree.o vptree.c  $(LIBS) $(CFLAGS)

sproullkd.o: sproullkd.c sproullkd.h eigens.o
	$(COMPILER) $(COPTS) sproullkd.o sproullkd.c  $(LIBS) $(CFLAGS)


## utils rosseler makes rosseler array, create is an eventual test script
rosseler: rosseler.c
	$(COMPILER) $(COPTS) rosseler rosseler.c $(LIBS) $(CPFLAGS)

create: create.c
	$(COMPILER) $(COPTS) create create.c $(LIBS) $(CPFLAGS)

test.o: test.c fileoper.o
	$(COMPILER) $(COPTS) test.o test.c $(LIBS) $(CPFLAGS) 

test: test.o fileoper.o ranger.so
	$(COMPILER) test.o fileoper.o libranger.so.1 -lm -L./libranger.so.1 -o test 


#voronoi/VorMain.o: voronoi/VorMain.c voronoi/edgelist.c voronoi/geometry.c voronoi/heap.c voronoi/memory.c voronoi/output.c voronoi/voronoi.c
#	$(COMPILER) $(COPTS) voronoi/VorMain.o voronoi/VorMain.c


# I don't know if this is good for anything yet
queryparser.o: queryparser.c # I think this is part of the visual bits only -SCL
	$(COMPILER) $(COPTS) queryparser.o queryparser.c  $(LIBS) $(CFLAGS) 


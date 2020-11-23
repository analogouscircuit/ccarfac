CC = gcc
CFLAGS = -g -Wall -O3 
LIBS = -lm -lsndfile

clib: library
test: carfac_test
default: library

carfac_test: carfac_test_main.o siggen.o carfac.o ctopy.o
	$(CC) -o carfac_test carfac_test_main.o siggen.o carfac.o ctopy.o $(LIBS)

carfac_test_main.o: carfac_test_main.c siggen.h ctopy.h
	$(CC) $(CFLAGS) $(LIBS) -c carfac_test_main.c 

ctopy.o: ctopy.c ctopy.h
	$(CC) $(CFLAGS) $(LIBS) -c ctopy.c

carfac.o: carfac.c carfac.h
	$(CC) $(CFLAGS) $(LIBS) -c carfac.c 

siggen.o: siggen.c siggen.h
	$(CC) $(CFLAGS) $(LIBS) -c siggen.c 

library: carfac.c carfac.h
	$(CC) -shared -fPIC carfac.c -o libcarfac.so

clean:
	$(RM) car*.o *~
	rm -r build pycarfac.c pycarfac*.so


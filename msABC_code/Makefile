CC=gcc
OBJECTS=ms.o data_sumstat.o rand1.o streec.o
CCFLAGS=-c -g -Wall
LIBS=-lm

msABC: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o msABC -g -lm
ms.o: ms.c
	$(CC) $(CCFLAGS) ms.c
data_sumstat.o: data_sumstat.c
	$(CC) $(CCFLAGS) data_sumstat.c
rand1.o: rand1.c
	$(CC) $(CCFLAGS) rand1.c
streec.o: streec.c
	$(CC) $(CCFLAGS) streec.c
clean:
	rm -f *.o
	rm -f msABC

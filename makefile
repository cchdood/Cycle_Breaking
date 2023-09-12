CC = g++
CFLAGS = -c

all	: bin/cb

bin/cb         : main.o
			$(CC) -o bin/cb main.o
main.o          : src/main.cpp
			$(CC) $(CFLAGS) src/main.cpp
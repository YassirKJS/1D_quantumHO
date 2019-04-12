CC= g++
CCTEST=cxxtestgen --error-printer
CFLAGS= -std=c++11 -Wall
TARGET= main
OBJS= main.cpp Calcul.cpp
BINDIR= ../bin/

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -larmadillo
	./*.out
	

clean:
	rm -f *.o *.out

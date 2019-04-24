CC = g++
CFLAGS = -std=c++11 -Wall
TESTGEN = cxxtestgen
TESTFLAGS = --error-printer
TESTS = CalculTest.h DataInterfaceTest.h Orthonormality.h
TESTSCPP = DataInterface.cpp Calcul.cpp Miscellaneous.cpp Orthonormality.cpp
OBJSTEST = DataInterface.o Calcul.o Miscellaneous.o Orthonormality.o
TARGET = main
TARGETTEST = test
OBJS = main.cpp Calcul.cpp DataInterface.cpp Miscellaneous.cpp Orthonormality.cpp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -larmadillo 
	./$(TARGET)

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< -c -o $@

$(TARGETTEST): $(TESTS)
	$(CC) $(CFLAGS) -c $(TESTSCPP) -larmadillo -llapack -lblas
	$(TESTGEN) $(TESTFLAGS) -o test.cpp $(TESTS)
	$(CC) $(CFLAGS) -o $@  test.cpp $(OBJSTEST) -larmadillo -llapack -lblas
	./$@

.PHONY: test clean

clean:
	rm -f $(TARGET) $(TARGETTEST) *.o t.cpp test.cpp *.orig

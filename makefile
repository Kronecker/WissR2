CFLAGS=
COMPILER=g++
all: B5A10

clean:
	rm *.out
B5A10: B5A10.cpp
	$(COMPILER) B5A10.cpp $(CFLAGS) -o B5_altSchwarz.out
CFLAGS=
COMPILER=g++
all: B4A8

clean:
	rm *.out
B4A8: 	B4A8.cpp
	$(COMPILER) B4A8.cpp $(CFLAGS) -o B4.out
B4A8Debug: 	B4A8.cpp
	$(COMPILER) B4A8.cpp -ggdb -o B4.out
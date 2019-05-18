numProc=4

all: B2A1 B2A2

clean:
	rm *.out
B2A1: B2A1.c
	mpicc B2A1.c -o B2A1.out
B2A1run: B2A1
	mpirun -n $(numProc) B2A1.out
B2A2: B2A2.c
	mpicc B2A2.c -oB2A2.out
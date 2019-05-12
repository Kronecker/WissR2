numProc=4

B2A1: B2A1.c
	mpicc B2A1.c -o B2A1.out
B2A1run: B2A1
	mpirun $(numProc) B2A1.out
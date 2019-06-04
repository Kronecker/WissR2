CFLAGS=-O3
all: B3A7Seriell B3A7Parallel B3A7FullAsync B3A7scripted

clean:
	rm *.out
B3A7Seriell:
	mpic++ B3A7.cpp -o B3seriell.out -DINNER_GRID_SIZE=1024 -DSAVE_MATRIX
B3A7Parallel:
	mpic++ B3A7.cpp -o B3parallel.out -DINNER_GRID_SIZE=1024 -DSAVE_MATRIX -DMPI_PARALLEL
B3A7FullAsync:
	mpic++ B3A7.cpp -o B3parallelAsync.out -DINNER_GRID_SIZE=1024 -DSAVE_MATRIX -DMPI_PARALLEL -DGO_FULL_ASYNC
B3A7scripted:
	mpic++ B3A7.cpp -o B3scripted.out -DINNER_GRID_SIZE=1024 -DMPI_PARALLEL -DSHOW_PROC_LAYOUT=0 -DGO_FULL_ASYNC

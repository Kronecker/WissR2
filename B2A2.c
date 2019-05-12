//
// Created by Looky on 12/05/2019.
//
#include <mpi.h>
#include <stdio.h>
#include <stlib.h>



#define N 500000


double* allocInitVecAUp(int startIndex, int numElements);
double* allocInitVecBDwn(int startIndex, int numElements);








double* allocInitVecAUp(int startIndex, int numElements) {
    double* vec;
    vec = (double*) malloc(sizeof(double) * numElements)
    for (int k=0; k<numElements; k++;) {
        vec[k]=k+1+startIndex;
    }
    return vec;
}


double* allocInitVecBDwn(int startIndex, int numElements) {
    double* vec;
    vec = (double*) malloc(sizeof(double) * numElements)
    for (int k=0; k<numElements; k++;) {
        vec[k]=N-k-startIndex;
    }
    return vec;
}



int main (int argc,char *argv[]) {
    printf("Hello world");
    int rank, size;
    MPI_INIT(&argc, char &argv[]);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *vecA, *vecB;

    vecA=allocInitVecAUp(10,2);
    vecB=allocInitVecBDwn(10,2);
    // Compute Scalprod


    // Communicate with master and send results
    if(rank==0) {
        // rcv

    }  else  {
        // send

    }





    free(vecA);
    free(vecB);


    MPI_Finalize();
    printf("Hello world");
}
































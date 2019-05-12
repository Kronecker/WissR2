//
// Created by Looky on 12/05/2019.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>



#define N 500000


double* allocInitVecAUp(int startIndex, int numElements);
double* allocInitVecBDwn(int startIndex, int numElements);








double* allocInitVecAUp(int startIndex, int numElements) {
    double* vec;
    vec = (double*) malloc(sizeof(double) * numElements);
    for (int k=0; k<numElements; k++) {
        vec[k]=k+1+startIndex;
    }
    return vec;
}


double* allocInitVecBDwn(int startIndex, int numElements) {
    double* vec;
    vec = (double*) malloc(sizeof(double) * numElements);
    for (int k=0; k<numElements; k++) {
        vec[k]=N-k-startIndex;
    }
    return vec;
}



int main (int argc,char *argv[]) {


    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *vecA, *vecB, sum=0;
    double end, start;
    int numELements;


    if(rank==0) {
      start  = MPI_Wtime();
    }

    numElements=N/size;
    // Allocate Vetors (Malloc) and init with series a_i=i+1  b_i=N-i
    vecA=allocInitVecAUp(rank*numElements,numElements);
    vecB=allocInitVecBDwn(rank*numElements,numELements);

    // Compute Scalprod

    for(k=0;k<numElements;k++) {
        printf("%f %f", vecA[k], vecB[k]);
        sum+=a[k]*b[k];
    }

    printf("I'm %d : ",rank,sum);



    // Communicate with master and send results
    if(rank==0) {
        // rcv

    }  else  {
        // send

    }





    free(vecA);
    free(vecB);


    if(rank==0)  {
        end = MPI_Wtime();
        printf("%.6f", end-start);
    }


    MPI_Finalize();


}





























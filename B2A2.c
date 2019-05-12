//
// Created by Looky on 12/05/2019.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


#define USE_MPI_REDUCE 1

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
    int numElements;


    if(rank==0) {
      start  = MPI_Wtime();
    }

    numElements=N/size;

    // Allocate Vectors (Malloc) and init with series a_i=i+1  b_i=N-i
    vecA=allocInitVecAUp(rank*numElements,numElements);
    vecB=allocInitVecBDwn(rank*numElements,numElements);

    // Compute Scalprod
    for(int k=0;k<numElements;k++) {
        sum+=vecA[k]*vecB[k];
    }



    MPI_Status status;
    double partsum=sum;

#if USE_MPI_REDUCE

    MPI_Reduce(&sum,&partsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    sum=partsum;
#else


    // Communicate with master and send results
    if(rank==0) {
        printf("PartSum %d : %f\n",0,sum);
        // rcv
        for(int k=1;k<size;k++) {
            MPI_Recv(&partsum, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
            sum+=partsum;
            printf("PartSum %d : %f\n",k,partsum);
        }
    }  else  {
        // send
        MPI_Send(&sum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }

#endif

    printf("GlobalSum : %f\n",sum);


    free(vecA);
    free(vecB);


    if(rank==0)  {
        end = MPI_Wtime();
        printf("%d nodes in %.6fs\n", size ,end-start);
    }


    MPI_Finalize();


}





























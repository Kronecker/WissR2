//
// Created by Looky on 12/05/2019.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef USE_MPI_REDUCE
    #define USE_MPI_REDUCE 1
#endif

#ifndef LAST_VECTOR_TAKES_REST
    #define LAST_VECTOR_TAKES_REST 1
#endif

#ifndef N
    #define N 500000
#endif

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

#if LAST_VECTOR_TAKES_REST
    if(rank==(size-1)) {
        vecA = allocInitVecAUp(rank * numElements, N-rank*numElements);
        vecB = allocInitVecBDwn(rank * numElements, N-rank*numElements);
        numElements=N-rank*numElements;
    } else {
        vecA = allocInitVecAUp(rank * numElements, numElements);
        vecB = allocInitVecBDwn(rank * numElements, numElements);
    }
#else  
    int restElements=N%size;
    if((rank)<restElements) { // take one for the team
        vecA = allocInitVecAUp(rank * numElements + rank, numElements+1);
        vecB = allocInitVecBDwn(rank * numElements + rank, numElements+1);
        numElements++;
    }  else  {               // get carried
        vecA = allocInitVecAUp(rank * numElements + restElements, numElements);
        vecB = allocInitVecBDwn(rank * numElements + restElements, numElements);
    }
#endif


    // Compute Scalprod
    for(int k=0;k<numElements;k++) {
        sum+=vecA[k]*vecB[k];
    }



    MPI_Status status;

    #if USE_MPI_REDUCE  // -------------------------------------------------------------

        double partsum=sum;
        // Use Reduce to send all results to Master and Sum all values
        MPI_Reduce(&partsum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    #else  // SEND & RECEIVE -----------------------------------------------------------

        double partsum;
        // Communicate with master and send results
        if(rank==0) {
            // rcv
            for(int k=1;k<size;k++) {
                MPI_Recv(&partsum, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
                sum+=partsum;
            }
        }  else  {
            // send
            MPI_Send(&sum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }

    #endif  // -------------------------------------------------------------------------




    free(vecA);
    free(vecB);


    if(rank==0)  {
        printf("GlobalSum : %f\n",sum);
        end = MPI_Wtime();

        #if USE_MPI_REDUCE
            printf("%d nodes with Reduce in %.6fs\n", size ,end-start);
        #else
            printf("%d nodes with Snd/Rcv in %.6fs\n", size ,end-start);
        #endif


    }


    MPI_Finalize();


}





























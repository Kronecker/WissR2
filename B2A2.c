//
// Created by Looky on 12/05/2019.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef USE_MPI_REDUCE
    #define USE_MPI_REDUCE 0
#endif

#ifndef LAST_VECTOR_TAKES_REST
    #define LAST_VECTOR_TAKES_REST 0
#endif

#ifndef SND_RCV_WITH_TREE  // only active if USE_MPI_REDUCE is 0
    #define SND_RCV_WITH_TREE 0
#endif

#ifndef SND_RCV_TREE_NODES   // only active if  SND_RCV_WITH_TREE is active
    #define SND_RCV_TREE_NODES 2
#endif


#ifndef N
    #define N 500000
#endif

double* allocInitVecAUp(int startIndex, int numElements);
double* allocInitVecBDwn(int startIndex, int numElements);




int main (int argc,char *argv[]) {


    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *vecA, *vecB, sum=0;
    double end, start; // MPI Timestamps
    int numElements;


    if(rank==0) {  // Master measures time
      start  = MPI_Wtime();
    }

    numElements=N/size;

    // Allocate Vectors (Malloc) and init with series a_i=i+1  b_i=N-i

#if LAST_VECTOR_TAKES_REST  // Push remaining elements to the last mpi node
    if(rank==(size-1)) {
        vecA = allocInitVecAUp(rank * numElements, N-rank*numElements);
        vecB = allocInitVecBDwn(rank * numElements, N-rank*numElements);
        numElements=N-rank*numElements;
    } else {
        vecA = allocInitVecAUp(rank * numElements, numElements);
        vecB = allocInitVecBDwn(rank * numElements, numElements);
    }
#else  // Evenly Distributed Vecs
    int restElements=N%size;
    if((rank)<restElements) { // take one for the team
        vecA = allocInitVecAUp(rank * numElements + rank, numElements+1);
        vecB = allocInitVecBDwn(rank * numElements + rank, numElements+1);
        numElements++;
    }  else  {               // get away with less
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
        // Use MPI_Reduce to send all results to Master and Sum all values
        MPI_Reduce(&partsum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    #else  // SEND & RECEIVE -----------------------------------------------------------

        #if SND_RCV_WITH_TREE // SndRcv with tree, SND_RCV_TREE_NODES determines number of nodes per step
            int sendingNodes=SND_RCV_TREE_NODES-1;

            int moduloActive=0;
            int moduloSndRCV;
            int groupSpan=SND_RCV_TREE_NODES;
            int innerGroupHops=1;
            double partsum;

            while(size>innerGroupHops) {
                if(!moduloActive) {       // ACTIVE
                    moduloSndRCV=rank%groupSpan;

                    if(!moduloSndRCV) {   // RECEIVE
                        if(size-rank < groupSpan) {
                            sendingNodes=(size-rank-1)/innerGroupHops;
                        }
                        for(int k=1;k<=sendingNodes;k++) {
                            MPI_Recv(&partsum, 1, MPI_DOUBLE, rank+k*innerGroupHops, 0, MPI_COMM_WORLD, &status);
                            sum+=partsum;
                        }
                        sendingNodes=SND_RCV_TREE_NODES-1;
                    } else {              // SEND
                        MPI_Send(&sum,1,MPI_DOUBLE,rank-moduloSndRCV,0,MPI_COMM_WORLD);
                    }
                } else {                // INACTIVE
                    break;
                }
                innerGroupHops=groupSpan;
                groupSpan*=SND_RCV_TREE_NODES;
                moduloActive=moduloSndRCV;
            }

        #else // All slaves send to Master 0


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


        #endif


    #endif  // -------------------------------------------------------------------------




    free(vecA);
    free(vecB);


    if(rank==0)  {
        end = MPI_Wtime();

        #if USE_MPI_REDUCE
                printf("%3d MPI nodes in %.1fsm using MPI Reduce", size ,(end-start)*1000);
        #else
            #if SND_RCV_WITH_TREE
                printf("%3d MPI nodes in %.1fms using Snd/Rcv Tree with %d nodes", size ,(end-start)*1000, SND_RCV_TREE_NODES);
            #else
                printf("%3d MPI nodes in %.1fms using Snd/Rcv Slave->Master", size ,(end-start)*1000);
            #endif

        #endif
        printf("   GlobalSum : %e\n",sum);
    }
    MPI_Finalize();

}


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























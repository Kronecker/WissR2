#include <stdio.h>
#include <mpi.h>

int main(int argc,char *argv[]){

    int rank,size;
    MPI_Init(&argc,&argv); // MPI starten
    MPI_Comm_size(MPI_COMM_WORLD,&size); // Anzahl der Prozesse ermitteln
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); // Nummer des aktuellen Prozesse

    if(rank==0){

        printf("Hello World from Master of %d Processes.\n",size);

    }else{

        printf("Hello World from Slave %d.\n",rank);
    }
    MPI_Finalize(); // MPI beenden
}
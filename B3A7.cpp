//
// Created by Looky on 02/06/2019.
//

//<editor-fold desc="Header">

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>


#ifndef INNER_GRID_SIZE
#define INNER_GRID_SIZE 1024
#endif
#ifndef SHOW_MATRIX
#define SHOW_MATRIX 0
#endif
#ifndef SAVE_MATRIX
#define SAVE_MATRIX 0
#endif
#ifndef MPI_PARALLEL
#define MPI_PARALLEL 0
#endif



// Seriell
int main_jacobiSeriell(int argc, char *argv[]);

double *initMatrixRightHandSideSeriell(int n, double h);

void jacobiIterationWithExtGridSeriell(int n, double *f, double *gridCurrentIteration,
                                       double *gridLastIteration, int numberOfIterations, double h);

double *prepareGridAndExtendWithBoundaryVals(int n, double boundary);



// MPI Parallel
int main_jacobiMPI(int argc, char *argv[]);



// Utility
void displayMyMatrix(double *matrix, int m, int n);

void saveMyMatrix(double *matrix, int m, int n, double h, int numberTask);


//</editor-fold>

int main(int argc, char *argv[]) {

#if MPI_PARALLEL
    main_jacobiMPI(argc, argv);
#else
    main_jacobiSeriell(argc, argv);
#endif

}
//<editor-fold desc="Jacobi Seriell">
int main_jacobiSeriell(int argc, char *argv[]) {

    int n = INNER_GRID_SIZE + 2;
    printf("%d\n", n);
    double *currentIteration, *lastIteration, *f;
    double h;

    h = 1. / (INNER_GRID_SIZE + 1.);

    currentIteration = prepareGridAndExtendWithBoundaryVals(n, 0);
    lastIteration = prepareGridAndExtendWithBoundaryVals(n, 0);

    f = initMatrixRightHandSideSeriell(n, h);

        jacobiIterationWithExtGridSeriell(n, f, currentIteration, lastIteration, 500, h);



#if SHOW_MATRIX
    displayMyMatrix(lastIteration,n,n);
#endif

#if SAVE_MATRIX
    saveMyMatrix(lastIteration,n,n,h,42);
#endif

    delete (f);
    delete (currentIteration);
    delete (lastIteration);
}

double *prepareGridAndExtendWithBoundaryVals(int n, double boundary) {

    double *matrix = new double[n * n]();

    // boundary values init (outer)
    for (int i = 0; i < n; i++) {
        matrix[i] = boundary;
        matrix[n * (n - 1) + i] = boundary;
    }
    for (int k = 1; k < n - 1; k++) { // iterate through blocks
        matrix[k * n] = boundary;
        matrix[(k + 1) * n - 1] = boundary;
    }

    return matrix;
}

double *initMatrixRightHandSideSeriell(int n,
                                       double h) {   // parts of the matrix correspond to boundary values and are irrelevant
    double *matrix = new double[n * n];
    double x;
    double y;
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            x = h * i;
            y = h * j;
            matrix[i * n + j] = x * (1 - x) + y * (1 - y);
            //         printf("<%f %f> %f\n",x,y,matrix[i*m+j]);
        }
    }
    return matrix;
}

void jacobiIterationWithExtGridSeriell(int n, double *f, double *gridCurrentIteration, double *gridLastIteration,
                                       int numberOfIterations, double h) {
    double *temp;

    // some caching
    double hSquare = h * h;
    double oneBy4 = 1 / 4.;
    int nm1 = n - 1;
    int i ,k, index;

    for (int l = 0; l < numberOfIterations; l++) {

        for (k = 1; k < nm1; k++) {
            for (i = 1; i < nm1; i++) {
                index = k * n + i;
                gridCurrentIteration[index] = oneBy4 * (hSquare * f[index]
                                                        + (gridLastIteration[index - n] +
                                                           gridLastIteration[index - 1] +
                                                           gridLastIteration[index + 1] +
                                                           gridLastIteration[index + n]));
            }
        }
        temp=gridCurrentIteration;
        gridCurrentIteration=gridLastIteration;
        gridLastIteration=temp;
    }


}
//</editor-fold>


int main_jacobiMPI(int argc, char *argv[]) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);




    int n = INNER_GRID_SIZE + 2;
    double *currentIteration, *lastIteration, *f;
    double h;

    h = 1. / (INNER_GRID_SIZE + 1.);

    // reverse trial division, find "most square" configuration
    // under the assumption that all procs should have to work, not necessary the best configuration
    // e.g. 17 scales very bad, [17 1], but with one proc idle, 16 scales perfect (square, [4 4])

    int divi=sqrt(size);
    while(divi>1) {
        if (!(size % divi))
            break;
        divi--;
    }

    int yProcs=divi, xProcs=size/divi;
    if(!rank) {
        printf("%d Procs: [%d in x | %d in y]\n",size, xProcs,yProcs);
        int l;
        for(int k=0;k<xProcs;k++) {
            for(l=0;l<yProcs;l++) {
                printf("%2d ",l+yProcs*k);
            }
            printf("\n");
        }
    }






    MPI_Finalize();

}





//<editor-fold desc="Utility">
void displayMyMatrix(double *matrix, int m, int n) {
    printf(" \n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            //printf("<%d %d %f>",i,j,matrix[i*m+j]);
            printf("%f ", matrix[i * m + j]);
        }
        printf(" \n");
    }
}

void saveMyMatrix(double *matrix, int m, int n, double h, int numberTask) {
    // h=1 for save indices
    std::ofstream myfile;

    if (numberTask == 0)
        myfile.open("./results_a.dat");
    else if (numberTask == 1)
        myfile.open("./results_b.dat");
    else
        myfile.open("./results_temp.dat");


    double x;
    double y;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            x = h * i;
            y = h * j;
            // printf("<%d %d %f>",x,y,matrix[i*m+j]);
            myfile << x << " " << y << " " << matrix[i * m + j] << "\n";
        }
        myfile << std::endl;
        // printf(" \n");
    }
    myfile.close();
}
//</editor-fold>




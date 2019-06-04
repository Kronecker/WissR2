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
#define SAVE_MATRIX 1
#endif
#ifndef MPI_PARALLEL
#define MPI_PARALLEL 0
#endif
#ifndef SHOW_PROC_LAYOUT
#define SHOW_PROC_LAYOUT 1
#endif
#ifndef JACOBI_ITERATIONS
#define JACOBI_ITERATIONS 2000
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

    double end, start; // MPI Timestamps

    start = MPI_Wtime();


    double *currentIteration, *lastIteration, *f;
    double h;

    h = 1. / (INNER_GRID_SIZE + 1.);

    currentIteration = prepareGridAndExtendWithBoundaryVals(n, 0);
    lastIteration = prepareGridAndExtendWithBoundaryVals(n, 0);

    f = initMatrixRightHandSideSeriell(n, h);

    jacobiIterationWithExtGridSeriell(n, f, currentIteration, lastIteration, JACOBI_ITERATIONS, h);


    end = MPI_Wtime();
    printf("%3d SERIELL node in %.1fms\n", 1, (end - start) * 1000);


#if SHOW_MATRIX
    displayMyMatrix(lastIteration,n,n);
#endif

#if SAVE_MATRIX
    saveMyMatrix(lastIteration, n, n, h, 42);
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
    int i, k, index;

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
        temp = gridCurrentIteration;
        gridCurrentIteration = gridLastIteration;
        gridLastIteration = temp;
    }


}
//</editor-fold>


double *initMatrixRightHandSideParallel(int startX, int numElX, int startY, int numElY, double h);

void
saveMyMatrixMPI(double *matrix, int m, int n, int offsetM, int offsetN, double h, int numberTask, int rank, int size,
                int nn_left, int nn_up, int nn_right, int nn_down);

double *prepareGridParallel(int numElX, int numElY);

void jacobiIterationWithExtGridParallel(int numElX, int numElY, double *f, double *gridCurrentIteration,
                                        double *gridLastIteration,
                                        int numberOfIterations, double h, int nn_left, int nn_up, int nn_right,
                                        int nn_down, int rank);

void commOverBoundaries(int nn_oneD, int nn_otherD, double *dataSend, double *dataRecv, MPI_Datatype dataType,
                        int numOfDatas, int rank);

int main_jacobiMPI(int argc, char *argv[]) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double end, start; // MPI Timestamps
    if (rank == 0) {  // Master measures time
        start = MPI_Wtime();
    }


    int n = INNER_GRID_SIZE + 2;
    double *currentIteration, *lastIteration, *f;
    double h;


    h = 1. / (INNER_GRID_SIZE + 1.);

    // reverse trial division, find "most square" configuration
    // under the assumption that all procs should have to work, not necessary the best configuration
    // e.g. 17 scales very bad, [17 1], but with one proc idle, 16 scales perfect (square, [4 4])

    int divi = sqrt(size);
    while (divi > 1) {
        if (!(size % divi))
            break;
        divi--;
    }

    // number of procs in each dim
    int yProcs = divi, xProcs = size / divi;


    // determin next neighbours
    int nn_left, nn_up, nn_right, nn_down;
    if (!(rank % yProcs)) {
        nn_left = -1;
    } else {
        nn_left = rank - 1;
    }

    if ((rank % yProcs) == (yProcs - 1)) {
        nn_right = -1;
    } else {
        nn_right = rank + 1;
    }

    if (rank < yProcs) {
        nn_up = -1;
    } else {
        nn_up = rank - yProcs;
    }

    if (rank > size - yProcs - 1) {
        nn_down = -1;
    } else {
        nn_down = rank + yProcs;
    }

    // evenly distribute grid points
    int numElY, numElX;
    int restEl;
    int startY, startX;

    numElY = INNER_GRID_SIZE / yProcs;
    restEl = INNER_GRID_SIZE % yProcs;

    if ((rank % yProcs) < restEl) {  // take one for the team
        startY = (rank % yProcs) * (numElY) + rank % yProcs;
        numElY++;
    } else {                  // get away with less
        startY = (rank % yProcs) * (numElY) + restEl;
    }

    numElX = INNER_GRID_SIZE / xProcs;
    restEl = INNER_GRID_SIZE % xProcs;

    if ((rank / yProcs) < restEl) {  // take one for the team
        startX = (rank / yProcs) * numElX + rank / yProcs;
        numElX++;
    } else {                   // get away with less
        startX = (rank / yProcs) * numElX + restEl;
    }

    // Add boundaries
    numElX += 2;
    numElY += 2;

    // Init Vectors with 0
    currentIteration = prepareGridParallel(numElX, numElY);
    lastIteration = prepareGridParallel(numElX, numElY);

    f = initMatrixRightHandSideParallel(startX, numElX, startY, numElY, h);


#if SHOW_PROC_LAYOUT
    if (!rank) {
        printf("%d Procs: [%d in x | %d in y]\n", size, xProcs, yProcs);
        int l;
        for (int k = 0; k < xProcs; k++) {
            for (l = 0; l < yProcs; l++) {
                printf("%2d ", l + yProcs * k);
            }
            printf("\n");
        }
    }
#endif


    jacobiIterationWithExtGridParallel(numElX, numElY, f, currentIteration, lastIteration,
                                       JACOBI_ITERATIONS, h, nn_left, nn_up, nn_right, nn_down, rank);


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {  // Master measures time
        end = MPI_Wtime();
        printf("%3d MPI nodes in %.1fms\n", size, (end - start) * 1000);
    }


#if SAVE_MATRIX
    saveMyMatrixMPI(lastIteration, numElX, numElY, startX, startY, h, 42, rank, size, nn_left, nn_up, nn_right,
                    nn_down);
#endif

    delete (f);
    delete (currentIteration);
    delete (lastIteration);


    MPI_Finalize();

}


/*
          for (int k = 0; k < size; k++) {
            if (rank == k) {
                saveMyMatrixOffset(f, numElX, numElY, startX, startY, h, 42);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

 */

void jacobiIterationWithExtGridParallel(int numElX, int numElY, double *f, double *gridCurrentIteration,
                                        double *gridLastIteration,
                                        int numberOfIterations, double h, int nn_left, int nn_up, int nn_right,
                                        int nn_down, int rank) {
    double *temp;
    MPI_Datatype column;
    MPI_Status status;

    MPI_Type_vector(numElX, 1, numElY, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);




    // some caching
    double hSquare = h * h;
    double oneBy4 = 1 / 4.;
    int xm1 = numElX - 1;
    int ym1 = numElY - 1;

    int i, k, index;

    for (int l = 0; l < numberOfIterations; l++) {

        for (k = 1; k < xm1; k++) {
            for (i = 1; i < ym1; i++) {
                index = k * numElY + i;
                gridCurrentIteration[index] = oneBy4 * (hSquare * f[index]
                                                        + (gridLastIteration[index - numElY] +
                                                           gridLastIteration[index - 1] +
                                                           gridLastIteration[index + 1] +
                                                           gridLastIteration[index + numElY]));
            }
        }


        commOverBoundaries(nn_left, nn_right, &gridCurrentIteration[1], &gridCurrentIteration[numElY - 1], column, 1,
                           rank);
        commOverBoundaries(nn_up, nn_down, &gridCurrentIteration[numElY], &gridCurrentIteration[(numElX - 1) * numElY],
                           MPI_DOUBLE, numElY, rank);
        commOverBoundaries(nn_right, nn_left, &gridCurrentIteration[numElY - 2], &gridCurrentIteration[0], column, 1,
                           rank);
        commOverBoundaries(nn_down, nn_up, &gridCurrentIteration[(numElX - 2) * numElY], &gridCurrentIteration[0],
                           MPI_DOUBLE, numElY, rank);

        temp = gridCurrentIteration;
        gridCurrentIteration = gridLastIteration;
        gridLastIteration = temp;
    }

    MPI_Type_free(&column);
}

void commOverBoundaries(int nn_oneD, int nn_otherD, double *dataSend, double *dataRecv, MPI_Datatype dataType,
                        int numOfDatas, int rank) {

    MPI_Status status;

    if (rank % 2) {
        if (nn_oneD >= 0) {
            MPI_Send(dataSend, numOfDatas, dataType, nn_oneD, 0, MPI_COMM_WORLD);
        }
    } else {
        if (nn_otherD >= 0) {
            MPI_Recv(dataRecv, numOfDatas, dataType, nn_otherD, 0, MPI_COMM_WORLD, &status);
        }
    }

    if (!(rank % 2)) {
        if (nn_oneD >= 0) {
            MPI_Send(dataSend, numOfDatas, dataType, nn_oneD, 0, MPI_COMM_WORLD);
        }
    } else {
        if (nn_otherD >= 0) {
            MPI_Recv(dataRecv, numOfDatas, dataType, nn_otherD, 0, MPI_COMM_WORLD, &status);
        }
    }


}


double *initMatrixRightHandSideParallel(int startX, int numElX, int startY, int numElY, double h) {
    double *matrix = new double[numElX * numElY];
    double x;
    double y;
    int j;
    for (int i = 0; i < numElX; i++) {

        for (j = 0; j < numElY; j++) {
            x = h * (i + startX);
            y = h * (j + startY);
            matrix[i * numElY + j] = x * (1 - x) + y * (1 - y);
            //printf("<%f %f> %f\n", x, y, matrix[i * numElY + j]);
        }
    }
    return matrix;
}


double *prepareGridParallel(int numElX, int numElY) {

    double *matrix = new double[numElX * numElY]();

    return matrix;
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

void
saveMyMatrixMPI(double *matrix, int m, int n, int offsetM, int offsetN, double h, int numberTask, int rank, int size,
                int nn_left, int nn_up, int nn_right, int nn_down) {
    // h=1 for save indices
    std::ofstream myfile;
    double x;
    double y;
    int startI = 0, startJ = 0;


    for (int k = 0; k < size; k++) {
        if (rank == k) {
            if (!rank) {
                if (numberTask == 0)
                    myfile.open("./results_a.dat", std::ios::trunc);
                else if (numberTask == 1)
                    myfile.open("./results_b.dat", std::ios::trunc);
                else
                    myfile.open("./results_temp.dat", std::ios::trunc);
            } else {
                if (numberTask == 0)
                    myfile.open("./results_a.dat", std::ios::app);
                else if (numberTask == 1)
                    myfile.open("./results_b.dat", std::ios::app);
                else
                    myfile.open("./results_temp.dat", std::ios::app);
            }


            if (!nn_left < 0) {   // ignore ghost boundaries
                startJ = 1;
            }
            if (!nn_up < 0) {
                startI = 1;
            }
            if (!nn_right < 0) {
                n--;
            }
            if (!nn_down < 0) {
                m--;
            }


            for (int i = startI; i < m; i++) {
                for (int j = startJ; j < n; j++) {
                    x = h * (i + offsetM);
                    y = h * (j + offsetN);
                    // printf("<%d %d %f>",x,y,matrix[i*m+j]);
                    myfile << x << " " << y << " " << matrix[i * n + j] << "\n";
                }
                myfile << std::endl;
                // printf(" \n");
            }
            myfile.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


}


//</editor-fold>



/*
 *
 *
 *         // Comm Left
        if(rank%2) {
            if (nn_left >= 0) {
                printf("%d sending to %d\n",rank,nn_left);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_left, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_right >= 0) {
                printf("%d receiving from %d\n",rank,nn_right);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_right, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_left > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_left, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_right > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_right, 0, MPI_COMM_WORLD, &status);
            }
        }



        // Comm Up
        if(rank%2) {
            if (nn_up >= 0) {
                printf("%d sending to %d\n",rank,nn_up);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_up, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_down >= 0) {
                printf("%d receiving from %d\n",rank,nn_down);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_down, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_up > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_up, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_down > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_down, 0, MPI_COMM_WORLD, &status);
            }
        }


        // Comm Right
        if(rank%2) {
            if (nn_right >= 0) {
                printf("%d sending to %d\n",rank,nn_right);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left >= 0) {
                printf("%d receiving from %d\n",rank,nn_left);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_right > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        // Comm Down

        // Comm Right
        if(rank%2) {
            if (nn_right >= 0) {
                printf("%d sending to %d\n",rank,nn_right);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left >= 0) {
                printf("%d receiving from %d\n",rank,nn_left);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_right > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }
        // Comm Left
        if(rank%2) {
            if (nn_left >= 0) {
                printf("%d sending to %d\n",rank,nn_left);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_left, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_right >= 0) {
                printf("%d receiving from %d\n",rank,nn_right);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_right, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_left > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_left, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_right > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_right, 0, MPI_COMM_WORLD, &status);
            }
        }



        // Comm Up
        if(rank%2) {
            if (nn_up >= 0) {
                printf("%d sending to %d\n",rank,nn_up);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_up, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_down >= 0) {
                printf("%d receiving from %d\n",rank,nn_down);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_down, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_up > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_up, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_down > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_down, 0, MPI_COMM_WORLD, &status);
            }
        }


        // Comm Right
        if(rank%2) {
            if (nn_right >= 0) {
                printf("%d sending to %d\n",rank,nn_right);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left >= 0) {
                printf("%d receiving from %d\n",rank,nn_left);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_right > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        // Comm Down

        // Comm Right
        if(rank%2) {
            if (nn_right >= 0) {
                printf("%d sending to %d\n",rank,nn_right);
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left >= 0) {
                printf("%d receiving from %d\n",rank,nn_left);
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }

        if(!(rank%2)) {
            if (nn_right > 0) {
                MPI_Send(&gridCurrentIteration[1], 1, column, nn_right, 0, MPI_COMM_WORLD);
            }
        } else {
            if (nn_left > 0) {
                MPI_Recv(&gridCurrentIteration[numElY - 1], 1, column, nn_left, 0, MPI_COMM_WORLD, &status);
            }
        }


 */
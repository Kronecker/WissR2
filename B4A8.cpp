//
// Created by Looky on 14.06.2019.
//

#include <iostream>
#include <math.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>



#ifndef INNER_GRID_SIZE
#define INNER_GRID_SIZE 1024
#endif
#ifndef SHOW_MATRIX
#define SHOW_MATRIX 0
#endif
#ifndef SAVE_MATRIX
#define SAVE_MATRIX 1
#endif
#ifndef MAX_ITER
#define MAX_ITER 1000
#endif
#ifndef PARALLEL
#define PARALLEL 1
#endif








// Seriell
void saveMyMatrix(double *matrix, int m, int n, double h, int numberTask);
double* prepareGridAndExtendWithBoundaryVals(int n);
double* initMatrixRightHandSideSeriell(int n, double h);
double scalProd(double* vecA, double* vecB, int length);
void matrixMult5PointStar(double* rightSideVec, double* resultVec, int n);
void fMA (double* vecA, double scalar, double* vecB, double *resultVec, int length);
void displayMyMatrix(double *matrix, int m, int n);



// Parallel
saveMyMatrixMPI(double *matrix, int m, int n, int offsetM, int offsetN, double h, int numberTask, int rank, int size,
                int nn_left, int nn_up, int nn_right, int nn_down)





int main() {
    std::cout << "Hello there!" << std::endl;
#if PARALLEL
    main_Parallel()
#else
    main_Seriell();
#endif

}

int main_Seriell() {

    int n = INNER_GRID_SIZE + 2;
    int nSquare=n*n;

    double *currentIteration,*fHsquare;
    double h;


    h = 1. / (INNER_GRID_SIZE + 1.);

    printf("Prep Grid\n");
    currentIteration = prepareGridAndExtendWithBoundaryVals(n);

   printf("Init right hand side\n");
   fHsquare = initMatrixRightHandSideSeriell(n, h);

    double *conjVec, *residual;
    double alpha, beta;
    double *cacheATimesConjVec, cacheConjVecTimesATimesConjVec;


    conjVec=new double[n*n];
    residual=new double[n*n];
    cacheATimesConjVec=new double[n*n];
    printf("Init residual and conj. direction\n");
    // Start is: A*x(0)-b ....set x(0)=0  -> r(0)=-b  p(0)=r(0)
    for(int k=0;k<nSquare;k++) {
        conjVec[k]=residual[k]=-fHsquare[k];
    }

    printf("Start Iteration\n\n");
    for(int k=0;k<MAX_ITER;k++) {

        matrixMult5PointStar(conjVec, cacheATimesConjVec, n);

        cacheConjVecTimesATimesConjVec=scalProd(conjVec,cacheATimesConjVec,nSquare);

        alpha=scalProd(conjVec,residual,nSquare)/cacheConjVecTimesATimesConjVec;

        fMA(currentIteration,alpha,conjVec,currentIteration,nSquare);

        fMA(residual, -alpha, cacheATimesConjVec, residual, nSquare);

        beta=scalProd(residual,cacheATimesConjVec, nSquare)/cacheConjVecTimesATimesConjVec;

        fMA(residual, -beta, conjVec, conjVec, nSquare);
    }

#if SAVE_MATRIX
    saveMyMatrix(currentIteration, n, n, h,1);

#endif



   delete(residual);
   delete(conjVec);
   delete(fHsquare);
   delete(currentIteration);



    return 0;
}

double* prepareGridAndExtendWithBoundaryVals(int n)  {
    double *matrix = new double[n * n]();

    return matrix;
}

void fMA (double* vecA, double scalar, double* vecB, double *resultVec, int length) {
    for(int k=0;k<length;k++) {
        resultVec[k] = vecA[k] + scalar * vecB[k];
    }
}


double *initMatrixRightHandSideSeriell(int n,
                                       double h) {   // parts of the matrix correspond to boundary values and are irrelevant
    double *matrix = new double[n * n];
    double x;
    double y;
    double hSquare=h*h;

    for (int i = 1; i < n-1; i++) {

        for (int j = 1; j < n-1; j++) {
            x = h * i;
            y = h * j;
            matrix[i * n + j] = (x * (1 - x) + y * (1 - y))*hSquare;
            //         printf("<%f %f> %f\n",x,y,matrix[i*m+j]);
        }
    }
    return matrix;
}

double scalProd(double* vecA, double* vecB, int length) {
    double sum=0;
    for(int k=0;k<length;k++)
        sum+=vecA[k]*vecB[k];
    return sum;
}

void matrixMult5PointStar(double* rightSideVec, double* resultVec, int n) {
    int nm1=n-1;
    int i,k, index;


    k=0;
    for (i = 1; i < nm1; i++) {
        index = k * n + i;
        resultVec[index] =  0;
    }
    k=nm1;
    for (i = 1; i < nm1; i++) {
        index = k * n + i;
        resultVec[index] =  0;
    }
    i=0;
    for (k = 1; k < nm1; k++) {
            index = k * n + i;
            resultVec[index] =  0;
    }
    i=nm1;
    for (k = 1; k < nm1; k++) {
        index = k * n + i;
        resultVec[index] =  0;
    }

    for (k = 1; k < nm1; k++) {
        for (i = 1; i < nm1; i++) {
            index = k * n + i;

            resultVec[index] =  4 * rightSideVec[index] -
                                rightSideVec[index - n] -
                                rightSideVec[index - 1] -
                                rightSideVec[index + 1] -
                                rightSideVec[index + n];

        }
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








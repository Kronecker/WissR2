//
// Created by looky on 24.06.2019.
//


#include <iostream>
#include <math.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#ifndef INNER_GRID_SIZE
#define INNER_GRID_SIZE 101
#endif
#ifndef INNER_GRID_SIZE_LEFT_SIDE
#define INNER_GRID_SIZE_LEFT_SIDE 51
#endif
#ifndef SAVE_MATRIX
#define SAVE_MATRIX 1
#endif
#ifndef JACOBI_ITERATIONS
#define JACOBI_ITERATIONS 200
#endif
#ifndef SCHWARZ_ITERATIONS
#define SCHWARZ_ITERATIONS 50
#endif



void saveMyMatrix(double *matrix, int m, int n, double h, int numberTask);



void gsIteration(int nx, int ny, double *f, double *gridCurrentIteration,
                 int numberOfIterations, double h);
void pushBoundariesToRight(int nx ,int ny, double *gridLeft,double *gridRight);
void pushBoundariesToLeft(int nx ,int ny, double *gridRight,double *gridLeft);


int main() {


    if (INNER_GRID_SIZE < INNER_GRID_SIZE_LEFT_SIDE ) {
        std::cout << "Grid size is too small" << std::endl;
        return 0;
    }
    if (INNER_GRID_SIZE > 2*INNER_GRID_SIZE_LEFT_SIDE ) {
        std::cout << "Overlapp is too small" << std::endl;
        return 0;
    }

    int nx = INNER_GRID_SIZE_LEFT_SIDE + 2;
    int ny = INNER_GRID_SIZE + 2;

    int numel = nx * ny;

    double *gridLeft, *gridRight;
    double h;


    gridLeft = new double[numel];

    gridRight = new double[numel];



    h = 1. / (INNER_GRID_SIZE + 1.);


    double* rhSide_Left = new double[numel];
    double* rhSide_Right = new double[numel];


    for (int k = 0; k < numel; k++) {
        rhSide_Left[k] = 1;
        rhSide_Right[k] = 1;
    }




    for(int k=0;k<SCHWARZ_ITERATIONS;k++) {

        gsIteration(nx, ny, rhSide_Left, gridLeft, JACOBI_ITERATIONS, h);
        pushBoundariesToRight(nx, ny, gridLeft, gridRight);
        gsIteration(nx, ny, rhSide_Right, gridRight, JACOBI_ITERATIONS, h);
        pushBoundariesToLeft(nx, ny, gridRight, gridLeft);

    }



#if SAVE_MATRIX
    saveMyMatrix(gridLeft, ny, nx, h, 0);
    saveMyMatrix(gridRight, ny, nx, h, 1);
#endif


    delete (rhSide_Left);
    delete (rhSide_Right);
    delete (gridLeft);
    delete (gridRight);

}


void pushBoundariesToRight(int nx ,int ny, double *gridLeft,double *gridRight) {
    for(int k=0;k<ny;k++) {
        gridRight[k*nx]=gridLeft[(nx-1)+(k*nx)];
    }
}

void pushBoundariesToLeft(int nx ,int ny, double *gridRight,double *gridLeft) {
    for(int k=0;k<ny;k++) {
        gridLeft[(nx)+(k*nx)]=gridRight[k*nx];
    }
}



void gsIteration(int nx, int ny, double *f, double *gridCurrentIteration,
                                       int numberOfIterations, double h) {
    // some caching
    double hSquare = h * h;
    double oneBy4 = 1 / 4.;
    int nxm1 = nx - 1;
    int nym1 = ny - 1;
    int i, k, index;

    for (int l = 0; l < numberOfIterations; l++) {

        for (k = 1; k < nym1; k++) {
            for (i = 1; i < nxm1; i++) {
                index = k * nx + i;
                gridCurrentIteration[index] = oneBy4 * (hSquare * f[index]
                                                        + (gridCurrentIteration[index - nx] +
                                                           gridCurrentIteration[index - 1] +
                                                           gridCurrentIteration[index + 1] +
                                                           gridCurrentIteration[index + nx]));
            }
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
            myfile << x << " " << y << " " << matrix[i * n + j] << "\n";
        }
        myfile << std::endl;
        // printf(" \n");
    }
    myfile.close();
}



















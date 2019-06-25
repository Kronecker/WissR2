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
#ifndef SOR_ITERATIONS
#define SOR_ITERATIONS 200
#endif
#ifndef SCHWARZ_ITERATIONS
#define SCHWARZ_ITERATIONS 50
#endif
#ifndef RELAXATION_FACTOR
#define RELAXATION_FACTOR 1.75
#endif



void saveMyMatrix(double *matrix, int m, int n, double h, int numberTask);



void sorIteration(int nx, int ny, double *f, double *gridCurrentIteration,
                 int numberOfIterations, double h);
void pushBoundariesToRight(int nx ,int ny, double *gridLeft,double *gridRight, int overlap);
void pushBoundariesToLeft(int nx ,int ny, double *gridRight,double *gridLeft, int overlap);
void saveSlice (double* slice, int n, int create, int nx);



int main(int argc,char *argv[]) {


    int innerGridSizeLeftSide=INNER_GRID_SIZE_LEFT_SIDE;
    int schwarz_iterations=SCHWARZ_ITERATIONS;

    if(argc>1)  {

        innerGridSizeLeftSide=atoi(argv[1]);
        if (INNER_GRID_SIZE < innerGridSizeLeftSide ) {
            std::cout << "StartParam: Grid size is too small" << std::endl;
            return 0;
        }
        if (INNER_GRID_SIZE > 2*innerGridSizeLeftSide ) {
            std::cout << "StartParam: Overlapp is too small" << std::endl;
            return 0;
        }

        if(argc>2)  {


            schwarz_iterations=atoi(argv[2]);

        }



    } else {


        if (INNER_GRID_SIZE < INNER_GRID_SIZE_LEFT_SIDE) {
            std::cout << "MakroParam: Grid size is too small" << std::endl;
            return 0;
        }
        if (INNER_GRID_SIZE > 2 * INNER_GRID_SIZE_LEFT_SIDE) {
            std::cout << "MakroParam: Overlapp is too small" << std::endl;
            return 0;
        }

    }

    std::cout<<"Alternierender Schwarz mit Nx = "<<  innerGridSizeLeftSide  << " und " << schwarz_iterations << " Schwarz-Iterationen."<<std::endl;


    int nx = innerGridSizeLeftSide + 2;
    int ny = INNER_GRID_SIZE + 2;
    int overlap=2*innerGridSizeLeftSide-INNER_GRID_SIZE;

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


    double *slice=new double[schwarz_iterations];

    for(int k=0;k<schwarz_iterations;k++) {

        sorIteration(nx, ny, rhSide_Left, gridLeft, SOR_ITERATIONS, h);
        pushBoundariesToRight(nx, ny, gridLeft, gridRight, overlap);
        sorIteration(nx, ny, rhSide_Right, gridRight, SOR_ITERATIONS, h);
        pushBoundariesToLeft(nx, ny, gridRight, gridLeft, overlap);
        slice[k]=gridRight[ny/2*nx+overlap/2];
    }
    saveSlice(slice,schwarz_iterations,0,innerGridSizeLeftSide);


#if SAVE_MATRIX
    saveMyMatrix(gridLeft, ny, nx, h, 0);
    saveMyMatrix(gridRight, ny, nx, h, 1);
#endif


    delete (rhSide_Left);
    delete (rhSide_Right);
    delete (gridLeft);
    delete (gridRight);

}


void pushBoundariesToRight(int nx ,int ny, double *gridLeft,double *gridRight, int overlap) {
    for(int k=0;k<ny;k++) {
        gridRight[k*nx]=gridLeft[(nx-3-overlap)+(k*nx)];  // -3 == -1 for boundary, -2 for inner<->outer
    }
}

void pushBoundariesToLeft(int nx ,int ny, double *gridRight,double *gridLeft, int overlap) {
    for(int k=0;k<ny;k++) {
        gridLeft[(nx-1)+(k*nx)]=gridRight[overlap+1+k*nx];
    }
}



void sorIteration(int nx, int ny, double *f, double *gridCurrentIteration,
                                       int numberOfIterations, double h) {
    // some caching
    double hSquare = h * h;
    double oneBy4 = 1 / 4.;
    int nxm1 = nx - 1;
    int nym1 = ny - 1;
    int i, k, index;
    double relaxa=RELAXATION_FACTOR;  // relaxation factor

    for (int l = 0; l < numberOfIterations; l++) {

        for (k = 1; k < nym1; k++) {
            for (i = 1; i < nxm1; i++) {
                index = k * nx + i;
                gridCurrentIteration[index] = oneBy4 * relaxa* (hSquare * f[index]
                                                        + (gridCurrentIteration[index - nx] +
                                                           gridCurrentIteration[index - 1] +
                                                           gridCurrentIteration[index + 1] +
                                                           gridCurrentIteration[index + nx]))
                                               + gridCurrentIteration[index]*(1-relaxa);
            }
        }
    }


}

void saveSlice (double* slice, int n, int create, int nx) {

    std::ofstream myfile;

        if(create) {
            myfile.open("./slice.dat",std::ios::trunc);
        }else {
            myfile.open("./slice.dat",std::ios::app);
        }

        myfile << "\"Nx = " << nx << "\""<<std::endl;
        for (int j = 0; j < n; j++) {

            myfile << j+1 << "\t" << slice[j] <<std::endl;
        }
        myfile<<std::endl<<std::endl;


    myfile.close();

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



















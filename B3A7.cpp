//
// Created by Looky on 02/06/2019.
//



#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef INNER_GRID_SIZE
    #define INNER_GRID_SIZE 1024
#endif


double* initMatrixRightHandSideSeriell(int n, double h  );
void jacobiIterationWithExtGridSeriell(int n, double *f, double valBoundary, int* numberOfIterations, double h);
void displayMyMatrix(double* matrix, int m,int n);
void saveMyMatrix(double* matrix, int m,int n, double h, int numberTask);

int main_jacobiSeriell (int argc,char *argv[]) {

    int n=INNER_GRID_SIZE+2;
    double* matrix

  //  initMatrixRightHandSide

    matrix=prepareGridAndExtendWithBoundaryVals(n,1);


    displayMyMatrix(matrix,n,n);


}

double* prepareGridAndExtendWithBoundaryVals(int n, double boundary) {

    double* matrix=new double[n*n]();

    // boundary values init (outer)
    for(int i=0;i<n;i++) {
        matrix[i]=boundary;
        matrix[n*(n-1)+i]=boundary;
    }
    for(int k=1;k<n-1;k++) { // iterate through blocks
        matrix[k*n]=boundary;
        matrix[(k+1)*n-1]=boundary;
    }

    return matrix;
}




void jacobiIterationWithExtGridSeriell(int n, double *f, double* gridActualIteration, double* gridLastIteration, int* numberOfIterations, double h) {






}



double* initMatrixRightHandSide(int n, double h  ) {
    double*matrix=new double[n*n];
    double x;
    double y;
    for (int i=0;i<n;i++) {

        for (int j=0;j<n;j++) {
            x=h*i;
            y=h*j;
            matrix[i*n+j]=x*(1-x)+y*(1-y);
            //         printf("<%f %f> %f\n",x,y,matrix[i*m+j]);
        }
    }
    return matrix;
}




void displayMyMatrix(double* matrix, int m,int n) {
    printf(" \n");
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            //printf("<%d %d %f>",i,j,matrix[i*m+j]);
            printf("%f ",matrix[i*m+j]);
        }
        printf(" \n");
    }
}

void saveMyMatrix(double* matrix, int m,int n, double h, int numberTask) {
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
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            x=h*i;
            y=h*j;
            // printf("<%d %d %f>",x,y,matrix[i*m+j]);
            myfile<<x<<" "<<y<<" "<<matrix[i*m+j]<<"\n";
        }
        myfile<<std::endl;
        // printf(" \n");
    }
    myfile.close();
}


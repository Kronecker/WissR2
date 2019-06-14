//
// Created by grabiger on 14.06.2019.
//

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


saveMyMatrix(double *matrix, int m, int n, double h, int numberTask);

int main() {
    std::cout << "Hello, World!" << std::endl;




    return 0;
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

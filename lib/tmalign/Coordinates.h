//
// Created by Martin Steinegger on 1/15/21.
//

#ifndef STRUCCLUST_COORDINATES_H
#define STRUCCLUST_COORDINATES_H
#include "simd.h"

struct Coordinates{
    Coordinates(int size){
        x =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
        y =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
        z =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
        allocated = true;
    }
    Coordinates(){
        allocated = false;
    }
    ~Coordinates(){
        if(allocated == true){
            free(x);
            free(y);
            free(z);
        }
    }
    bool allocated;
    float * x;
    float * y;
    float * z;
};
#endif //STRUCCLUST_COORDINATES_H

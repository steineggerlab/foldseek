//
// Created by Martin Steinegger on 1/15/21.
//

// #ifndef STRUCCLUST_COORDINATES_H
// #define STRUCCLUST_COORDINATES_H
// #include "simd.h"
// #include <cstring>

// struct Coordinates{
    // Coordinates(int size) : size(size) {
    //     x =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //     y =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //     z =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //     allocated = true;
    // }
    // Coordinates(){
        // allocated = false;
    // }
    // ~Coordinates(){
        // if(allocated == true){
        // free(x);
        // free(y);
        // free(z);
        // }
    // }
    // bool allocated;
    // float * x;
    // float * y;
    // float * z;
    // int size;

    // void reallocate(int newsize){
    //     if (allocated == false) {
    //         x =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //         y =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //         z =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    //         allocated = true;
    //         size = newsize;
    //     } else {
    //         if (newsize > size) {
    //             float* new_x =(float*) mem_align(ALIGN_FLOAT, (newsize+VECSIZE_FLOAT)*sizeof(float));
    //             float* new_y =(float*) mem_align(ALIGN_FLOAT, (newsize+VECSIZE_FLOAT)*sizeof(float));
    //             float* new_z =(float*) mem_align(ALIGN_FLOAT, (newsize+VECSIZE_FLOAT)*sizeof(float));

    //             if (x) memcpy(new_x, x, size * sizeof(float));
    //             if (y) memcpy(new_y, y, size * sizeof(float));
    //             if (z) memcpy(new_z, z, size * sizeof(float));

    //             free(x);
    //             free(y);
    //             free(z);

    //             x = new_x;
    //             y = new_y;
    //             z = new_z;
    //             size = newsize;
    //         }
    //     }
    // }
// };
// #endif //STRUCCLUST_COORDINATES_H
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
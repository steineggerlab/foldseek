#pragma once
#include <iostream>
#include <string>

struct Atom {
    float x;
    float y;
    float z;
    char structure='x';        // 'x' : default, 'H' : helix, 'S' : sheet

    // 기능 2: pLDDT (B-factor column)
    float bfactor = 0.0f;

    // 기능 1: interface region
    bool is_interface = false;

    // 기능 4: UTMatrix 정렬 구조
    bool is_aligned = false;

    // 기능 5: MSA conservation score (-1 = 미설정)
    float conservation_score = -1.0f;

    // 기능 6: 잔기 정보 hover
    int residue_number = -1;
    std::string residue_name;

    Atom(float x_, float y_, float z_) : x(x_), y(y_), z(z_), structure{'x'} {}
    Atom(float x_, float y_, float z_, char c) : x(x_), y(y_), z(z_), structure{c} {}
    Atom() : x(0), y(0), z(0), structure{'x'} {}

    void set_position(float x_, float y_, float z_) {
        x = x_;
        y = y_;
        z = z_;
        return;
    }

    float* get_position() const {
        static float coords[3]; 
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
        return coords;
    }

    void set_structure(char c){
        if (c == 'x' || c == 'H' || c == 'S'){
            structure = c;
            return;
        }
    }

    char get_structure() const {
        return structure;
    }

    void print_position() {
        std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
        return;
    }
};

#include <iostream>
#include "PulchraWrapper.h"

PulchraWrapper::PulchraWrapper(){}

PulchraWrapper::~PulchraWrapper(){
    for (size_t i=0;i<xcaVec.size();i++){
      free(xcaVec.data()[i]);
      free(nVec.data()[i]);
      free(cVec.data()[i]);
    }
}

void PulchraWrapper::rebuildBackbone(Vec3 * ca, Vec3 * n, Vec3 * c, char * ami, size_t chainLen){
    double ** _xca;
    double ** _c;
    double ** _n;

    // Ca coordinate arrays need to be longer than chainLen for pulchra
    if (chainLen+2*offset > xcaVec.size()){
        size_t prevSize =xcaVec.size(); 

        xcaVec.resize(chainLen+2*offset);
        nVec.resize(chainLen+2*offset);
        cVec.resize(chainLen+2*offset);

        _xca = xcaVec.data();
        _c = nVec.data();
        _n = cVec.data();

        for (size_t i=prevSize;i<xcaVec.size();i++){
          _xca[i] = (double*)calloc(sizeof(double)*3,1);
        }
        for (size_t i=prevSize;i<xcaVec.size();i++){
          _n[i] = (double*)calloc(sizeof(double)*3,1);
        }
        for (size_t i=prevSize;i<xcaVec.size();i++){
          _c[i] = (double*)calloc(sizeof(double)*3,1);
        }

    } else {
        _xca = xcaVec.data();
        _c = nVec.data();
        _n = cVec.data();
    }

    // Load data into pulchra datastructure
    for (size_t i=0;i<chainLen;i++){
      _xca[i+offset][0] = ca[i].x;
      _xca[i+offset][1] = ca[i].y;
      _xca[i+offset][2] = ca[i].z;
    }
    for (size_t i=0;i<chainLen;i++){
      _n[i][0] = n[i].x;
      _n[i][1] = n[i].y;
      _n[i][2] = n[i].z;
    }
    for (size_t i=0;i<chainLen;i++){
      _c[i][0] = c[i].x;
      _c[i][1] = c[i].y;
      _c[i][2] = c[i].z;
    }
   
    // Run pulchra
    pulchra_rebuild_backbone(_xca, _n, _c, ami, chainLen); 
    
    // Write reconstructed coords. back
    for (size_t i=0;i<chainLen;i++){
      n[i].x = _n[i][0];
      n[i].y = _n[i][1];
      n[i].z = _n[i][2];
    }
    for (size_t i=0;i<chainLen;i++){
      c[i].x = _c[i][0];
      c[i].y = _c[i][1];
      c[i].z = _c[i][2];
    }
}

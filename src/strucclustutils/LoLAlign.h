#ifndef LoLAlign
#define LoLAlign

#include "SubstitutionMatrix.h"
#include "IndexReader.h"
#include "DBReader.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>

class LoLAlign{
public:
    typedef struct {
        float start_anchor_go;
        float start_anchor_ge;
        unsigned int start_anchor_temp;
        unsigned int start_anchor_length;
        float lol_go;
        float lol_ge;
        float lol_min_m;
        int lol_temp;
    } lol_params;

    void LoLAlign();
    void ~LoLAlign();
    void LoLscore_matrix(int d_ij, int d_kl, int seq_dist);



private:


};

#endif // LoLAligner




#include <math.h>
#include "basic_fun.h"

void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    //parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    D0_MIN=0.5; 
    dcu0=4.25;                       //update 3.85-->4.25
 
    Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
    if (Lnorm<=19)                    //update 15-->19
        d0=0.168;                   //update 0.5-->0.168
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    D0_MIN=0.5; 
 
    Lnorm=len;            //normaliz TMscore by this in searching
    if (Lnorm<=21) d0=0.5;          
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;   
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    d0=d_s;          
    Lnorm=len;            //normaliz TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;  
}

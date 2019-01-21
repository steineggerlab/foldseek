#include "TMalign.h"

/* entry function for se */
int se_main(
    double **xa, double **ya, const char *seqx, const char *seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const double Lnorm_ass, const double d0_scale,
    const bool a_opt, const bool u_opt, const bool d_opt)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double t[3]={0,0,0};  // dummy translation vection
    double u[3][3]={{1,0,0},{0,1,0},{0,0,1}}; // dummy rotation matrix
    double **score;       // Input score table for dynamic programming
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  

    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y

    /***********************/
    /* allocate memory     */
    /***********************/
    int minlen = min(xlen, ylen);
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);
    int *invmap          = new int[ylen+1];
    for(int i=0; i<ylen; i++) invmap[i]=-1;

    /* set d0 */
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
        score_d8, d0, d0_search, dcu0); // set score_d8
    parameter_set4final(xlen, D0_MIN, Lnorm,
        score_d8, d0B, d0_search, dcu0); // set d0B
    parameter_set4final(ylen, D0_MIN, Lnorm,
        score_d8, d0A, d0_search, dcu0); // set d0A
    if (a_opt)
        parameter_set4final((xlen+ylen)*0.5, D0_MIN, Lnorm,
            score_d8, d0a, d0_search, dcu0); // set d0a
    if (u_opt)
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            score_d8, d0u, d0_search, dcu0); // set d0u

    /* perform alignment by NW */
    NWDP_TM(score, path, val, xa, ya, xlen, ylen,
        t, u, d0*d0, 0, invmap);

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    for(int i=0,j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xa[i][0], &ya[j][0]));
            if (d <= score_d8)
            {
                m1[k]=i;
                m2[k]=j;
                k++;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);

    /* extract aligned sequence */
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    int kk=0, i_old=0, j_old=0;
    d=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';                    
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(dist(&xa[m1[k]][0], &ya[m2[k]][0]));
        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        kk++;  
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }    
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);

    /* free memory */
    delete [] invmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);
    return 0; // zero for no exception
}

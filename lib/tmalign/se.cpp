#include "se.h"

using namespace std;

void print_extra_help()
{
    cout <<
"Additional options: \n"
"    -dir     Perform all-against-all alignment extraction among the list\n"
"             of PDB chains listed by 'chain_list' under 'chain_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ se -dir chain1_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to perform alignment extraction from a list of\n"
"             PDB chains listed by 'chain1_list' under 'chain1_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ se -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain2 to perform alignment extraction from a list of\n"
"             PDB chains listed by 'chain2_list' under 'chain2_folder'\n"
"             $ se chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
"    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"             add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: (default in the first C++ TMalign) end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"    -outfmt  Output format\n"
"             0: (default) full output\n"
"             1: fasta format compact output\n"
"             2: tabular format very compact output\n"
"\n"
"    -infmt1  Input format for chain1\n"
"    -infmt2  Input format for chain2\n"
"             0: (default) PDB format\n"
"             1: SPICKER format\n"
"             2: xyz format\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    cout <<
"Extract sequence alignment from a pair of superposed structures.\n"
"\n"
"Usage: se PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -h    print the full help message.\n"
"\n"
"    (Options -u, -a, -d won't change the final structure alignment)\n"
"\n"
"Example usages:\n"
"    se PDB1.pdb PDB2.pdb\n"
"    se PDB1.pdb PDB2.pdb -u 100 -d 5.0 -a T\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = "";
    string yname       = "";
    double Lnorm_ass, d0_scale;

    bool A_opt = false; // marker for whether structure A is specified
    bool B_opt = false; // marker for whether structure B is specified
    bool h_opt = false; // print full help message
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    int    infmt1_opt=0;     // PDB format for chain_1
    int    infmt2_opt=0;     // PDB format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set

    int nameIdx = 0;
    for(int i = 1; i < argc; i++)
    {
        if ( (!strcmp(argv[i],"-u") || 
              !strcmp(argv[i],"-L")) && i < (argc-1) )
        {
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-a") && i < (argc-1) )
        {
            if (!strcmp(argv[i + 1], "T"))      a_opt=true;
            else if (!strcmp(argv[i + 1], "F")) a_opt=false;
            else PrintErrorAndQuit("Wrong value for option -a! It should be T or F");
            i++;
        }
        else if ( !strcmp(argv[i],"-d") && i < (argc-1) )
        {
            d0_scale = atof(argv[i + 1]); d_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if ( !strcmp(argv[i],"-infmt1") && i < (argc-1) )
        {
            infmt1_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-infmt2") && i < (argc-1) )
        {
            infmt2_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
        {
            ter_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-split") && i < (argc-1) )
        {
            split_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-atom") && i < (argc-1) )
        {
            atom_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir") && i < (argc-1) )
        {
            dir_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir1") && i < (argc-1) )
        {
            dir1_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir2") && i < (argc-1) )
        {
            dir2_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-suffix") && i < (argc-1) )
        {
            suffix_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-outfmt") && i < (argc-1) )
        {
            outfmt_opt=atoi(argv[i + 1]); i++;
        }
        else
        {
            if (nameIdx == 0)
            {
                xname=argv[i]; B_opt = true;
            }
            else if (nameIdx == 1)
            {
                yname=argv[i]; A_opt = true;
            }
            nameIdx++;
        }
    }

    if(!B_opt || (!A_opt && dir_opt.size()==0) || (A_opt && dir_opt.size()))
        if (h_opt) print_help(h_opt);

    if (!A_opt && dir_opt.size()==0) 
        PrintErrorAndQuit("Please provide structure A");
    if (!B_opt) 
        PrintErrorAndQuit("Please provide structure B");

    if (suffix_opt.size() && dir1_opt.size()==0 && dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if (dir_opt.size() && (dir1_opt.size() || dir2_opt.size()))
        PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! atom name must have 4 characters, including space.");
    if (u_opt && Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    
    /* parse file list */
    if (dir1_opt.size()==0 && dir_opt.size()==0)
        chain1_list.push_back(xname);
    else
    {
        ifstream fp(xname.c_str());
        if (! fp.is_open())
        {
            char message[5000];
            sprintf(message, "Can not open file: %s\n", xname.c_str());
            PrintErrorAndQuit(message);
        }
        string line;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            chain1_list.push_back(dir_opt+dir1_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    if (dir_opt.size())
    {
        for (int i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    }
    else if (dir2_opt.size()==0)
        chain2_list.push_back(yname);
    else
    {
        ifstream fp(yname.c_str());
        if (! fp.is_open())
        {
            char message[5000];
            sprintf(message, "Can not open file: %s\n", yname.c_str());
            PrintErrorAndQuit(message);
        }
        string line;
        while (fp.good())
        {
            getline(fp, line);
            if (! line.size()) continue;
            chain2_list.push_back(dir2_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    if (outfmt_opt==2)
        cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;

    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    int    *xresno, *yresno;   // residue number for fragment gapless threading
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    double t0[3]={0,0,0};
    double u0[3][3]={{1,0,0},{0,1,0},{0,0,1}};

    /* loop over file names */
    for (int i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            resi_vec1, 0, ter_opt, infmt1_opt, atom_opt, split_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (int chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            xresno = new int[xlen];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, xresno);

            for (int j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2,
                        chainID_list2, resi_vec2, 0, ter_opt,
                        infmt2_opt, atom_opt, split_opt);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (int chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    yresno = new int[ylen];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy, yresno);

                    /* declare variable specific to this pair of TMalign */
                    double TM1, TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;

                    /* entry function for structure alignment */
                    se_main(
                        xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, Lnorm_ass, d0_scale, a_opt, u_opt, d_opt);

                    /* print result */
                    output_results(
                        xname.substr(dir1_opt.size()).c_str(),
                        yname.substr(dir2_opt.size()).c_str(),
                        chainID_list1[chain_i].c_str(),
                        chainID_list2[chain_j].c_str(),
                        xlen, ylen, t0, u0, TM1, TM2, 
                        TM3, TM4, TM5, rmsd0, d0_out,
                        seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, n_ali, L_ali, TM_ali, rmsd_ali,
                        TM_0, d0_0, d0A, d0B,
                        Lnorm_ass, d0_scale, d0a, d0u, 
                        "", outfmt_opt, ter_opt, "",
                        false, false, a_opt, u_opt, d_opt);

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    delete [] yresno;
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (int chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    resi_vec2.clear();
                    chainID_list2.clear();
                }
            } // j
            DeleteArray(&xa, xlen);
            delete [] seqx;
            delete [] xresno;
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        resi_vec1.clear();
        chainID_list1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (int chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
    }
    chain1_list.clear();
    chain2_list.clear();
    return 0;
}

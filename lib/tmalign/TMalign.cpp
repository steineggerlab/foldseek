/*
===============================================================================
   This is a re-implementation of TM-align algorithm in C/C++. The code was
   is written by Jianyi Yang and later updated by Jianjie Wu at The Yang Zhang
   lab, Department of Computational Medicine and Bioinformatics, University of
   Michigan, 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218. Please report bugs
   and questions to zhng@umich.edu

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of
     the Software. It is provided "as is" without express or implied
     warranty.
   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by Jianyi Yang
   2016/05/21: Several updates of this program were made by Jianjie Wu, including
              (1) fixed several compiling bugs
              (2) made I/O of C/C++ version consistent with the Fortran version
              (3) added outputs including full-atom and ligand structures
              (4) added options of '-i', '-I' and '-m'
   2016/05/25: Fixed a bug on PDB file reading
   2018/06/04: Fixed a bug in PDB file with negative residue number. Added
               options -fast, -dir1, -dir2, -suffix, -atom, -ter, -outfmt.
               Re-write the file reading function to reduce the number of
               times a PDB file need to be read.
   2018/07/27: Added the -byresi option for TM-score superposition without
               re-alignment as in TMscore and TMscore -c
   2018/08/07: Added the -dir option
   2018/08/14: Added the -split option
   2018/08/16: Added the -infmt1, -infmt2 options.
               TMalign can now read .gz and .bz2 compressed files.
===============================================================================
*/
#include "TMalign.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" *****************************************************************************\n"
" * TM-align (Version 20180816): A protein structural alignment algorithm     *\n"
" * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *\n"
" * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *\n"
" *****************************************************************************"
    << endl;
}

void print_extra_help()
{
    cout <<
"Additional options: \n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
"    -dir     Perform all-against-all alignment among the list of PDB\n"
"             chains listed by 'chain_list' under 'chain_folder'. Note\n"
"             that the slash is necessary.\n"
"             $ TMalign -dir chain1_folder/ chain_list\n"
"\n"
"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"             under 'chain1_folder'. Note that the slash is necessary.\n"
"             $ TMalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"             under 'chain2_folder'\n"
"             $ TMalign chain1 -dir2 chain2_folder/ chain2_list\n"
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
"            -1: full output, but without version or citation information\n"
"\n"
"    -byresi  Whether to align two structures by residue index.\n"
"             0: (default) do not align by residue index\n"
"             1: (same as TMscore program) align by residue index\n"
"             2: (same as TMscore -c, should be used with -ter 1)\n"
"                align by residue index and chain ID\n"
"             3: (similar to TMscore -c, should be used with -ter 1)\n"
"                align by residue index and order of chain\n"
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
    print_version();
    cout <<
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length (the same as -L)\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n"
"\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n"
"\n"
"    -i    Ask TM-align to start with an alignment, specified in fasta\n"
"          file 'align.txt'\n"
"\n"
"    -I    Ask TM-align to stick the alignment to 'align.txt'\n"
"\n"
"    -m    Output TM-align rotation matrix\n"
"\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n"
"\n"
"    -o    output the superposition of chain1 to TM.sup\n"
"          $ TMalign chain1 chain2 -o TM.sup\n"
"          To view superimposed full-atom structures:\n"
"          $ pymol TM.sup chain2\n"
"\n"
"    -v    print the version of TM-align\n"
"\n"
"    -h    print the full help message, including options not available\n"
"          in standard TM-align program\n"
"\n"
"    (Options -u, -a, -d, -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    TMalign PDB1.pdb PDB2.pdb\n"
"    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0\n"
"    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup\n"
"    TMalign PDB1.pdb PDB2.pdb -i align.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -m matrix.txt\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = "";
    string yname       = "";
    string fname_super = ""; // file name for superposed structure
    string fname_lign  = ""; // file name for user alignment
    string fname_matrix= ""; // file name for output matrix
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool A_opt = false; // marker for whether structure A is specified
    bool B_opt = false; // marker for whether structure B is specified
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool i_opt = false; // flag for -i, with user given initial alignment
    bool I_opt = false; // flag for -I, stick to user given alignment
    bool o_opt = false; // flag for -o, output superposed structure
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    int    infmt1_opt=0;     // PDB format for chain_1
    int    infmt2_opt=0;     // PDB format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    int    byresi_opt=0;     // set -byresi to 0
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set

    int nameIdx = 0;
    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-o") && i < (argc-1) )
        {
            fname_super = argv[i + 1];     o_opt = true; i++;
        }
        else if ( (!strcmp(argv[i],"-u") || 
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
        else if ( !strcmp(argv[i],"-v") )
        {
            v_opt = true;
        }
        else if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if ( !strcmp(argv[i],"-i") && i < (argc-1) )
        {
            fname_lign = argv[i + 1];      i_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-m") && i < (argc-1) )
        {
            fname_matrix = argv[i + 1];    m_opt = true; i++;
        }// get filename for rotation matrix
        else if (!strcmp(argv[i], "-I") && i < (argc-1) )
        {
            fname_lign = argv[i + 1];      I_opt = true; i++;
        }
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
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
        else if ( !strcmp(argv[i],"-byresi") && i < (argc-1) )
        {
            byresi_opt=atoi(argv[i + 1]); i++;
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
    {

        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
    }

    if (!A_opt && dir_opt.size()==0) 
        PrintErrorAndQuit("Please provide structure A");
    if (!B_opt) 
        PrintErrorAndQuit("Please provide structure B");

    if (suffix_opt.size() && dir1_opt.size()==0 && dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir1 or -dir2 is set");
    if ((dir_opt.size() || dir1_opt.size() || dir2_opt.size()))
    {
        if (m_opt || o_opt)
            PrintErrorAndQuit("-m or -o cannot be set with -dir, -dir1 or -dir2");
        else if (dir_opt.size() && (dir1_opt.size() || dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! atom name must have 4 characters, including space.");

    if (i_opt && I_opt)
        PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
    if (u_opt && Lnorm_ass<=0)
        PrintErrorAndQuit("Wrong value for option -u!  It should be >0");
    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d!  It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (byresi_opt)
    {
        if (i_opt || I_opt)
            PrintErrorAndQuit("-byresi 1 or 2 cannot be used with -i or -I");
        if (byresi_opt<0 || byresi_opt>3)
            PrintErrorAndQuit("-byresi can only be 0, 1, 2 or 3");
        if (split_opt)
            PrintErrorAndQuit("-byresi >0 cannot be used with -split >0");
    }
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    

    /* read initial alignment file from 'align.txt' */
    string basename = string(argv[0]);
    int idx = basename.find_last_of("\\");
    basename = basename.substr(0, idx + 1);
    if (i_opt || I_opt)// Ask TM-align to start with an alignment,specified in fasta file 'align.txt'
    {
        if (fname_lign == "")
            PrintErrorAndQuit("Please provide a file name for option -i!");
        // open alignment file
        int n_p = 0;// number of structures in alignment file
        string line;

        string fullpath = basename + fname_lign;
        ifstream fileIn(fullpath.c_str());
        if (fileIn.is_open())
        {
            while (fileIn.good())
            {
                getline(fileIn, line);
                if (line.compare(0, 1, ">") == 0)// Flag for a new structure
                {
                    if (n_p >= 2) break;
                    sequence.push_back("");
                    n_p++;
                }
                else if (n_p > 0 && line!="") sequence.back()+=line;
            }
            fileIn.close();
        }
        else
            PrintErrorAndQuit("ERROR! Alignment file does not exist.");

        if (n_p < 2)
            PrintErrorAndQuit("ERROR: Fasta format is wrong, two proteins should be included.");
        if (sequence[0].size() != sequence[1].size())
            PrintErrorAndQuit("ERROR! FASTA file is wrong. The length in alignment should be equal respectively to the two aligned proteins.");
        if (I_opt)
        {
            int aligned_resNum=0;
            for (int i=0;i<sequence[0].size();i++) 
                aligned_resNum+=(sequence[0][i]!='-' && sequence[1][i]!='-');
            if (aligned_resNum<3)
                PrintErrorAndQuit("ERROR! Superposition is undefined for <3 aligned residues.");
        }
    }

    if (byresi_opt) I_opt=true;

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

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
    int    *secx, *secy;       // for the secondary structure 
    int    *xresno, *yresno;   // residue number for fragment gapless threading
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    /* loop over file names */
    for (int i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            resi_vec1, byresi_opt, ter_opt, infmt1_opt, atom_opt, split_opt);
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
            else if (xlen<=5)
            {
                cerr<<"Sequence is too short <=5!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new int[xlen];
            xresno = new int[xlen];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, xresno);
            make_sec(xa, xlen, secx); // secondary structure assignment

            for (int j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2,
                        chainID_list2, resi_vec2, byresi_opt, ter_opt,
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
                    else if (ylen<=5)
                    {
                        cerr<<"Sequence is too short <=5!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    yresno = new int[ylen];
                    secy = new int[ylen];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy, yresno);
                    make_sec(ya, ylen, secy);

                    if (byresi_opt)
                    {
                        sequence.clear();
                        sequence.push_back("");
                        sequence.push_back("");
                        extract_aln_from_resi(sequence,seqx,seqy,
                            resi_vec1,resi_vec2,byresi_opt);
                    }

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
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
                    TMalign_main(
                        xa, ya, xresno, yresno, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, I_opt, a_opt, u_opt, d_opt, fast_opt);

                    /* print result */
                    if (outfmt_opt==0) print_version();
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
                        (m_opt?fname_matrix+chainID_list1[chain_i]:"").c_str(),
                        outfmt_opt, ter_opt, 
                        (o_opt?fname_super+chainID_list1[chain_i]:"").c_str(),
                        i_opt, I_opt, a_opt, u_opt, d_opt);

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    delete [] secy;
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
            delete [] secx;
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
    sequence.clear();

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("Total CPU time is %5.2f seconds\n", diff);
    return 0;
}

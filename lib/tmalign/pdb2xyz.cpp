#include "basic_fun.h"

using namespace std;

void print_help()
{
    cout <<
"Converting PDB file(s) into xyz format.\n"
"\n"
"Usage: pdb2xyz pdb.pdb > ca.xyz\n"
"\n"
"    -dir     Convert all chains listed by 'chain_list' under 'chain_folder'.\n"
"             Note that the slash is necessary.\n"
"             $ pdb2xyz -dir chain_folder/ chain_list\n"
"\n"
"    -suffix  (Only when -dir is set, default is empty)\n"
"             add file name suffix to files listed by chain_list\n"
"\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: (default) TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: end of file\n"
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: (default) treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: treat each chain as a seperate chain (-ter should be <=1)\n"
    <<endl;
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    /**********************/
    /*    get argument    */
    /**********************/
    string xname     = "";
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    vector<string> chain_list; // only when -dir1 is set

    int nameIdx = 0;
    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
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
        else if ( !strcmp(argv[i],"-suffix") && i < (argc-1) )
        {
            suffix_opt=argv[i + 1]; i++;
        }
        else xname=argv[i];
    }

    if(xname.size()==0||xname=="-h") print_help();

    if (suffix_opt.size() && dir_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir is set");
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! atom name must have 4 characters, including space.");
    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");

    /* parse file list */
    if (dir_opt.size()==0)
        chain_list.push_back(xname);
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
            chain_list.push_back(dir_opt+Trim(line)+suffix_opt);
        }
        fp.close();
        line.clear();
    }

    /* declare previously global variables */
    vector<vector<string> >PDB_lines; // text of chain
    vector<string> chainID_list;      // list of chainID1
    vector<string> resi_vec;          // residue index for chain
    int    i;                         // file index
    int    l;                         // residue index
    int    chain_i;                   // chain index
    int    xlen;                      // chain length
    int    xchainnum;                 // number of chains in a PDB file

    /* loop over file names */
    for (i=0;i<chain_list.size();i++)
    {
        xname=chain_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines, chainID_list,
            resi_vec, 0, ter_opt, 0, atom_opt, split_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines[chain_i].size();
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            
            cout<<xlen<<'\n'<<xname.substr(dir_opt.size(),
                xname.size()-dir_opt.size()-suffix_opt.size())
                <<chainID_list[chain_i];
            for (l=0;l<PDB_lines[chain_i].size();l++)
            {
                cout<<'\n'<<AAmap(PDB_lines[chain_i][l].substr(17,3))<<' '
                    <<PDB_lines[chain_i][l].substr(30,8)<<' '
                    <<PDB_lines[chain_i][l].substr(38,8)<<' '
                    <<PDB_lines[chain_i][l].substr(46,8);
            }
            cout<<endl;

            PDB_lines[chain_i].clear();
        } // chain_i
        xname.clear();
        PDB_lines.clear();
        resi_vec.clear();
    } // i
    chain_list.clear();
    return 0;
}

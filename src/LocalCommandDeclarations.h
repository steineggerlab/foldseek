#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int structcreatedb(int argc, const char **argv, const Command& command);
extern int structuresearch(int argc, const char** argv, const Command &command);
extern int structureindex(int argc, const char** argv, const Command &command);
extern int structurecluster(int argc, const char** argv, const Command &command);
extern int easystructuresearch(int argc, const char** argv, const Command &command);
extern int easystructurecluster(int argc, const char** argv, const Command &command);
extern int tmalign(int argc, const char** argv, const Command &command);
extern int aln2tmscore(int argc, const char** argv, const Command &command);
extern int structurealign(int argc, const char** argv, const Command &command);
extern int samplemulambda(int argc, const char** argv, const Command &command);
extern int structureconvertalis(int argc, const char** argv, const Command &command);
extern int structureto3didescriptor(int argc, const char** argv, const Command &command);
extern int structurerbh(int argc, const char** argv, const Command &command);
extern int structureeasyrbh(int argc, const char** argv, const Command &command);
extern int structureungappedalign(int argc, const char** argv, const Command &command);
extern int convert2pdb(int argc, const char** argv, const Command &command);
extern int compressca(int argc, const char** argv, const Command &command);
extern int scoremultimer(int argc, const char **argv, const Command& command);
extern int filtermultimer(int argc, const char **argv, const Command& command);
extern int easymultimercluster(int argc, const char** argv, const Command &command);
extern int multimercluster(int argc, const char** argv, const Command &command);
extern int easymultimersearch(int argc, const char **argv, const Command &command);
extern int createmultimerreport(int argc, const char **argv, const Command &command);
extern int expandmultimer(int argc, const char **argv, const Command &command);
extern int multimersearch(int argc, const char **argv, const Command &command);
extern int makepaddeddb(int argc, const char **argv, const Command& command);
extern int result2structprofile(int argc, const char **argv, const Command& command);
extern int createstructsubdb(int argc, const char **argv, const Command& command);
extern int fwbw(int argc, const char **argv, const Command& command);
extern int lolalign(int argc, const char **argv, const Command& command);
#endif

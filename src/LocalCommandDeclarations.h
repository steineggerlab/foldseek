#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int strucclust(int argc, const char** argv, const Command &command);
extern int structuresearch(int argc, const char** argv, const Command &command);
extern int structcreatedb(int argc, const char** argv, const Command &command);
extern int structkmermatcher(int argc, const char** argv, const Command &command);
extern int tmalign(int argc, const char** argv, const Command &command);
extern int aln2tmscore(int argc, const char** argv, const Command &command);
#endif

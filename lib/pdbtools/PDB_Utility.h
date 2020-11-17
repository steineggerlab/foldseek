#pragma once

//-----extern PDB_switch ------// 
// transer CHAR to 63
extern int CHAIN_to_INT62(char code);
// transer 63 to CHAR
extern int INT62_to_CHAIN(int code);
// transer 26 to 20
extern int AA26_to_AA20(int aa26);
// transer 20 to 26
extern int AA20_to_AA26(int aa20);
//get the side chain size of given amino acid
extern int AA26_sidechain_size(int aa26);
// three-digit to one-digit
extern char Three2One_III(const char *input);
// one-digit to three-digit
extern const char* One2Three_III(char c);

//-----extern SideChain_switch ------// 
// hashing
extern int PDB_atom_name_hashing(const char *atom_name);
// for backbone, atom_name to index
extern int backbone_atom_name_encode(const char *atom_name);
// for backbone, index to atom_name
extern const char *backbone_atom_name_decode(int index);
// for sidechain, atom_name to index
extern int sidechain_atom_name_encode(const char *atom_name, char amino);
// for sidechain, index to atom_name
extern const char *sidechain_atom_name_decode(int index, char amino);

#include "PDB_Utility.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

//------ Digit_Switch -----//
//DIGIT_trans 
int CHAIN_to_INT62(char code) 
{ 
	switch(code) 
	{ 
		case 'A': return 0; 
		case 'B': return 1; 
		case 'C': return 2; 
		case 'D': return 3; 
		case 'E': return 4; 
		case 'F': return 5; 
		case 'G': return 6; 
		case 'H': return 7; 
		case 'I': return 8; 
		case 'J': return 9; 
		case 'K': return 10; 
		case 'L': return 11; 
		case 'M': return 12; 
		case 'N': return 13; 
		case 'O': return 14; 
		case 'P': return 15; 
		case 'Q': return 16; 
		case 'R': return 17; 
		case 'S': return 18; 
		case 'T': return 19; 
		case 'U': return 20; 
		case 'V': return 21; 
		case 'W': return 22; 
		case 'X': return 23; 
		case 'Y': return 24; 
		case 'Z': return 25; 
		case 'a': return 26; 
		case 'b': return 27; 
		case 'c': return 28; 
		case 'd': return 29; 
		case 'e': return 30; 
		case 'f': return 31; 
		case 'g': return 32; 
		case 'h': return 33; 
		case 'i': return 34; 
		case 'j': return 35; 
		case 'k': return 36; 
		case 'l': return 37; 
		case 'm': return 38; 
		case 'n': return 39; 
		case 'o': return 40; 
		case 'p': return 41; 
		case 'q': return 42; 
		case 'r': return 43; 
		case 's': return 44; 
		case 't': return 45; 
		case 'u': return 46; 
		case 'v': return 47; 
		case 'w': return 48; 
		case 'x': return 49; 
		case 'y': return 50; 
		case 'z': return 51; 
		case '0': return 52; 
		case '1': return 53; 
		case '2': return 54; 
		case '3': return 55; 
		case '4': return 56; 
		case '5': return 57; 
		case '6': return 58; 
		case '7': return 59; 
		case '8': return 60; 
		case '9': return 61;
		case ' ': return 62; 	
		default:return -1; 
	} 
}
int INT62_to_CHAIN(int code) 
{ 
	switch(code) 
	{ 
		case 0: return 'A'; 
		case 1: return 'B'; 
		case 2: return 'C'; 
		case 3: return 'D'; 
		case 4: return 'E'; 
		case 5: return 'F'; 
		case 6: return 'G'; 
		case 7: return 'H'; 
		case 8: return 'I'; 
		case 9: return 'J'; 
		case 10:return 'K';  
		case 11:return 'L';  
		case 12:return 'M';  
		case 13:return 'N';  
		case 14:return 'O';  
		case 15:return 'P';  
		case 16:return 'Q';  
		case 17:return 'R';  
		case 18:return 'S';  
		case 19:return 'T';  
		case 20:return 'U';  
		case 21:return 'V';  
		case 22:return 'W';  
		case 23:return 'X';  
		case 24:return 'Y';  
		case 25:return 'Z';  
		case 26:return 'a';  
		case 27:return 'b';  
		case 28:return 'c';  
		case 29:return 'd';  
		case 30:return 'e';  
		case 31:return 'f';  
		case 32:return 'g';  
		case 33:return 'h';  
		case 34:return 'i';  
		case 35:return 'j';  
		case 36:return 'k';  
		case 37:return 'l';  
		case 38:return 'm';  
		case 39:return 'n';  
		case 40:return 'o';  
		case 41:return 'p';  
		case 42:return 'q';  
		case 43:return 'r';  
		case 44:return 's';  
		case 45:return 't';  
		case 46:return 'u';  
		case 47:return 'v';  
		case 48:return 'w';  
		case 49:return 'x';  
		case 50:return 'y';  
		case 51:return 'z';  
		case 52:return '0';  
		case 53:return '1';  
		case 54:return '2';  
		case 55:return '3';  
		case 56:return '4';  
		case 57:return '5';  
		case 58:return '6';  
		case 59:return '7';  
		case 60:return '8';  
		case 61:return '9'; 
		case 62:return ' ';  	
		default:return -1; 
	} 
}

//------- Switch -------// 
//AMI_trans 
int AA20_to_AA26(int code) 
{ 
	switch(code) 
	{ 
		case 0: return  0;//A 
		case 1: return  2;//C 
		case 2: return  3;//D 
		case 3: return  4;//E 
		case 4: return  5;//F 
		case 5: return  6;//G 
		case 6: return  7;//H 
		case 7: return  8;//I 
		case 8: return 10;//K 
		case 9: return 11;//L 
		case 10:return 12;//M 
		case 11:return 13;//N 
		case 12:return 15;//P 
		case 13:return 16;//Q 
		case 14:return 17;//R 
		case 15:return 18;//S 
		case 16:return 19;//T 
		case 17:return 21;//V 
		case 18:return 22;//W 
		case 19:return 24;//Y 
		case 20:return 25;//Z 
		default:return -1; 
	} 
}
int AA26_to_AA20(int amino)
{
	switch (amino)
	{
		case 0:return 0;	//A 
		case 1:return 20;	//B 
		case 2:return 1;	//C 
		case 3:return 2;	//D 
		case 4:return 3;	//E 
		case 5:return 4;	//F 
		case 6:return 5;	//G 
		case 7:return 6;	//H 
		case 8:return 7;	//I 
		case 9:return 20;	//J 
		case 10:return 8;	//K 
		case 11:return 9;	//L 
		case 12:return 10;	//M 
		case 13:return 11;	//N 
		case 14:return 20;	//O 
		case 15:return 12;	//P 
		case 16:return 13;	//Q 
		case 17:return 14;	//R 
		case 18:return 15;	//S 
		case 19:return 16;	//T 
		case 20:return 20;	//U 
		case 21:return 17;	//V 
		case 22:return 18;	//W 
		case 23:return 20;	//X 
		case 24:return 19;	//Y 
		case 25:return 20;	//Z 
		default:return -1;
	}
}
//SideChain_trans 
int AA26_sidechain_size(int amino)
{
	switch (amino)
	{
		case 0:return 1;	//A 
		case 1:return 1;	//B 
		case 2:return 2;	//C 
		case 3:return 4;	//D 
		case 4:return 5;	//E 
		case 5:return 7;	//F 
		case 6:return 1;	//G -> note: pseudo CB(=CA)
		case 7:return 6;	//H 
		case 8:return 4;	//I 
		case 9:return 1;	//J 
		case 10:return 5;	//K 
		case 11:return 4;	//L 
		case 12:return 4;	//M 
		case 13:return 4;	//N 
		case 14:return 1;	//O 
		case 15:return 3;	//P 
		case 16:return 5;	//Q 
		case 17:return 7;	//R 
		case 18:return 2;	//S 
		case 19:return 3;	//T 
		case 20:return 1;	//U 
		case 21:return 3;	//V 
		case 22:return 10;	//W 
		case 23:return 1;	//X 
		case 24:return 8;	//Y 
		case 25:return 1;	//Z 
		default:return -1;
	}
}

//---------- One_To_Three -----------//
char Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(26.0,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}
const char* One2Three_III(char c)
{
	//switch
	switch(c)
	{
		case 'A':return "ALA";
		case 'R':return "ARG";
		case 'N':return "ASN";
		case 'D':return "ASP";
		case 'C':return "CYS";
		case 'Q':return "GLN";
		case 'E':return "GLU";
		case 'G':return "GLY";
		case 'H':return "HIS";
		case 'I':return "ILE";
		case 'L':return "LEU";
		case 'K':return "LYS";
		case 'M':return "MET";
		case 'F':return "PHE";
		case 'P':return "PRO";
		case 'S':return "SER";
		case 'T':return "THR";
		case 'W':return "TRP";
		case 'Y':return "TYR";
		case 'V':return "VAL";
		default:return "UNK";
	}
}

//----------- SideChain_Related ------------//
//function (input three-digit atom)
int PDB_atom_name_hashing(const char *atom)
{
	int i;
	int result=0;
	int pos;
	int len=(int)strlen(atom);
	int rellen;
	// modify by wangsheng, 2011.04.30
	rellen=len<3?len:3;
	// modify by shaomingfu, 2011.04.12
	if(len!=3)len=3;
	for(i=0;i<rellen;i++)
	{
		if(atom[i]==' ')pos=0;
		else if(atom[i]>='A'&&atom[i]<='Z')pos=atom[i]-'A';
		else if(atom[i]>='0'&&atom[i]<='9')pos=atom[i]-'0';
		else pos=0;
		result+=(int)(pow(26.0,1.0*(len-i-1)))*pos;
	}
	return result;
}
//backbone 
int backbone_atom_name_encode(const char *atom)
{
	int pos=PDB_atom_name_hashing(atom);
	switch(pos)
	{
		case 8788:return 0;
		case 1352:
		{
			int len=(int)strlen(atom);
			if(len<2)return -1;
			char add=atom[1]; //second atom
			switch(add)
			{
				case 'A':return 1;
				case ' ':return 2;
				default:return -1;
			}
		}
		case 9464:return 3;
		default:return -1;
	}
}
const char *backbone_atom_name_decode(int pos)
{
	switch (pos)
	{
		case 0:return "N  ";
		case 1:return "CA ";
		case 2:return "C  ";
		case 3:return "O  ";
		default:return 0;
	}
}
//sidechain (ARNDCQEGHILKMFPSTWYVZ) 
int sidechain_atom_name_encode(const char *atom, char amino)
{
	int pos = PDB_atom_name_hashing(atom);
	switch(amino)
	{
		case 'A':
		{
			switch(pos)
			{
				case 1378:return 0;
				default:return -1;
			}
		}
		case 'R':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 8892:return 3;
				case 2002:return 4;
				case 8971:return 5;
				case 8972:return 6;
				default:return -1;
			}
		}
		case 'N':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 9543:return 2;
				case 8868:return 3;
				default:return -1;
			}
		}
		case 'D':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 9543:return 2;
				case 9544:return 3;
				default:return -1;
			}
		}
		case 'C':
		{
			switch(pos)
			{
				case 1378 :return 0;
				case 12324:return 1;
				default:return -1;
			}
		}
		case 'Q':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 9569:return 3;
				case 8894:return 4;
				default:return -1;
			}
		}
		case 'E':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 9569:return 3;
				case 9570:return 4;
				default:return -1;
			}
		}
		case 'G':
		{
			switch(pos)
			{
				case 1378:return 0; //note: pseudo CB(=CA)
				default:return -1;
			}
		}
		case 'H':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 8867:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 8894:return 5;
				default:return -1;
			}
		}
		case 'I':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1509:return 1;
				case 1510:return 2;
				case 1431:return 3;
				default:return -1;
			}
		}
		case 'L':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				default:return -1;
			}
		}
		case 'K':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 1456:return 3;
				case 9438:return 4;
				default:return -1;
			}
		}
		case 'M':
		{
			switch(pos)
			{
				case 1378 :return 0;
				case 1508 :return 1;
				case 12246:return 2;
				case 1456 :return 3;
				default:return -1;
			}
		}
		case 'F':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 1458:return 5;
				case 2002:return 6;
				default:return -1;
			}
		}
		case 'P':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				default:return -1;
			}
		}
		case 'S':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 9620:return 1;
				default:return -1;
			}
		}
		case 'T':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 9621:return 1;
				case 1510:return 2;
				default:return -1;
			}
		}
		case 'W':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 8893:return 4;
				case 1458:return 5;
				case 1459:return 6;
				case 2004:return 7;
				case 2005:return 8;
				case 1536:return 9;
				default:return -1;
			}
		}
		case 'Y':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 1458:return 5;
				case 2002:return 6;
				case 9646:return 7;
				default:return -1;
			}
		}
		case 'V':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1509:return 1;
				case 1510:return 2;
				default:return -1;
			}
		}
		case 'X':  //only consider CB!! (alanine model)
		{
			switch(pos)
			{
				case 1378:return 0;
				default:return -1;
			}
		}
		default:return -1;
	}
}
const char *sidechain_atom_name_decode(int pos, char amino)
{
	switch(amino)
	{
		case 'A':
		{
			switch(pos)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		case 'R':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "NE ";
				case 4:return "CZ ";
				case 5:return "NH1";
				case 6:return "NH2";
				default:return 0;
			}
		}
		case 'N':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "OD1";
				case 3:return "ND2";
				default:return 0;
			}
		}
		case 'D':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "OD1";
				case 3:return "OD2";
				default:return 0;
			}
		}
		case 'C':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "SG ";
				default:return 0;
			}
		}
		case 'Q':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "OE1";
				case 4:return "NE2";
				default:return 0;
			}
		}
		case 'E':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "OE1";
				case 4:return "OE2";
				default:return 0;
			}
		}
		case 'G':
		{
			switch(pos)   //note: pseudo CB (=CA)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		case 'H':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "ND1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "NE2";
				default:return 0;
			}
		}
		case 'I':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG1";
				case 2:return "CG2";
				case 3:return "CD1";
				default:return 0;
			}
		}
		case 'L':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				default:return 0;
			}
		}
		case 'K':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "CE ";
				case 4:return "NZ ";
				default:return 0;
			}
		}
		case 'M':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "SD ";
				case 3:return "CE ";
				default:return 0;
			}
		}
		case 'F':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "CE2";
				case 6:return "CZ ";
				default:return 0;
			}
		}
		case 'P':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				default:return 0;
			}
		}
		case 'S':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "OG ";
				default:return 0;
			}
		}
		case 'T':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "OG1";
				case 2:return "CG2";
				default:return 0;
			}
		}
		case 'W':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "NE1";
				case 5:return "CE2";
				case 6:return "CE3";
				case 7:return "CZ2";
				case 8:return "CZ3";
				case 9:return "CH2";
				default:return 0;
			}
		}
		case 'Y':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "CE2";
				case 6:return "CZ ";
				case 7:return "OH ";
				default:return 0;
			}
		}
		case 'V':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG1";
				case 2:return "CG2";
				default:return 0;
			}
		}
		case 'X':  //only consider CB!! (alanine model)
		{
			switch(pos)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		default:return 0;
	}
}

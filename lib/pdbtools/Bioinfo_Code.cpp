#include "Bioinfo_Code.h"


//--------------constructor--------------------//
Bioinfo_Code::Bioinfo_Code(int version)
{
	//init
	NewArray2D(&Gen_CLESUM,3*CLE_NUM,3*CLE_NUM);
	Gen_CLESUM_Create(version);
	//parameter
	AMI_Weight=5.0;
	CLE_Weight=0.5;
}
Bioinfo_Code::~Bioinfo_Code(void)
{
	//dele
	DeleteArray2D(&Gen_CLESUM,3*CLE_NUM);
}


//===================================================================================================//
//=================//
//--Ori_CLESUM----//
//===============//
//origin_CLESUM(from FSSP)
int Ori_CLESUM[CLE_NUM][CLE_NUM]={
{ 73, 20, 13, -17, -25, -20, -6, -45, -31,-23,-19,-11, -2, 10, 25, 35, 16,0}, //A 
{ 20, 51,  7,  13,  15,   7, 13, -96, -74,-57,-50,-12,-13,-11,-12, 42, 12,0}, //B 
{ 13,  7, 53,  21,   3,  20, -4, -77, -56,-43,-33,  0,-12, -5,  3,  4, 29,0}, //C 
{-17, 13, 21,  52,  22,  22,-31,-124,-105,-88,-81,-22,-49,-44,-42,-10, 14,0}, //D 
{-25, 15,  3,  22,  36,  26,-22,-127,-108,-93,-84,-21,-47,-43,-48, -5, -6,0}, //E 
{-20,  7, 20,  22,  26,  50, -5,-107, -88,-73,-69,-16,-33,-32,-30,  0,  3,0}, //F 
{ -6, 13, -4, -31, -22,  -5, 69, -51, -34,-21,-13, 29, 21, -8, -1,  5,  8,0}, //G 
{-45,-96,-77,-124,-127,-107,-51,  23,  18, 13,  5,-62, -4,-34,-55,-60,-87,0}, //H 
{-31,-74,-56,-105,-108, -88,-34,  18,  23, 16, 21,-41,  1,-11,-34,-49,-62,0}, //I 
{-23,-57,-43, -88, -93, -73,-21,  13,  16, 37, 13,-32, 16, -2,-24,-34,-44,0}, //J 
{-19,-50,-33, -81, -84, -69,-13,   5,  21, 13, 49, -1, 12, 28,  5,-36,-24,0}, //K 
{-11,-12,  0, -22, -21, -16, 29, -62, -41,-32, -1, 74,  5,  8, -4,-12, 26,0}, //L 
{ -2,-13,-12, -49, -47, -33, 21,  -4,   1, 16, 12,  5, 61,  7,  5,  8, -7,0}, //M 
{ 10,-11, -5, -44, -43, -32, -8, -34, -11, -2, 28,  8,  7, 90, 15, -3, 32,0}, //N 
{ 25,-12,  3, -42, -48, -30, -1, -55, -34,-24,  5, -4,  5, 15,104,  4,-13,0}, //O 
{ 35, 42,  4, -10,  -5,   0,  5, -60, -49,-34,-36,-12,  8, -3,  4, 66,  7,0}, //P 
{ 16, 12, 29,  14,  -6,   3,  8, -87, -62,-44,-24, 26, -7, 32,-13,  7, 90,0}, //Q 
{  0,  0,  0,   0,   0,   0,  0,   0,   0,  0,  0,  0,  0,  0,  0,  0,  0,0}};//R  
// A   B   C    D    E    F   G    H    I   J   K   L   M   N   O   P   Q R


//===================================================================================================//
//=================//
//--Ori_BLOSUM----//(BLOSUM_62)
//===============//
int Ori_BLOSUM[AMI_NUM][AMI_NUM]={
{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5 },  //A
{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5 },  //R
{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5 },  //N
{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5 },  //D
{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5 },  //C
{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5 },  //Q
{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5 },  //E
{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5 },  //H
{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5 },  //I
{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5 },  //L
{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5 },  //K
{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5 },  //M
{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5 },  //F
{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5 },  //P
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5 },  //T
{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5 },  //W
{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5 },  //Y
{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z

//BLOSUM_Mapping//--------------ARNDCQEGHILKMFPSTWYVZ
int Blo_AA_Map[AMI_NUM]=
{ 0,19, 4, 3, 6, 13,7, 8, 9, 17,11,10,12,2, 18,14,5, 1, 15,16,20};
//A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R   S  T  Z
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20

//Ori_Mapping//-----------------AVCDEFGHIWKLMNYPQRSTZ
int Ori_AA_Map[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

//Ori_Anti_Mapping/
int Ori_Anti_Map[AMI_NUM]=
{ 0,21,2,3,4,5,6,7,8,22,10,11,12,13,24,15,16,17,18,19,25};
//A  V C D E F G H I  W  K  L  M  N  Y  P  Q  R  S  T  Z
//0  1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 17 18 19 20

const char Ami_Clus[AMI_NUM+1]="000110110010010111112";
//------------------------------AVCDEFGHIWKLMNYPQRSTZ



//===================================================================================================//
//=================//
//--Gen_CLESUM----//
//===============//
//>[0][0]->Cluster
int Ori_Gen_CLESUM_00[CLE_NUM][CLE_NUM]={
{  85,   21,   17,  -22,  -23,  -24,  -10,  -36,  -26,  -19,  -25,  -24,   -6,   13,   32,   42,   11,  0}, 
{  21,   59,   10,   14,   24,   10,   31,  -99,  -82,  -64,  -61,   -3,  -12,   -7,   13,   45,   15,  0}, 
{  17,   10,   62,   26,    9,   29,   -7,  -87,  -67,  -49,  -45,   -8,  -24,  -17,   14,    7,   36,  0}, 
{ -22,   14,   26,   63,   30,   28,  -36, -138, -111,  -98,  -95,  -20,  -66,  -67,  -31,  -10,   33,  0}, 
{ -23,   24,    9,   30,   44,   34,  -22, -138, -126, -110,  -98,  -24,  -58,  -56,  -22,    3,   -5,  0}, 
{ -24,   10,   29,   28,   34,   61,   -4, -131, -105,  -92, -106,  -24,  -45,  -46,    5,    3,    3,  0}, 
{ -10,   31,   -7,  -36,  -22,   -4,   80,  -48,  -35,  -22,  -10,   40,   25,  -11,   13,    4,    2,  0}, 
{ -36,  -99,  -87, -138, -138, -131,  -48,   30,   26,   18,   14,  -58,    7,  -16,  -37,  -57, -110,  0}, 
{ -26,  -82,  -67, -111, -126, -105,  -35,   26,   30,   20,   28,  -37,    9,    3,  -18,  -47,  -67,  0}, 
{ -19,  -64,  -49,  -98, -110,  -92,  -22,   18,   20,   40,   18,  -33,   24,   -2,   -9,  -31,  -41,  0}, 
{ -25,  -61,  -45,  -95,  -98, -106,  -10,   14,   28,   18,   59,   10,   16,   40,  -14,  -39,  -19,  0}, 
{ -24,   -3,   -8,  -20,  -24,  -24,   40,  -58,  -37,  -33,   10,   86,    9,    8,   -8,  -21,   30,  0}, 
{  -6,  -12,  -24,  -66,  -58,  -45,   25,    7,    9,   24,   16,    9,   71,   10,   11,   13,   -5,  0}, 
{  13,   -7,  -17,  -67,  -56,  -46,  -11,  -16,    3,   -2,   40,    8,   10,  108,   25,   -1,   27,  0}, 
{  32,   13,   14,  -31,  -22,    5,   13,  -37,  -18,   -9,  -14,   -8,   11,   25,  121,   36,   29,  0}, 
{  42,   45,    7,  -10,    3,    3,    4,  -57,  -47,  -31,  -39,  -21,   13,   -1,   36,   76,   13,  0}, 
{  11,   15,   36,   33,   -5,    3,    2, -110,  -67,  -41,  -19,   30,   -5,   27,   29,   13,  109,  0}, 
{   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0}}; 

//>[0][1]->Cluster                                                                                            
int Ori_Gen_CLESUM_01[CLE_NUM][CLE_NUM]={                                                   
{  66,   18,    7,  -17,  -21,  -18,   -8,  -46,  -30,  -28,  -25,  -10,   -4,    3,    0,   32,   12,  0}, 
{   2,   34,   -3,   11,   12,    4,   -5, -112,  -80,  -65,  -48,  -13,  -17,   -9,  -35,   21,   16,  0}, 
{   2,    6,   42,   23,    8,   21,   -5,  -87,  -65,  -50,  -40,   -4,  -16,  -13,  -22,   -1,   15,  0}, 
{ -31,   -1,    2,   34,   14,   14,  -37, -139, -121, -101,  -93,  -29,  -54,  -59,  -57,  -27,   -5,  0}, 
{ -37,   -2,  -14,    9,   25,   14,  -31, -140, -119, -103,  -94,  -26,  -54,  -46,  -69,  -22,  -14,  0}, 
{ -30,   -3,    0,    9,   18,   33,  -13, -126, -108,  -84,  -77,  -24,  -41,  -41,  -55,  -16,  -11,  0}, 
{  -8,    5,   -8,  -32,  -24,   -9,   52,  -60,  -40,  -31,  -19,   17,    9,   -8,   -4,    4,    5,  0}, 
{ -49,  -95,  -74, -114, -116,  -95,  -54,   14,    9,    3,   -2,  -65,  -13,  -48,  -62,  -62,  -94,  0}, 
{ -37,  -75,  -57, -100,  -95,  -79,  -35,   10,   13,    6,   12,  -47,   -5,  -25,  -43,  -50,  -66,  0}, 
{ -26,  -54,  -42,  -84,  -81,  -64,  -20,   10,   13,   27,   11,  -35,   11,  -12,  -35,  -34,  -50,  0}, 
{ -24,  -53,  -36,  -79,  -80,  -62,  -15,   -1,   13,    3,   36,   -8,    4,    4,    4,  -36,  -37,  0}, 
{ -19,  -19,   -4,  -25,  -23,  -15,   16,  -69,  -44,  -38,  -15,   65,   -3,  -10,   -7,  -19,   11,  0}, 
{  -7,  -21,  -12,  -45,  -43,  -30,   12,  -11,   -7,    6,    4,   -3,   47,   -1,   -3,    1,  -13,  0}, 
{   7,  -19,  -10,  -30,  -44,  -20,  -17,  -25,   -9,   -3,   23,   12,   12,   76,    9,   -6,   14,  0}, 
{  26,    4,    9,  -29,  -19,    2,   21,  -35,  -21,   -8,    8,   10,   13,   15,   67,   17,   -9,  0}, 
{  26,   35,    1,   -3,   -1,    5,    0,  -67,  -59,  -40,  -40,   -8,    4,   -6,  -12,   53,    5,  0}, 
{   6,   -1,   30,   20,   -6,    2,    6,  -86,  -64,  -47,  -37,   25,   -1,   16,  -21,    0,   70,  0}, 
{   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0}}; 

//>[1][0]->Cluster
int Ori_Gen_CLESUM_10[CLE_NUM][CLE_NUM]={
{  66,    2,    2,  -31,  -37,  -30,   -8,  -49,  -37,  -26,  -24,  -19,   -7,    7,   26,   26,    6,  0}, 
{  18,   34,    6,   -1,   -2,   -3,    5,  -95,  -75,  -54,  -53,  -19,  -21,  -19,    4,   35,   -1,  0}, 
{   7,   -3,   42,    2,  -14,    0,   -8,  -74,  -57,  -42,  -36,   -4,  -12,  -10,    9,    1,   30,  0}, 
{ -17,   11,   23,   34,    9,    9,  -32, -114, -100,  -84,  -79,  -25,  -45,  -30,  -29,   -3,   20,  0}, 
{ -21,   12,    8,   14,   25,   18,  -24, -116,  -95,  -81,  -80,  -23,  -43,  -44,  -19,   -1,   -6,  0}, 
{ -18,    4,   21,   14,   14,   33,   -9,  -95,  -79,  -64,  -62,  -15,  -30,  -20,    2,    5,    2,  0}, 
{  -8,   -5,   -5,  -37,  -31,  -13,   52,  -54,  -35,  -20,  -15,   16,   12,  -17,   21,    0,    6,  0}, 
{ -46, -112,  -87, -139, -140, -126,  -60,   14,   10,   10,   -1,  -69,  -11,  -25,  -35,  -67,  -86,  0}, 
{ -30,  -80,  -65, -121, -119, -108,  -40,    9,   13,   13,   13,  -44,   -7,   -9,  -21,  -59,  -64,  0}, 
{ -28,  -65,  -50, -101, -103,  -84,  -31,    3,    6,   27,    3,  -38,    6,   -3,   -8,  -40,  -47,  0}, 
{ -25,  -48,  -40,  -93,  -94,  -77,  -19,   -2,   12,   11,   36,  -15,    4,   23,    8,  -40,  -37,  0}, 
{ -10,  -13,   -4,  -29,  -26,  -24,   17,  -65,  -47,  -35,   -8,   65,   -3,   12,   10,   -8,   25,  0}, 
{  -4,  -17,  -16,  -54,  -54,  -41,    9,  -13,   -5,   11,    4,   -3,   47,   12,   13,    4,   -1,  0}, 
{   3,   -9,  -13,  -59,  -46,  -41,   -8,  -48,  -25,  -12,    4,  -10,   -1,   76,   15,   -6,   16,  0}, 
{   0,  -35,  -22,  -57,  -69,  -55,   -4,  -62,  -43,  -35,    4,   -7,   -3,    9,   67,  -12,  -21,  0}, 
{  32,   21,   -1,  -27,  -22,  -16,    4,  -62,  -50,  -34,  -36,  -19,    1,   -6,   17,   53,    0,  0}, 
{  12,   16,   15,   -5,  -14,  -11,    5,  -94,  -66,  -50,  -37,   11,  -13,   14,   -9,    5,   70,  0}, 
{   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0}}; 

//>[1][1]->Cluster
int Ori_Gen_CLESUM_11[CLE_NUM][CLE_NUM]={
{  75,   31,   20,   -4,  -16,  -11,   -1,  -46,  -30,  -21,  -12,   -3,    4,   14,   31,   38,   21,  0}, 
{  31,   64,   13,   25,   20,   14,    7,  -84,  -64,  -49,  -43,  -12,   -6,  -12,   -4,   55,   10,  0}, 
{  20,   13,   60,   33,   12,   28,    1,  -69,  -48,  -37,  -24,    7,   -4,    0,   10,    9,   33,  0}, 
{  -4,   25,   33,   61,   30,   30,  -17, -108,  -89,  -73,  -64,  -14,  -34,  -31,  -32,    0,   20,  0}, 
{ -16,   20,   12,   30,   45,   32,  -11, -113,  -93,  -79,  -66,  -11,  -32,  -35,  -37,    3,    3,  0}, 
{ -11,   14,   28,   30,   32,   58,    5,  -92,  -72,  -62,  -52,   -5,  -21,  -25,  -24,    7,   13,  0}, 
{  -1,    7,    1,  -17,  -11,    5,   80,  -45,  -28,  -13,   -8,   35,   33,   -5,   -4,    9,   14,  0}, 
{ -46,  -84,  -69, -108, -113,  -92,  -45,   30,   24,   19,    7,  -60,   -5,  -31,  -55,  -57,  -76,  0}, 
{ -30,  -64,  -48,  -89,  -93,  -72,  -28,   24,   31,   24,   27,  -38,    3,   -6,  -31,  -44,  -57,  0}, 
{ -21,  -49,  -37,  -73,  -79,  -62,  -13,   19,   24,   45,   19,  -27,   21,    3,  -22,  -32,  -41,  0}, 
{ -12,  -43,  -24,  -64,  -66,  -52,   -8,    7,   27,   19,   57,    5,   20,   39,    6,  -33,  -16,  0}, 
{  -3,  -12,    7,  -14,  -11,   -5,   35,  -60,  -38,  -27,    5,   76,   12,   15,   -4,   -8,   33,  0}, 
{   4,   -6,   -4,  -34,  -32,  -21,   33,   -5,    3,   21,   20,   12,   69,   11,    9,   13,   -4,  0}, 
{  14,  -12,    0,  -31,  -35,  -25,   -5,  -31,   -6,    3,   39,   15,   11,   94,   16,   -1,   39,  0}, 
{  31,   -4,   10,  -32,  -37,  -24,   -4,  -55,  -31,  -22,    6,   -4,    9,   16,  107,    7,  -14,  0}, 
{  38,   55,    9,    0,    3,    7,    9,  -57,  -44,  -32,  -33,   -8,   13,   -1,    7,   74,    8,  0}, 
{  21,   10,   33,   20,    3,   13,   14,  -76,  -57,  -41,  -16,   33,   -4,   39,  -14,    8,   95,  0}, 
{   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0}}; 




//=================Gen_CLESUM(create_related)=============//
void Bioinfo_Code::Gen_CLESUM_Create(int Version)
{
	if(Version==1)Gen_CLESUM_Create_Ori(); //Ori_Version//
	else Gen_CLESUM_Create_Glo();          //Glo_Version//
}
void Bioinfo_Code::Gen_CLESUM_Create_Ori(void)  //Ori_Version//
{
	int ii,jj;
	int i,j;
	for(ii=0;ii<3;ii++)for(jj=0;jj<3;jj++)
		for(i=0;i<CLE_NUM;i++)
			for(j=0;j<CLE_NUM;j++)Gen_CLESUM[ii*CLE_NUM+i][jj*CLE_NUM+j]=Ori_CLESUM[i][j];
}
void Bioinfo_Code::Gen_CLESUM_Create_Glo(void) //Glo_Version//
{
	int i,j;
	for(i=0;i<CLE_NUM;i++)for(j=0;j<CLE_NUM;j++)Gen_CLESUM[i][j]=Ori_Gen_CLESUM_00[i][j];
	for(i=0;i<CLE_NUM;i++)for(j=0;j<CLE_NUM;j++)Gen_CLESUM[i][j+CLE_NUM]=Ori_Gen_CLESUM_01[i][j];
	for(i=0;i<CLE_NUM;i++)for(j=0;j<CLE_NUM;j++)Gen_CLESUM[i+CLE_NUM][j]=Ori_Gen_CLESUM_10[i][j];
	for(i=0;i<CLE_NUM;i++)for(j=0;j<CLE_NUM;j++)Gen_CLESUM[i+CLE_NUM][j+CLE_NUM]=Ori_Gen_CLESUM_11[i][j];
	for(i=0;i<CLE_NUM;i++)
	{
		for(j=0;j<CLE_NUM;j++)
		{
			Gen_CLESUM[i][j+2*CLE_NUM]=Ori_CLESUM[i][j];
			Gen_CLESUM[i+CLE_NUM][j+2*CLE_NUM]=Ori_CLESUM[i][j];
			Gen_CLESUM[i+2*CLE_NUM][j]=Ori_CLESUM[i][j];
			Gen_CLESUM[i+2*CLE_NUM][j+CLE_NUM]=Ori_CLESUM[i][j];
			Gen_CLESUM[i+2*CLE_NUM][j+2*CLE_NUM]=Ori_CLESUM[i][j];
		}
	}	
}


//====================================//
//--PART_0:CLESUM(related)-----------//__PART_0=>START
//==================================//
//---- ami/cle transform -----//
void Bioinfo_Code::AMI_transform(const char *AMI,int *out)
{
	int i;
	int len=(int)strlen(AMI);
	for(i=0;i<len;i++)out[i]=Blo_AA_Map[Ori_AA_Map[AMI[i]-'A']];
}
void Bioinfo_Code::CLE_transform(const char *CLE,int *out)
{
	int i;
	int len=(int)strlen(CLE);
	for(i=0;i<len;i++)out[i]=CLE[i]-'A';
}
void Bioinfo_Code::AMI_CLE_transform_Ori(const char *AMI,const char *CLE,int *out)
{
#ifdef DEBUG
	if(strlen(AMI)!=strlen(CLE))
	{
		fprintf(stderr,"AMI_Len =!= CLE_Len =>BREAK!!\n");
		fprintf(stderr,"[%s][%s]\n",AMI,CLE);
		exit(-1);
	}
#endif
	int i;
	int len=(int)strlen(CLE);
	for(i=0;i<len;i++)out[i]=(CLE[i]-'A')+CLE_NUM*(Blo_AA_Map[Ori_AA_Map[AMI[i]-'A']]);
}
//---calculation_related-------// AMI+CLE->GLO
void Bioinfo_Code::AMI_CLE_transform(const char *AMI,const char *CLE,int *GEN)
{
#ifdef DEBUG
	if(strlen(AMI)!=strlen(CLE))
	{
		fprintf(stderr,"AMI_Len =!= CLE_Len =>BREAK!!\n");
		fprintf(stderr,"[%s][%s]\n",AMI,CLE);
		exit(-1);
	}
#endif
	int i;
	int len=(int)strlen(AMI);	
	for(i=0;i<len;i++)GEN[i]=(CLE[i]-'A')+CLE_NUM*(Ami_Clus[Ori_AA_Map[AMI[i]-'A']]-'0');
}
//---calculation_related-------// CLE->GLO
void Bioinfo_Code::AMI_CLE_transform(const char *CLE,int *GEN)
{
	int i;
	int len=(int)strlen(CLE);
	for(i=0;i<len;i++)GEN[i]=(CLE[i]-'A')+2*CLE_NUM;
}

//==================== Universal_Version ===================//
int Bioinfo_Code::Universal_Calc(int ii,int jj,int len,
	int *IN1,int *IN2,int mollen1,int mollen2,int InType)
{
#ifdef DEBUG
	overrange_debug_test(ii+len-1,mollen1);
	overrange_debug_test(jj+len-1,mollen2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
	int k;
	int score=0;
	switch(InType) 
	{ 
		case 1:   //-> CLE
		{
			for(k=0;k<len;k++)score+=Ori_CLESUM[IN1[ii+k]][IN2[jj+k]];
			return score;
		}
		case 2:   //-> GEN
		{
			for(k=0;k<len;k++)score+=Gen_CLESUM[IN1[ii+k]][IN2[jj+k]];
			return score;
		}
		case 3:   //-> AMI
		{
			for(k=0;k<len;k++)score+=Ori_BLOSUM[IN1[ii+k]][IN2[jj+k]];
			return score;
		}
		case 4:   //-> CLE+AMI
		{
			int a1,a2,c1,c2;
			for(k=0;k<len;k++)
			{
				a1=IN1[ii+k]/CLE_NUM;
				c1=IN1[ii+k]%CLE_NUM;
				a2=IN2[jj+k]/CLE_NUM;
				c2=IN2[jj+k]%CLE_NUM;
				score+=(int)(AMI_Weight*Ori_BLOSUM[a1][a2]+CLE_Weight*Ori_CLESUM[c1][c2]);
			}
			return score;
		}
		default:return 0;
	}
}
int Bioinfo_Code::Universal_Check(int *IN,int start,int len,int moln,int InType)
{
	int i;
	int rellen;
	int wsmax,curlen;
	rellen=len;

	//init_check
	if(start<=1)return 1;
	if(len<=2)return 1;
	if(start+len>=moln)rellen-=(start+len-moln);
	//check_tail
	if(start+rellen-2<0)return 1;
	//InType
	if(InType==3)return 1; //AMI return [ami==3]
	if((IN[start+rellen-1]%CLE_NUM==CLE_NUM-1)&&(IN[start+rellen-2]%CLE_NUM==CLE_NUM-1))return 0;
	//check_continuous_'R'
	wsmax=0;
	curlen=0;
	for(i=0;i<rellen;i++)
	{
		if(IN[i+start]%CLE_NUM==CLE_NUM-1)curlen++;
		else
		{
			if(curlen>wsmax)wsmax=curlen;
			curlen=0;
		}
	}
	if(curlen>wsmax)wsmax=curlen;
	if(wsmax>=3)return 0;
	//final
	return 1;
}
//--------- Universal Calculate SFP -----------//
void Bioinfo_Code::Universal_SFP(vector <SFP_Record> &SFP_List,int winlen,int ic,int step,
						         int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num)
{
	int i,j;
	int score;
	int temp[4];
	int wscount;
	int wsmem;
	int len;
	int head1;
	int head2;
	int head,tail;
	int count;

	//init
	for(i=0;i<4;i++)temp[i]=INT_MIN_NUM;
	score=INT_MIN_NUM;
	SFP_Record SFP;
//	SFP_List.clear();
	//process
	count=-1;
	for(i=0;i<mollen1+mollen2;i++)
	{
		if(i<mollen2)
		{
			head1=0;
			head2=mollen2-i;
			len=i;
		}
		else if(i>mollen1)
		{
			head1=i-mollen2;
			head2=0;
			len=mollen1+mollen2-i;
		}
		else
		{
			head1=i-mollen2;
			head2=0;
			len=mollen2;
		}
		if(len<winlen) continue;


		//__071122__//
		//temp_init
		{			
			temp[0]=INT_MIN_NUM;
			wscount=0;
			wsmem=0;
		}

		//calc_score
		score=Universal_Calc(head1,head2,winlen,IN1,IN2,mollen1,mollen2,InType);
		if(score>ic)
		{
			if(Universal_Check(IN1,head1,winlen,mollen1,InType)==1 && 
				Universal_Check(IN2,head2,winlen,mollen2,InType)==1)
			{
				temp[0]=score;
				temp[1]=head1;
				temp[2]=head2;
				temp[3]=winlen;
				wscount++;
				wsmem=0;
			}
		}

		//search_value
		for(j=1;j<len-winlen+1;j++)
		{
			head=Universal_Calc(head1-1+j,head2-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			tail=Universal_Calc(head1+winlen-1+j,head2+winlen-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			score=score-head+tail;
			//[1]whether there appears gap || check whther is full
			if((j!=(wsmem+1) && temp[0]!=INT_MIN_NUM)||wscount>step)
			{
				//push_back
				SFP.score=temp[0];
				if(isChange==0)
				{
					SFP.ii=temp[1];
					SFP.jj=temp[2];
				}
				else
				{
					SFP.ii=temp[2];
					SFP.jj=temp[1];
				}
				SFP.winlen=temp[3];
				count++;
				SFP_List[count]=SFP;
				//clear
				temp[0]=INT_MIN_NUM;
				wscount=0;
				wsmem=j;
			}
			//[2]check the neo_record
			if(score>ic)
			{
				if(Universal_Check(IN1,head1+j,winlen,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+j,winlen,mollen2,InType)==1)
				{
					if(score>temp[0])
					{
						temp[0]=score;
						temp[1]=head1+j;
						temp[2]=head2+j;
						temp[3]=winlen;
					}
					wscount++;
					wsmem=j;
				}
			}
		}//end of FOR(j)

		//final_check_2
		if(temp[0]!=INT_MIN_NUM)
		{
			//push_back
			SFP.score=temp[0];
			if(isChange==0)
			{
				SFP.ii=temp[1];
				SFP.jj=temp[2];
			}
			else
			{
				SFP.ii=temp[2];
				SFP.jj=temp[1];
			}
			SFP.winlen=temp[3];
			count++;
			SFP_List[count]=SFP;
		}
	}//end of FOR(i)
	num=count;
}
//create SFPs of two different length [shaving-include] (mollen1 must bigger than mollen2)
//step=0 degenerate to [shaving-exclude]
void Bioinfo_Code::Universal_SFP_II(vector <SFP_Record> &SFP_List1,vector <SFP_Record> &SFP_List2,int winlen1,int winlen2,int ic1,int ic2,
									int step1,int step2,int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num1,int &num2)
{
	int i,j,k;
	int score1,score2;
	int temp1[4],temp2[4];
	int wscount1,wscount2;
	int wsmem1,wsmem2;
	int len;
	int head1;
	int head2;
	int head,tail1,tail2;
	int count1,count2;

	//init
	for(i=0;i<4;i++)
	{
		temp1[i]=INT_MIN_NUM;
		temp2[i]=INT_MIN_NUM;
	}
	score1=INT_MIN_NUM;
	score2=INT_MIN_NUM;
	SFP_Record SFP;
//	SFP_List1.clear();
//	SFP_List2.clear();
	//start
	count1=-1;
	count2=-1;
	for(i=0;i<mollen1+mollen2;i++)
	{
		if(i<mollen2)
		{
			head1=0;
			head2=mollen2-i;
			len=i;
		}
		else if(i>mollen1)
		{
			head1=i-mollen2;
			head2=0;
			len=mollen1+mollen2-i;
		}
		else
		{
			head1=i-mollen2;
			head2=0;
			len=mollen2;
		}
		if(len<winlen1) continue;


		//__071122__//
		///temp_init
		{
			temp1[0]=INT_MIN_NUM;
			temp2[0]=INT_MIN_NUM;
			wscount1=0;  
			wscount2=0;  
			wsmem1=0;
			wsmem2=0;
		}

		//calc_score
		score1=Universal_Calc(head1,head2,winlen1,IN1,IN2,mollen1,mollen2,InType);
		if(score1>ic1)
		{
			if(Universal_Check(IN1,head1,winlen1,mollen1,InType)==1 && 
				Universal_Check(IN2,head2,winlen1,mollen2,InType)==1)
			{
				temp1[0]=score1;
				temp1[1]=head1;
				temp1[2]=head2;
				temp1[3]=winlen1;
				wscount1++;
				wsmem1=0;
			}
		}
		if(len>=winlen2)
		{
			score2=Universal_Calc(head1,head2,winlen2,IN1,IN2,mollen1,mollen2,InType);
			if(score2>ic2)
			{
				if(Universal_Check(IN1,head1,winlen2,mollen1,InType)==1 && 
					Universal_Check(IN2,head2,winlen2,mollen2,InType)==1)
				{
					temp2[0]=score2;
					temp2[1]=head1;
					temp2[2]=head2;
					temp2[3]=winlen2;
					wscount2++;
					wsmem2=0;
				}
			}
		}

		//search_value
		for(j=1;j<len-winlen2+1;j++)
		{
			head=Universal_Calc(head1-1+j,head2-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			tail1=Universal_Calc(head1+winlen1-1+j,head2+winlen1-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			tail2=Universal_Calc(head1+winlen2-1+j,head2+winlen2-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			score1=score1-head+tail1;
			score2=score2-head+tail2;
			//[1]whether there appears gap || check whther is full
			if((j!=(wsmem1+1) && temp1[0]!=INT_MIN_NUM)||wscount1>step1)
			{
				//push_back
				SFP.score=temp1[0];
				if(isChange==0)
				{
					SFP.ii=temp1[1];
					SFP.jj=temp1[2];
				}
				else
				{
					SFP.ii=temp1[2];
					SFP.jj=temp1[1];
				}
				SFP.winlen=temp1[3];
				count1++;
				SFP_List1[count1]=SFP;
				//clear
				temp1[0]=INT_MIN_NUM;
				wscount1=0;
				wsmem1=j;
			}
			if((j!=(wsmem2+1) && temp2[0]!=INT_MIN_NUM)||wscount2>step2)
			{
				//push_back
				SFP.score=temp2[0];
				if(isChange==0)
				{
					SFP.ii=temp2[1];
					SFP.jj=temp2[2];
				}
				else
				{
					SFP.ii=temp2[2];
					SFP.jj=temp2[1];
				}
				SFP.winlen=temp2[3];
				count2++;
				SFP_List2[count2]=SFP;
				//clear
				temp2[0]=INT_MIN_NUM;
				wscount2=0;
				wsmem2=j;
			}
			//[2]check the neo_record
			if(score1>ic1)
			{
				if(Universal_Check(IN1,head1+j,winlen1,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+j,winlen1,mollen2,InType)==1)
				{
					if(score1>temp1[0])
					{
						temp1[0]=score1;
						temp1[1]=head1+j;
						temp1[2]=head2+j;
						temp1[3]=winlen1;
					}
					wscount1++;
					wsmem1=j;
				}
			}
			if(score2>ic2)
			{
				if(Universal_Check(IN1,head1+j,winlen2,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+j,winlen2,mollen2,InType)==1)
				{
					if(score2>temp2[0])
					{
						temp2[0]=score2;
						temp2[1]=head1+j;
						temp2[2]=head2+j;
						temp2[3]=winlen2;
					}
					wscount2++;
					wsmem2=j;
				}
			}
		}//end of FOR(j)


		//final_check_2
		if(temp2[0]!=INT_MIN_NUM)
		{
			//push_back
			SFP.score=temp2[0];
			if(isChange==0)
			{
				SFP.ii=temp2[1];
				SFP.jj=temp2[2];
			}
			else
			{
				SFP.ii=temp2[2];
				SFP.jj=temp2[1];
			}
			SFP.winlen=temp2[3];
			count2++;
			SFP_List2[count2]=SFP;
		}
		for(k=j;k<len-winlen1+1;k++)
		{
			head=Universal_Calc(head1-1+k,head2-1+k,1,IN1,IN2,mollen1,mollen2,InType);
			tail1=Universal_Calc(head1+winlen1-1+k,head2+winlen1-1+k,1,IN1,IN2,mollen1,mollen2,InType);
			score1=score1-head+tail1;
			//[1]whether there appears gap || check whther is full
			if((k!=(wsmem1+1) && temp1[0]!=INT_MIN_NUM)||wscount1>step1)
			{
				//push_back
				SFP.score=temp1[0];
				if(isChange==0)
				{
					SFP.ii=temp1[1];
					SFP.jj=temp1[2];
				}
				else
				{
					SFP.ii=temp1[2];
					SFP.jj=temp1[1];
				}
				SFP.winlen=temp1[3];
				count1++;
				SFP_List1[count1]=SFP;
				//clear
				temp1[0]=INT_MIN_NUM;
				wscount1=0;
				wsmem1=k;
			}
			//[2]check the neo_record
			if(score1>ic1)
			{
				if(Universal_Check(IN1,head1+k,winlen1,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+k,winlen1,mollen2,InType)==1)
				{
					if(score1>temp1[0])
					{
						temp1[0]=score1;
						temp1[1]=head1+k;
						temp1[2]=head2+k;
						temp1[3]=winlen1;
					}
					wscount1++;
					wsmem1=k;
				}
			}
		}//end of FOR(k)

		//final_check_1
		if(temp1[0]!=INT_MIN_NUM)
		{
			//push_back
			SFP.score=temp1[0];
			if(isChange==0)
			{
				SFP.ii=temp1[1];
				SFP.jj=temp1[2];
			}
			else
			{
				SFP.ii=temp1[2];
				SFP.jj=temp1[1];
			}
			SFP.winlen=temp1[3];
			count1++;
			SFP_List1[count1]=SFP;
		}
	}//end of FOR(i)
	num1=count1+1;
	num2=count2+1;
}

//================================ Seed Generation ================================//
int Bioinfo_Code::Seed_Explosion(int temp[4],int ic,int limi,int &bk,int &fw,
								 int *IN1,int *IN2,int mollen1,int mollen2,int InType)
{
	int i;
	int score;
	int value;
	int pos1,pos2;
	int v1,v2;
	int s1,s2;
	int ori_sco;
	int ii,jj;
	int len;

	//init
	bk=0;
	fw=0;
	score=0;
	ori_sco=temp[0];
	ii=temp[1];
	jj=temp[2];
	len=temp[3];
	s1=0;
	s2=0;
	v1=1; //default:valid
	v2=1; //default:valid
	for(i=1;i<=limi;i++)
	{
		if(v1==0&&v2==0)break;
		//backward
		if(v1==0)goto forw;

		//[1]check position
		pos1=ii-i;
		pos2=jj-i;
		if(pos1<0||pos2<0)
		{
			v1=0;
			goto forw;
		}
		//[2]check broken
		if(InType!=2) //AMI return [ami==2]
		{
			if((IN1[pos1]%CLE_NUM==CLE_NUM-1)&&(IN2[pos2]%CLE_NUM==CLE_NUM-1))
			{
				v1=0;
				goto forw;
			}
		}
		//[3]calc score
		value=Universal_Calc(pos1,pos2,1,IN1,IN2,mollen1,mollen2,InType);
		if((ori_sco+value+s1+s2)<ic*(len+bk+fw))
		{
			v1=0;
			goto forw;
		}
		score+=value;
		s1+=value;
		bk++;

forw:
		//forward
		if(v2==0)continue;

		//[1]check position
		pos1=ii+len-1+i;
		pos2=jj+len-1+i;
		if(pos1>=mollen1||pos2>=mollen2)
		{
			v2=0;
			continue;
		}		
		//[2]check broken
		if(InType!=2) //AMI return [ami==2]
		{
			if((IN1[pos1]%CLE_NUM==CLE_NUM-1)&&(IN2[pos2]%CLE_NUM==CLE_NUM-1))
			{
				v2=0;
				continue;
			}
		}
		//[3]calc score
		value=Universal_Calc(pos1,pos2,1,IN1,IN2,mollen1,mollen2,InType);
		if((ori_sco+value+s1+s2)<ic*(len+bk+fw))
		{
			v2=0;
			continue;
		}
		score+=value;
		s2+=value;
		fw++;
	}

	return score;
}

//--------------------- GEN_version ----------------//
//create SFPs with SEED-EXPLOSION (mollen1 must bigger than mollen2)
void Bioinfo_Code::Universal_SFP_Seed(vector <SFP_Record> &SFP_List,int winlen,int ic,int extend,int step,
									   int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num)
{
	int i,j;
	int score;
	int temp[4];
	int count;
	int wscount;
	int wsmem;
	int len;
	int head1;
	int head2;
	int head,tail;
	int len2;
	int forw,back;
	int value;
	int cur_sco,cur_pos;
	int isFirst1,isFirst2;
	int isLast;


	//init
	SFP_Record SFP;
//	SFP_List.clear();
	count=-1;
	//start
	for(i=0;i<mollen1+mollen2;i++)
	{
		if(i<mollen2)
		{
			head1=0;
			head2=mollen2-i;
			len=i;
		}
		else if(i>mollen1)
		{
			head1=i-mollen2;
			head2=0;
			len=mollen1+mollen2-i;
		}
		else
		{
			head1=i-mollen2;
			head2=0;
			len=mollen2;
		}
		if(len<winlen) continue;


		//__071122__//
		//temp_init
		{			
			temp[0]=INT_MIN_NUM;
			wscount=0;
			wsmem=0;
			isFirst1=1;
			isFirst2=1;
			isLast=0;
		}

		//calc_score
		score=Universal_Calc(head1,head2,winlen,IN1,IN2,mollen1,mollen2,InType);
		if(score>ic*winlen)
		{
			if(Universal_Check(IN1,head1,winlen,mollen1,InType)==1 && 
				Universal_Check(IN2,head2,winlen,mollen2,InType)==1)
			{
				temp[0]=score;
				temp[1]=head1;
				temp[2]=head2;
				temp[3]=winlen;
				wscount++;
				wsmem=0;
			}
		}

		//search_value
		for(j=1;j<len-winlen+1;j++)
		{
			head=Universal_Calc(head1-1+j,head2-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			tail=Universal_Calc(head1+winlen-1+j,head2+winlen-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			score=score-head+tail;
			//[1]whether there appears gap || check whther is full
			if((j!=(wsmem+1) && temp[0]!=INT_MIN_NUM)||wscount>step)
			{
ws_last:
				//SEED-EXPLOSION[1]
				value=Seed_Explosion(temp,ic,extend,back,forw,IN1,IN2,mollen1,mollen2,InType);
				value+=temp[0];
				len2=temp[3]+back+forw;
				if(isFirst1==1)isFirst1=0;
				else
				{
					//init_substitute
					cur_sco=SFP_List[count].score;
					if(isChange==0)cur_pos=SFP_List[count].ii;
					else cur_pos=SFP_List[count].jj;
					//check_substitute
					if(abs(temp[1]-cur_pos)<step)
					{
						if(temp[0]>cur_sco)
						{
							//QQ1_substitute
							SFP_List[count].score=value;
							if(isChange==0)
							{
								SFP_List[count].ii=temp[1]-back;
								SFP_List[count].jj=temp[2]-back;
							}
							else
							{
								SFP_List[count].jj=temp[1]-back;
								SFP_List[count].ii=temp[2]-back;
							}
							SFP_List[count].winlen=len2;
							goto ws_next;
						}
						else goto ws_next;
					}
				}

				//QQ1_add
				//push_back
				SFP.score=value;
				if(isChange==0)
				{
					SFP.ii=temp[1]-back;
					SFP.jj=temp[2]-back;
				}
				else
				{
					SFP.ii=temp[2]-back;
					SFP.jj=temp[1]-back;
				}
				SFP.winlen=len2;
				count++;
				SFP_List[count]=SFP;


ws_next:
				//temp_init
				temp[0]=INT_MIN_NUM;
				wscount=0;
				wsmem=j;
				if(isLast==1)goto ws_last_bk;
			}

			//[2]check the neo_record
			if(score>ic*winlen)
			{
				if(Universal_Check(IN1,head1+j,winlen,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+j,winlen,mollen2,InType)==1)
				{
					if(score>temp[0])
					{
						temp[0]=score;
						temp[1]=head1+j;
						temp[2]=head2+j;
						temp[3]=winlen;
					}
					wscount++;
					wsmem=j;
				}
			}
		}//end of FOR(j)


		//final_check_2
		if(temp[0]!=INT_MIN_NUM)
		{
			isLast=1;
			goto ws_last;
ws_last_bk:
			;
		}
	}//end of FOR(i)
	num=count+1;
}
//create SFPs with SEED-EXPLOSION (mollen1 must bigger than mollen2)
void Bioinfo_Code::Universal_SFP_Seed_II(vector <SFP_Record> &SFP_List1,vector <SFP_Record> &SFP_List2,
										  int winlen1,int winlen2,int ic1,int ic2,int extend1,int extend2,int step1,int step2,
										  int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num1,int &num2)   
{
	int i,j;
	int score;
	int count1,count2;
	int temp[4];
	int wscount;
	int wsmem;
	int len;
	int head1;
	int head2;
	int head,tail;
	int len2;
	int forw,back;
	int value;
	int cur_sco,cur_pos;
	int isFirst1,isFirst2;
	int isLast;

	//init
	SFP_Record SFP;
//	SFP_List1.clear();
//	SFP_List2.clear();
	//start
	count1=-1;
	count2=-1;
	for(i=0;i<mollen1+mollen2;i++)
	{
		if(i<mollen2)
		{
			head1=0;
			head2=mollen2-i;
			len=i;
		}
		else if(i>mollen1)
		{
			head1=i-mollen2;
			head2=0;
			len=mollen1+mollen2-i;
		}
		else
		{
			head1=i-mollen2;
			head2=0;
			len=mollen2;
		}
		if(len<winlen1) continue;


		//__071122__//
		//temp_init
		{			
			temp[0]=INT_MIN_NUM;
			wscount=0;
			wsmem=0;
			isFirst1=1;
			isFirst2=1;
			isLast=0;
		}

		//calc_score
		score=Universal_Calc(head1,head2,winlen1,IN1,IN2,mollen1,mollen2,InType);
		if(score>ic1*winlen1)
		{
			if(Universal_Check(IN1,head1,winlen1,mollen1,InType)==1 && 
				Universal_Check(IN2,head2,winlen1,mollen2,InType)==1)
			{
				temp[0]=score;
				temp[1]=head1;
				temp[2]=head2;
				temp[3]=winlen1;
				wscount++;
				wsmem=0;
			}
		}

		//search_value
		for(j=1;j<len-winlen1+1;j++)
		{
			head=Universal_Calc(head1-1+j,head2-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			tail=Universal_Calc(head1+winlen1-1+j,head2+winlen1-1+j,1,IN1,IN2,mollen1,mollen2,InType);
			score=score-head+tail;
			//[1]whether there appears gap || check whther is full
			if((j!=(wsmem+1) && temp[0]!=INT_MIN_NUM)||wscount>step1)
			{
ws_last:
				//SEED-EXPLOSION[1]
				value=Seed_Explosion(temp,ic1,extend1,back,forw,IN1,IN2,mollen1,mollen2,InType);
				value+=temp[0];
				len2=temp[3]+back+forw;
				if(isFirst1==1)isFirst1=0;
				else
				{
					//init_substitute
					cur_sco=SFP_List1[count1].score;
					if(isChange==0)cur_pos=SFP_List1[count1].ii;
					else cur_pos=SFP_List1[count1].jj;

					//check_substitute
					if(abs(temp[1]-cur_pos)<step1)
					{
						if(temp[0]>cur_sco)
						{
							//QQ1_substitute
							SFP_List1[count1].score=value;
							if(isChange==0)
							{
								SFP_List1[count1].ii=temp[1]-back;
								SFP_List1[count1].jj=temp[2]-back;
							}
							else
							{
								SFP_List1[count1].jj=temp[1]-back;
								SFP_List1[count1].ii=temp[2]-back;
							}
							SFP_List1[count1].winlen=len2;
							goto ws_next;
						}
						else goto ws_next;
					}
				}

				//QQ1_add
				//push_back
				SFP.score=value;
				if(isChange==0)
				{
					SFP.ii=temp[1]-back;
					SFP.jj=temp[2]-back;
				}
				else
				{
					SFP.ii=temp[2]-back;
					SFP.jj=temp[1]-back;
				}
				SFP.winlen=len2;
				count1++;
				SFP_List1[count1]=SFP;

ws_next:
				//SEED-EXPLOSION[2]
				if(temp[0]<ic2*temp[3])goto ws_mext;
				value=Seed_Explosion(temp,ic2,extend2,back,forw,IN1,IN2,mollen1,mollen2,InType);
				value+=temp[0];
				len2=temp[3]+back+forw;
				if(len2<winlen2)goto ws_mext;
				if(isFirst2==1)isFirst2=0;
				else
				{
					//init_substitute
					cur_sco=SFP_List2[count2].score;
					if(isChange==0)cur_pos=SFP_List2[count2].ii;
					else cur_pos=SFP_List2[count2].jj;
					//check_substitute
					if(abs(temp[1]-back-cur_pos)<step2)
					{
						if(value>cur_sco)
						{
							//QQ2_substitute
							SFP_List2[count2].score=value;
							if(isChange==0)
							{
								SFP_List2[count2].ii=temp[1]-back;
								SFP_List2[count2].jj=temp[2]-back;
							}
							else
							{
								SFP_List2[count2].jj=temp[1]-back;
								SFP_List2[count2].ii=temp[2]-back;
							}
							SFP_List2[count2].winlen=len2;
							goto ws_mext;
						}
						else goto ws_mext;
					}
				}

				//QQ2_add
				//push_back
				SFP.score=value;
				if(isChange==0)
				{
					SFP.ii=temp[1]-back;
					SFP.jj=temp[2]-back;
				}
				else
				{
					SFP.ii=temp[2]-back;
					SFP.jj=temp[1]-back;
				}
				SFP.winlen=len2;
				count2++;
				SFP_List2[count2]=SFP;


ws_mext:
				//temp_init
				temp[0]=INT_MIN_NUM;
				wscount=0;
				wsmem=j;
				if(isLast==1)goto ws_last_bk;
			}

			//[2]check the neo_record
			if(score>ic1*winlen1)
			{
				if(Universal_Check(IN1,head1+j,winlen1,mollen1,InType)==1 && 
					Universal_Check(IN2,head2+j,winlen1,mollen2,InType)==1)
				{
					if(score>temp[0])
					{
						temp[0]=score;
						temp[1]=head1+j;
						temp[2]=head2+j;
						temp[3]=winlen1;
					}
					wscount++;
					wsmem=j;
				}
			}
		}//end of FOR(j)


		//final_check_2
		if(temp[0]!=INT_MIN_NUM)
		{
			isLast=1;
			goto ws_last;
ws_last_bk:
			;
		}
	}//end of FOR(i)
	//final
	num1=count1+1;
	num2=count2+1;
}

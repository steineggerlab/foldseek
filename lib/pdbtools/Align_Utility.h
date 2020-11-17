#pragma once
#include <vector>
using namespace std;

//--------- type_def -----------//_Align_Record
class Align_Record
{
public:
	Align_Record(void){main_sco=-1.0,RMSD=-1.0,lali=0;};
	~Align_Record(void){};
	void operator=(const Align_Record & align_record)
	{
		alignment=align_record.alignment;
		score=align_record.score;
		rotmat=align_record.rotmat;
		index=align_record.index;
		main_sco=align_record.main_sco;
		RMSD=align_record.RMSD;
		TMsco=align_record.TMsco;
		lali=align_record.lali;
	}
	bool operator < (const Align_Record &m)const
	{
		return main_sco < m.main_sco;
	}
public:
	vector <int> alignment; //default: record ali2 (mol2 is fixed!!)
	vector <double> score;  //default: record all scores
	vector <double> rotmat; //default: mol2 is fixed!!
	vector <int> index;     //default: NULL
	double main_sco;        //default: TMscore
	double RMSD;
	double TMsco;
	int lali;
};

//--------- type_def -----------//_SFP_Record
class SFP_Record
{
public:
	SFP_Record(void){score=-1,ii=-1,jj=-1,winlen=0;};
	~SFP_Record(void){};
	void operator=(const SFP_Record & sfp_record)
	{
		score=sfp_record.score;
		ii=sfp_record.ii;
		jj=sfp_record.jj;
		winlen=sfp_record.winlen;
	}
	bool operator < (const SFP_Record &m)const
	{
		return score < m.score;
	}
public:
	int score;
	int ii;
	int jj;
	int winlen;
};

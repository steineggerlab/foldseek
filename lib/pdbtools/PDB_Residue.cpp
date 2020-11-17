#include "Backbone_Sidechain.h"
#include "Utility.h"
#include "PDB_Residue.h"
#include "PDB_Utility.h"
#include <iostream>
using namespace std;

//-------- PDB_Residue_Define ---------// 
PDB_Residue::PDB_Residue(void)
{
	AA_residue = 'X';
	PDB_residue_number.assign("A     "); //0: chain_id; 1-4: residue number; 5: insert code
	hydro_rec.clear();
}
PDB_Residue::~PDB_Residue(void)
{
}
PDB_Residue::PDB_Residue(const PDB_Residue & pdb)
{
	this->AA_residue=pdb.AA_residue;
	this->PDB_residue_number.assign(pdb.PDB_residue_number);
	this->backbone=pdb.backbone;
	this->sidechain=pdb.sidechain;
	this->hydro_rec=pdb.hydro_rec;
}
PDB_Residue & PDB_Residue::operator =(const PDB_Residue & pdb)
{
	if(this==&pdb)return *this;
	this->AA_residue=pdb.AA_residue;
	this->PDB_residue_number.assign(pdb.PDB_residue_number);
	this->backbone=pdb.backbone;
	this->sidechain=pdb.sidechain;
	this->hydro_rec=pdb.hydro_rec;
	return *this;
}

//--------- init --------// 
int PDB_Residue::PDB_residue_backbone_initialize(char c)
{
	int num;
	num=AA26_to_AA20(c-'A');
	if(num==-1)return -1;
	if(num==20)AA_residue='X';
	else AA_residue=c;
	backbone.backbone_sidechain_initialize(4);
	return 0;
}
int PDB_Residue::PDB_residue_sidechain_initialize(char c)
{
	int num;
	num=AA26_to_AA20(c-'A');
	if(num==-1)return -1;
	if(num==20)AA_residue='X';
	else AA_residue=c;
	num=AA26_sidechain_size(c-'A');
	sidechain.backbone_sidechain_initialize(num);
	return 0;
}
//-------- check ----------// 
int PDB_Residue::PDB_residue_backbone_check(int check)
{
	int i;
	int num;
	num=this->backbone.get_atom_totnum();
	if(num<check)return 0;
	for(i=0;i<check;i++)if(this->backbone.get_part_index(i)==0)return -1*(i+1);
	return 1;
}
int PDB_Residue::PDB_residue_CB_check(void)
{
	int num;
	num=this->sidechain.get_atom_totnum();
	if(num<1)return 0;
	if(this->sidechain.get_part_index(0)==0)return -1;
	return 1;
}
int PDB_Residue::PDB_residue_sidechain_check(void)
{
	int i;
	int num;
	char name;
	name=this->AA_residue;
	num=this->sidechain.get_atom_totnum();
	if(name<'A'||name>'Z')return 0;
	if(num!=AA26_sidechain_size(name-'A'))return 0;
	for(i=0;i<num;i++)if(this->sidechain.get_part_index(i)==0)return -1*(i+1);
	return 1;
}

//-------- Get and Set -------//
//[backbone]
//(atom_name mode)
int PDB_Residue::set_backbone_atom( const string &atom_name, const XYZ &atom, int numb, double R, double tmpr )
{
	return backbone.set_backbone_atom(atom_name,atom,numb,R,tmpr);
}
int PDB_Residue::get_backbone_atom( const string &atom_name, XYZ &atom, int & numb, double & R, double & tmpr ) const
{
	return backbone.get_backbone_atom(atom_name,atom,numb,R,tmpr);
}
int PDB_Residue::get_backbone_atom( const string &atom_name, XYZ &atom ) const
{
	return backbone.get_backbone_atom(atom_name,atom);
}
int PDB_Residue::get_backbone_atom( const string &atom_name, int &numb ) const
{
	return backbone.get_backbone_atom(atom_name,numb);
}
//(index mode)
int PDB_Residue::set_backbone_atom( int index, const XYZ &atom, int numb, double R, double tmpr )
{
	return backbone.set_atom(index,atom,numb, R, tmpr);
}
int PDB_Residue::get_backbone_atom( int index, XYZ &atom, int & numb, double & R, double & tmpr ) const
{
	return backbone.get_atom(index,atom,numb, R, tmpr);
}
int PDB_Residue::get_backbone_atom( int index, XYZ &atom ) const
{
	return backbone.get_atom(index,atom);
}
int PDB_Residue::get_backbone_atom( int index, int &numb ) const
{
	return backbone.get_atom(index,numb);
}

//[sidechain]
//(atom_name mode)
int PDB_Residue::set_sidechain_atom( const string &atom_name, const XYZ &atom, int numb, double R, double tmpr )
{
	return sidechain.set_sidechain_atom(atom_name,AA_residue,atom,numb,R,tmpr);
}
int PDB_Residue::get_sidechain_atom( const string &atom_name, XYZ &atom, int & numb, double & R, double & tmpr ) const
{
	return sidechain.get_sidechain_atom(atom_name,AA_residue,atom,numb,R,tmpr);
}
int PDB_Residue::get_sidechain_atom( const string &atom_name, XYZ &atom ) const
{
	return sidechain.get_sidechain_atom(atom_name,AA_residue,atom);
}
int PDB_Residue::get_sidechain_atom( const string &atom_name, int &numb ) const
{
	return sidechain.get_sidechain_atom(atom_name,AA_residue,numb);
}
//(index mode)
int PDB_Residue::set_sidechain_atom( int index, const XYZ &atom, int numb, double R, double tmpr )
{
	return sidechain.set_atom(index,atom,numb, R, tmpr);
}
int PDB_Residue::get_sidechain_atom( int index, XYZ &atom, int & numb, double & R, double & tmpr ) const
{
	return sidechain.get_atom(index,atom,numb, R, tmpr);
}
int PDB_Residue::get_sidechain_atom( int index, XYZ &atom ) const
{
	return sidechain.get_atom(index,atom);
}
int PDB_Residue::get_sidechain_atom( int index, int &numb ) const
{
	return sidechain.get_atom(index,numb);
}

//[set and get] (others)
int PDB_Residue::get_backbone_totnum(void) const
{
	return backbone.get_atom_totnum();
}
int PDB_Residue::get_sidechain_totnum(void) const
{
	return sidechain.get_atom_totnum();
}
int PDB_Residue::get_backbone_part_index(int i) const
{
	return backbone.get_part_index(i);
}
int PDB_Residue::get_sidechain_part_index(int i) const
{
	return sidechain.get_part_index(i);
}
int PDB_Residue::set_AA(char c)
{
	AA_residue = c;
	return 0;
}
char PDB_Residue::get_AA() const
{
	return AA_residue;
}
int PDB_Residue::set_PDB_residue_number(const string & res_num)
{
	int len=(int)res_num.length();
	if(len!=6)
	{
		cout << "WARINING!! -> PDB_Residue::set_PDB_residue_number -> invalid length: "<< len << endl;
		cout << "res_num: " << res_num << endl;
		return -1;
	}
	PDB_residue_number.assign(res_num);
	return 0;
}
int PDB_Residue::get_PDB_residue_number(string &res_num) const
{
	res_num.assign(PDB_residue_number);
	return 0;
}

//=========== PDB_Residue to XYZ_array + sidechain_count ==========//__130830__//
//XYZ_array MUST have at least N,CA,C,O,CB five atoms! (for glycine, CB=CA)
//sidechain_count records the additional sidechain atoms!
int PDB_Residue::get_XYZ_array(XYZ *XYZ_array,int &sidechain_count)
{
	sidechain_count=0;
	int i;
	XYZ xyz;
	int retv;
	//assign CA
	retv=PDB_residue_backbone_check();
	if(retv!=1)return retv;
	for(i=0;i<4;i++)
	{
		backbone.get_atom(i,xyz);
		XYZ_array[i]=xyz;
	}
	//assign CB
	retv=PDB_residue_CB_check();
	if(retv!=1)
	{
//		backbone.get_atom(1,xyz);
//		XYZ_array[4]=xyz;
	}
	else
	{
		if(AA_residue!='G')
		{
			sidechain.get_atom(0,xyz);
			XYZ_array[4]=xyz;
			sidechain_count++;
		}
	}
	//assign sidechain
	int number=get_sidechain_totnum();
	for(i=1;i<number;i++)
	{
		if(AA_residue!='G')
		{
			retv=get_sidechain_part_index(i);
			if(retv==0)continue;
			sidechain.get_atom(i,xyz);
			XYZ_array[4+i]=xyz;
			sidechain_count++;
		}
	}
	//return
	return 1;
}

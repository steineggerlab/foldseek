#include "Backbone_Sidechain.h"
#include "Utility.h"
#include "PDB_Utility.h"
#include <string>
using namespace std;

//============================// 
//-------- Backbone_Sidechain --------// 
//default: num=10, because the maximal sidechain is 10 for W
Backbone_Sidechain::Backbone_Sidechain(void)
{
	atom_totnum=0;
}  
Backbone_Sidechain::~Backbone_Sidechain(void)
{
}
Backbone_Sidechain::Backbone_Sidechain(const Backbone_Sidechain & amino)
{
	this->atom_totnum = amino.atom_totnum;
	for (int i = 0; i < this->atom_totnum; i++)
	{
		  this->part_index[i] = amino.part_index[i];
		  this->atom_index[i] = amino.atom_index[i];
		  this->atom_serial_number[i] = amino.atom_serial_number[i]; //dbu 2010/12/12
		  this->R_factor[i] = amino.R_factor[i];
		  this->temperature[i] = amino.temperature[i];
	}
}
Backbone_Sidechain & Backbone_Sidechain::operator =(const Backbone_Sidechain & amino)
{
	if(this==&amino)return *this;
	this->atom_totnum = amino.atom_totnum;
	for (int i = 0; i < this->atom_totnum; i++)
	{
		  this->part_index[i] = amino.part_index[i];
		  this->atom_index[i] = amino.atom_index[i];
		  this->atom_serial_number[i] = amino.atom_serial_number[i]; //dbu 2010/12/12
		  this->R_factor[i] = amino.R_factor[i];
		  this->temperature[i] = amino.temperature[i];
	}
	return *this;
}

//--------- basic function -------// 
int Backbone_Sidechain::backbone_sidechain_initialize(int number,XYZ a,int numb,double R,double tmpr)
{
	if(number<=0||number>10)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::initialize. number = %d\n", number);
		return -1;
	}
	atom_totnum=number;
	for(int i=0;i<number;i++)
	{
		part_index[i]=0;
		atom_index[i]=a;
		atom_serial_number[i]=numb;
		R_factor[i]=R;
		temperature[i]=tmpr;
	}
	return 0;
}
Backbone_Sidechain *Backbone_Sidechain::born(void)
{
	Backbone_Sidechain *ami = new Backbone_Sidechain;
	ami->atom_totnum=this->atom_totnum;
	for (int i = 0; i < this->atom_totnum; i++)
	{
		  ami->part_index[i] = this->part_index[i];
		  ami->atom_index[i] = this->atom_index[i];
		  ami->atom_serial_number[i] = this->atom_serial_number[i]; //dbu 2010/12/12
		  ami->R_factor[i] = this->R_factor[i];
		  ami->temperature[i] = this->temperature[i];
	} 
	return ami;
}

//------------ other functions ------------//
//--- normal_atom ---//
int Backbone_Sidechain::set_atom(int index, const XYZ &atom, int numb, double R, double tmpr)
{
	if(index<0||index>=atom_totnum)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::set_atom. index = %d, totnum = %d\n", index, atom_totnum);
		return -1;
	}
	part_index[index] = 1;
	atom_index[index] = atom;
	atom_serial_number[index] = numb;
	R_factor[index] = R;
	temperature[index] = tmpr;
	return 0;
}
int Backbone_Sidechain::get_atom(int index, XYZ & atom, int & numb, double & R, double & tmpr) const
{
	if(index<0||index>=atom_totnum)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom. index = %d, totnum = %d\n", index, atom_totnum);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom. part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	numb = atom_serial_number[index];
	R = R_factor[index];
	tmpr = temperature[index];
	return 0;
}
int Backbone_Sidechain::get_atom(int index, XYZ &atom ) const  //-> return atom_XYZ
{
	if(index<0||index>=atom_totnum)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom (XYZ). index = %d, totnum = %d\n", index, atom_totnum);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom (XYZ). part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	return 0;
}
int Backbone_Sidechain::get_atom(int index, int &numb ) const  //-> return atom_serial_num
{
	if(index<0||index>=atom_totnum)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom (numb). index = %d, totnum = %d\n", index, atom_totnum);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_atom (numb). part_index[%d]==0\n", index);
		return -1;
	}
	numb = atom_serial_number[index];
	return 0;
}

//--- backbone_atom ---//
int Backbone_Sidechain::set_backbone_atom(const string &atom_name, const XYZ &atom, int numb, double R, double tmpr)
{
	int index = backbone_atom_name_encode(atom_name.c_str());
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::set_backbone_atom. atom_name = %s, index = %d\n", 
			atom_name.c_str(), index);
		return -1;
	}
	part_index[index] = 1;
	atom_index[index] = atom;
	atom_serial_number[index] = numb;
	R_factor[index] = R;
	temperature[index] = tmpr;
	return 0;
}
int Backbone_Sidechain::get_backbone_atom(const string &atom_name, XYZ & atom, int & numb, double & R, double & tmpr) const
{
	int index = backbone_atom_name_encode(atom_name.c_str());
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom. atom_name = %s, index = %d\n", 
			atom_name.c_str(), index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom. part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	numb = atom_serial_number[index];
	R = R_factor[index];
	tmpr = R_factor[index];
	return 0;
}
int Backbone_Sidechain::get_backbone_atom(const string &atom_name, XYZ & atom ) const  //-> return atom_XYZ
{
	int index = backbone_atom_name_encode(atom_name.c_str());
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom (XYZ). atom_name = %s, index = %d\n", 
			atom_name.c_str(), index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom (XYZ). part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	return 0;
}
int Backbone_Sidechain::get_backbone_atom(const string &atom_name, int & numb ) const  //-> return atom_serial_num
{
	int index = backbone_atom_name_encode(atom_name.c_str());
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom (numb). atom_name = %s, index = %d\n", 
			atom_name.c_str(), index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_backbone_atom (numb). part_index[%d]==0\n", index);
		return -1;
	}
	numb = atom_serial_number[index];
	return 0;
}

//--- sidechain_atom ---//
int Backbone_Sidechain::set_sidechain_atom(const string &atom_name, char amino, const XYZ &atom, int numb, double R, double tmpr)
{
	int index = sidechain_atom_name_encode(atom_name.c_str(), amino);
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::set_sidechain_atom. atom_name = %s, amino = %c, index = %d\n", 
			atom_name.c_str(), amino, index);
		return -1;
	}
	part_index[index] = 1;
	atom_index[index] = atom;
	atom_serial_number[index] = numb;
	R_factor[index] = R;
	temperature[index] = tmpr;
	return 0;
}
int Backbone_Sidechain::get_sidechain_atom(const string &atom_name, char amino, XYZ & atom, int & numb, double & R, double & tmpr) const
{
	int index = sidechain_atom_name_encode(atom_name.c_str(), amino);
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom. atom_name = %s, amino = %c, index = %d\n", 
			atom_name.c_str(), amino, index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom. part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	numb = atom_serial_number[index];
	R = R_factor[index];
	tmpr = temperature[index];
	return 0;
}
int Backbone_Sidechain::get_sidechain_atom(const string &atom_name, char amino, XYZ & atom ) const  //-> return atom_XYZ
{
	int index = sidechain_atom_name_encode(atom_name.c_str(), amino);
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom (XYZ). atom_name = %s, amino = %c, index = %d\n", 
			atom_name.c_str(), amino, index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom (XYZ). part_index[%d]==0\n", index);
		return -1;
	}
	atom = atom_index[index];
	return 0;
}
int Backbone_Sidechain::get_sidechain_atom(const string &atom_name, char amino, int & numb ) const  //-> return atom_serial_num
{
	int index = sidechain_atom_name_encode(atom_name.c_str(), amino);
	if(index < 0)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom (numb). atom_name = %s, amino = %c, index = %d\n", 
			atom_name.c_str(), amino, index);
		return -1;
	}
	if(part_index[index]==0) 
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_sidechain_atom (numb). part_index[%d]==0\n", index);
		return -1;
	}
	numb = atom_serial_number[index];
	return 0;
}

//--- others ---//
int Backbone_Sidechain::get_part_index(int i) const
{
	if(i<0||i>=atom_totnum)
	{
		debug_info(DEBUG_WARNING, "@Backbone_Sidechain::get_part_index. index = %d, totnum = %d\n", i, atom_totnum);
		return -1;
	}
	return part_index[i];
}
int Backbone_Sidechain::get_atom_totnum() const
{
	return atom_totnum;
}

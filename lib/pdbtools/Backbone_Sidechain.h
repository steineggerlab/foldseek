#pragma once
#include "XYZ.h"
#include <string>
#include <cstring>
using namespace std;

//====class_Backbone_Sidechain====// 
class Backbone_Sidechain
{
public:
	Backbone_Sidechain(void); //maximal sidechain is 10 for W
	~Backbone_Sidechain(void);
	Backbone_Sidechain(const Backbone_Sidechain & amino);
	Backbone_Sidechain & operator =(const Backbone_Sidechain & amino);

//======== variables ========//
private:
	int atom_totnum;	        //total atom 
	int part_index[10];	        //max is TRP, atom index 
	XYZ atom_index[10];	        //max is TRP, atom record
	int atom_serial_number[10]; //atom serial number, dbu 2010/12/12
	double R_factor[10];	    //R_factor   (default: 1.0)
	double temperature[10];     //temprature (default: 20.0)

//======== functions ========//
public:
	//normal_function
	int backbone_sidechain_initialize(int number,XYZ a=0.0,int numb=0,double R=1.0,double tmpr=20.0 );
	Backbone_Sidechain *born(void);
	//atom
	int set_atom( int index, const XYZ &atom, int numb = 0, double R = 1.0, double tmpr = 20.0 );
	int get_atom( int index, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_atom( int index, XYZ &atom ) const;
	int get_atom( int index, int &numb ) const;
	// backbone_atom
	int set_backbone_atom( const string &atom_name, const XYZ &atom, int numb = 0, double R = 1.0, double tmpr = 20.0 );
	int get_backbone_atom( const string &atom_name, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_backbone_atom( const string &atom_name, XYZ &atom ) const;
	int get_backbone_atom( const string &atom_name, int &numb ) const;
	// sidechain_atom
	int set_sidechain_atom( const string &atom_name, char amino, const XYZ &atom, int numb = 0, double R = 1.0, double tmpr = 20.0 );
	int get_sidechain_atom( const string &atom_name, char amino, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_sidechain_atom( const string &atom_name, char amino, XYZ &atom ) const;
	int get_sidechain_atom( const string &atom_name, char amino, int &numb ) const;
	//others
	int get_part_index(int i) const;
	int get_atom_totnum() const;
};
//====class_Backbone_Sidechain====//over 

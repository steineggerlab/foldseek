#pragma once
#include <cstdio>
#include <math.h>

//====class_XYZ====// 
//=> Euclidean distance coordinate definition 
class XYZ
{
	//-------- pulic_data ---------// 
public:
	double X, Y, Z;
	//---- constructor ---//
public:
	XYZ(void);
	~XYZ(void);
	XYZ(double w);
	XYZ(double x, double y, double z);
	XYZ(const XYZ & xyz);
	XYZ & operator =(const XYZ & xyz);
	//------ operator_reload ------//
public:
	XYZ operator +(const XYZ & xyz) const;
	XYZ operator -(const XYZ & xyz) const;
	XYZ operator *(double i) const;
	XYZ operator /(double i) const;
	void operator +=(const XYZ & xyz);
	void operator -=(const XYZ & xyz);
	void operator *=(double i);
	void operator /=(double i);
	void operator =(double i);
	//--------- function ----------//
public:
	void XYZ2();
	double norm() const;
	void transform(XYZ & xyz, const double *rotmat);
	double distance(const XYZ & xyz) const;	//
	double distance_square(const XYZ & xyz) const;	//
	double coordinate_seperation(XYZ & xyz) const;	//
	void xyz2double(double t[3]) const;
	void double2xyz(double t[3]);
	void cross_product(const XYZ & xyz, XYZ & ans) const;
	//angle related
	double calc_bend_angle(const XYZ & xyz) const;
	double calc_tort_angle(const XYZ & b, const XYZ & c) const;
	double point_product(const XYZ & xyz) const;
	double mixed_product(const XYZ & b, const XYZ & c) const;
};
//====class_XYZ====//over 

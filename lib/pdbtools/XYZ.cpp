#include "XYZ.h"

//-------------- constructor --------------//
XYZ::XYZ(void)
{
	X = 0;
	Y = 0;
	Z = 0;
}
XYZ::~XYZ(void)
{
}
XYZ::XYZ(double w)
{
	X = w;
	Y = w;
	Z = w;
}
XYZ::XYZ(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}
XYZ::XYZ(const XYZ & xyz)
{
	this->X = xyz.X;
	this->Y = xyz.Y;
	this->Z = xyz.Z;
}
XYZ & XYZ::operator =(const XYZ & xyz)
{
	if(this==&xyz)return *this;
	this->X = xyz.X;
	this->Y = xyz.Y;
	this->Z = xyz.Z;
	return *this;
}

//---------------operator_reload ----------------------// 
XYZ XYZ::operator +(const XYZ & xyz) const
{
	XYZ temp;
	temp.X = X + xyz.X;
	temp.Y = Y + xyz.Y;
	temp.Z = Z + xyz.Z;
	return (temp);
}
XYZ XYZ::operator -(const XYZ & xyz) const
{
	XYZ temp;
	temp.X = X - xyz.X;
	temp.Y = Y - xyz.Y;
	temp.Z = Z - xyz.Z;
	return (temp);
}
XYZ XYZ::operator *(double i) const
{
	XYZ temp;
	temp.X = X * i;
	temp.Y = Y * i;
	temp.Z = Z * i;
	return (temp);
}
XYZ XYZ::operator /(double i) const
{
	XYZ temp;
	temp.X = X / i;
	temp.Y = Y / i;
	temp.Z = Z / i;
	return (temp);
}
void XYZ::operator +=(const XYZ & xyz)
{
	this->X += xyz.X;
	this->Y += xyz.Y;
	this->Z += xyz.Z;
}
void XYZ::operator -=(const XYZ & xyz)
{
	this->X -= xyz.X;
	this->Y -= xyz.Y;
	this->Z -= xyz.Z;
}
void XYZ::operator *=(double i)
{
	this->X *= i;
	this->Y *= i;
	this->Z *= i;
}
void XYZ::operator /=(double i)
{
	this->X /= i;
	this->Y /= i;
	this->Z /= i;
}
void XYZ::operator =(double i)
{
	this->X = i;
	this->Y = i;
	this->Z = i;
}

//---------------function ----------------------// 
void XYZ::XYZ2()
{
	X = X * X;
	Y = Y * Y;
	Z = Z * Z;
}
double XYZ::norm() const
{
	double ans = 0;
	ans += X * X;
	ans += Y * Y;
	ans += Z * Z;
	return sqrt(ans);
}
double XYZ::distance(const XYZ & xyz) const
{
	XYZ xyz_ = xyz;
	xyz_ -= *this;
	xyz_.XYZ2();
	double sum = xyz_.X + xyz_.Y + xyz_.Z;
	return (sum > 0.0 ? sqrt(sum) : 0.0);
}
double XYZ::distance_square(const XYZ & xyz)  const 
{
	XYZ xyz_ = xyz;
	xyz_ -= *this;
	xyz_.XYZ2();
	double sum = xyz_.X + xyz_.Y + xyz_.Z;
	return (sum);
}
double XYZ::coordinate_seperation(XYZ & xyz)  const 
{
	XYZ xyz_ = xyz;
	xyz_ -= *this;
	double max =
	(fabs(xyz_.X) > fabs(xyz_.Y) ? fabs(xyz_.X) : fabs(xyz_.Y));
	return (max > fabs(xyz_.Z) ? max : fabs(xyz_.Z));
}
void XYZ::xyz2double(double t[3]) const 
{
	t[0] = this->X;
	t[1] = this->Y;
	t[2] = this->Z;
} 
void XYZ::double2xyz(double t[3])
{
	this->X = t[0];
	this->Y = t[1];
	this->Z = t[2];
}
void XYZ::transform(XYZ & xyz, const double *rotmat)
{
	double dx,dy,dz;
	dx=this->X;
	dy=this->Y;
	dz=this->Z;
	xyz.X=dx*rotmat[0]+dy*rotmat[1]+dz*rotmat[2]+rotmat[9];
	xyz.Y=dx*rotmat[3]+dy*rotmat[4]+dz*rotmat[5]+rotmat[10];
	xyz.Z=dx*rotmat[6]+dy*rotmat[7]+dz*rotmat[8]+rotmat[11];
}
void XYZ::cross_product(const XYZ & xyz, XYZ & ans) const
{
	ans.X = this->Y * xyz.Z - this->Z * xyz.Y;
	ans.Y = -1 * this->X * xyz.Z + this->Z * xyz.X;
	ans.Z = this->X * xyz.Y - this->Y * xyz.X;
}
double XYZ::point_product(const XYZ & xyz) const
{
	return this->X * xyz.X + this->Y * xyz.Y + this->Z * xyz.Z;
}

//---- angle related ----//
double XYZ::calc_bend_angle(const XYZ & xyz) const 
{
	if(this->norm() < 1e-16) return 0;
	if(xyz.norm() < 1e-16) return 0;
	double cosine = point_product(xyz) / this->norm() / xyz.norm();
	return acos(cosine);
}
double XYZ::mixed_product(const XYZ & b, const XYZ & c) const
{
	XYZ t;
	cross_product(b, t);
	return t.point_product(c);
}
double XYZ::calc_tort_angle(const XYZ & b, const XYZ & c) const
{
	XYZ t1, t2;
	cross_product(b, t1);
	b.cross_product(c, t2);
	double angle = t1.calc_bend_angle(t2);
	double mp = mixed_product(b, c);
	if(mp >= 0) return angle;
	else return -1 * angle;
}


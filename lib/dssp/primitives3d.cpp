//
// Created by Martin Steinegger on 28.09.18.
//

#include "primitives3d.h"
// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// 3d routines

#include "mas.h"
#include <cstdint>
#include <valarray>
#include <cmath>


using namespace std;

const double
        kPI = 4 * std::atan(1.0);

// --------------------------------------------------------------------


// --------------------------------------------------------------------

double MPoint::Normalize()
{
    double length = mX * mX + mY * mY + mZ * mZ;
    if (length > 0)
    {
        length = sqrt(length);
        mX /= length;
        mY /= length;
        mZ /= length;
    }
    return length;
}

MPoint operator+(const MPoint& lhs, const MPoint& rhs)
{
    return MPoint(lhs.mX + rhs.mX, lhs.mY + rhs.mY, lhs.mZ + rhs.mZ);
}

MPoint operator-(const MPoint& lhs, const MPoint& rhs)
{
    return MPoint(lhs.mX - rhs.mX, lhs.mY - rhs.mY, lhs.mZ - rhs.mZ);
}

MPoint operator-(const MPoint& pt)
{
    return MPoint(-pt.mX, -pt.mY, -pt.mZ);
}

MPoint operator*(const MPoint& pt, double f)
{
    MPoint result(pt);
    result *= f;
    return result;
}

MPoint operator/(const MPoint& pt, double f)
{
    MPoint result(pt);
    result /= f;
    return result;
}

ostream& operator<<(ostream& os, const MPoint& pt)
{
    os << '(' << pt.mX << ',' << pt.mY << ',' << pt.mZ << ')';
    return os;
}

ostream& operator<<(ostream& os, const vector<MPoint>& pts)
{
    uint32_t n = pts.size();
    os << '[' << n << ']';

            for (const MPoint& pt : pts)
                {
                    os << pt;
                    if (n-- > 1)
                        os << ',';
                }

    return os;
}

// --------------------------------------------------------------------

double DihedralAngle(const MPoint& p1, const MPoint& p2, const MPoint& p3, const MPoint& p4)
{
    MPoint v12 = p1 - p2;	// vector from p2 to p1
    MPoint v43 = p4 - p3;	// vector from p3 to p4

    MPoint z = p2 - p3;		// vector from p3 to p2

    MPoint p = CrossProduct(z, v12);
    MPoint x = CrossProduct(z, v43);
    MPoint y = CrossProduct(z, x);

    double u = DotProduct(x, x);
    double v = DotProduct(y, y);

    double result = 360;
    if (u > 0 and v > 0)
    {
        u = DotProduct(p, x) / sqrt(u);
        v = DotProduct(p, y) / sqrt(v);
        if (u != 0 or v != 0)
            result = atan2(v, u) * 180 / kPI;
    }

    return result;
}

double CosinusAngle(const MPoint& p1, const MPoint& p2, const MPoint& p3, const MPoint& p4)
{
    MPoint v12 = p1 - p2;
    MPoint v34 = p3 - p4;

    double result = 0;

    double x = DotProduct(v12, v12) * DotProduct(v34, v34);
    if (x > 0)
        result = DotProduct(v12, v34) / sqrt(x);

    return result;
}

// --------------------------------------------------------------------


MPoint CenterPoints(vector<MPoint>& points)
{
    MPoint t;

            for (MPoint& pt : points)
                {
                    t.mX += pt.mX;
                    t.mY += pt.mY;
                    t.mZ += pt.mZ;
                }

    t.mX /= points.size();
    t.mY /= points.size();
    t.mZ /= points.size();

            for (MPoint& pt : points)
                {
                    pt.mX -= t.mX;
                    pt.mY -= t.mY;
                    pt.mZ -= t.mZ;
                }

    return t;
}

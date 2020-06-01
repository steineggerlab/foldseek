//
// Created by Martin Steinegger on 28.09.18.
//

#ifndef STRUCCLUST_MAS_H
#define STRUCCLUST_MAS_H


// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#if defined(_MSC_VER)

#ifndef P_WIN
#define P_WIN 1
#endif

// These are Microsoft Visual C++ special settings
// the iso646 file contains the C++ keywords that are
// otherwise not recognized.
#include <ciso646>
#define snprintf _snprintf

// Disable some warnings
#pragma warning (disable : 4996)
#pragma warning (disable : 4355)
#endif

#include <sstream>
#include <iostream>
#include <string>



// Code for amino acid sequences


// 22 real letters and 1 dummy (X is the dummy, B and Z are pseudo letters)
extern const char kResidues[]; // = "ACDEFGHIKLMNPQRSTVWYBZX";
extern const uint8_t kResidueNrTable[];

typedef std::string sequence;


//inline uint8_t ResidueNr(char inAA)
//{
//    int result = 23;
//
//    inAA |= 040;
//    if (inAA >= 'a' and inAA <= 'z')
//        result = kResidueNrTable[inAA - 'a'];
//
//    return result;
//}

//inline bool is_gap(char aa)
//{
//    return aa == ' ' or aa == '.' or aa == '-';
//}

//sequence encode(const std::string& s);
//std::string decode(const sequence& s);



#endif //STRUCCLUST_MAS_H

///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  2007-2008  Institute of Theoretical Physics,Academia Sinica
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this CLePAPS
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this CLePAPS Software into commercial products
or use for commercial purposes should contact Institute of Theoretical Physics,
Academia Sinica,Beijing 100190, China. 
EMAIL: zheng@itp.ac.cn, wangsheng@itp.ac.cn
 
IN NO EVENT SHALL ITP(INSTITUTE OF THEORETICAL PHYSICS) BE LIABLE TO ANY PARTY 
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS CLEPAPS SOFTWARE, EVEN IF ITP HAS 
BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
THE CLEPAPS SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND ITP HAS NO 
OBLIGATION TO PROVIDE MAINTENANCE,SUPPORT, UPDATES, ENHANCEMENTS, OR 
MODIFICATIONS.  ITP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF 
ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE CLEPAPS SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
*/

#pragma once
#include <vector>
#include "Confo_Lett.h"
#include "PDB_File.h"


//====class: PDB_To_CLE====//
//=> transfrom PDB_Format to CLE (or to XYZ_Format)
class PDB_To_CLE : virtual public PDB_File ,virtual public Confo_Lett
{
public:
	PDB_To_CLE(void);
	~PDB_To_CLE(void);

//--------- virtual_function ------------//
public:
	virtual void pdb_btb_ori(int head,int totnum,vector <XYZ> &r,vector <char> &CLE);
};
//====class: PDB_To_CLE====//over

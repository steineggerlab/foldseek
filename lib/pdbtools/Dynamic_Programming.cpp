#include "Dynamic_Programming.h"

//------- dyna_prog init --------//
double DP_M3[4][DYNA_PROG_MAXIMAL];
int DP_D3[4][DYNA_PROG_MAXIMAL];
int DP_align1[DYNA_PROG_LENGTH];
int DP_align2[DYNA_PROG_LENGTH];


//======================= the following is one-layer dynamic programming ====================//__110203__//
int Normal_Align_Dyna_Prog(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
						   vector<pair<int,int> > &alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int VERT  = 1;
	const int HORI  = 2;
	const int DIAG  = 3;
	//create V and D
	double *V=DP_M3[0];
	int *D=DP_D3[0];
	//init()
	V[0*DP_maximal+ 0] = 0;
	D[0*DP_maximal+ 0] = 0;
	for (i = 1; i < m; i++) 
	{
		V[i*DP_maximal+ 0] = V[(i-1)*DP_maximal+ 0] + GAP_HEAD;
		D[i*DP_maximal+ 0] = VERT;
	}

	for (j = 1; j < n; j++) 
	{
		V[0*DP_maximal+ j] = V[0*DP_maximal+ j-1] + GAP_HEAD;
		D[0*DP_maximal+ j] = HORI;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double v1,v2,v3;
	for (i = 1; i < m; i++) 
	{
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = 1; j < n; j++) 
		{
			v1 = V[cur_index1+ j-1] + GAP_OPEN;
			v2 = V[cur_index2+ j] + GAP_OPEN;
			v3 = V[cur_index2+ j-1] + score[cur_index3+ j-1];
			//double dist = distFunc(firstSeq[i-1], secondSeq[j-1]);
			//Point::squaredDistance(atom1, atom2);
			V[cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (V[cur_index1+ j] == v3) D[cur_index1+ j] = DIAG;
			else if (V[cur_index1+ j] == v2) D[cur_index1+ j] = VERT;
			else D[cur_index1+ j] = HORI;
		}
	}
	//build(ali);
	i = m - 1;
	j = n - 1;
	int count = 0;
	int matches = 0;
	int cur_case;
	while ( (i > 0) || (j > 0) ) 
	{
		cur_case=D[i*DP_maximal+ j];
		switch (cur_case)
		{
			case DIAG:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case VERT:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case HORI:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> normal_global: invalid direction D(" << i << ", " << j << ") = " 
				<< D[i*DP_maximal+ j] << endl;
				exit(-1);
		}
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=V[(m-1)*DP_maximal+ n-1];
	return matches;
}

//======================= the following is one-layer dynamic programming (for TM_align) ====================//__110230__//
int Normal_Align_Dyna_Prog_II(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
						      vector<pair<int,int> > &alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int VERT  = 1;
	const int HORI  = 2;
	const int DIAG  = 3;
	//create V and D
	double *V=DP_M3[0];
	int *D=DP_D3[0];
	//init()
	V[0*DP_maximal+ 0] = 0;
	D[0*DP_maximal+ 0] = 0;
	for (i = 1; i < m; i++) 
	{
		V[i*DP_maximal+ 0] = V[(i-1)*DP_maximal+ 0] + GAP_HEAD;
		D[i*DP_maximal+ 0] = VERT;
	}

	for (j = 1; j < n; j++) 
	{
		V[0*DP_maximal+ j] = V[0*DP_maximal+ j-1] + GAP_HEAD;
		D[0*DP_maximal+ j] = HORI;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double v1,v2,v3;
	for (i = 1; i < m; i++) 
	{
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = 1; j < n; j++) 
		{
			v1 = V[cur_index1+ j-1];
			if(D[cur_index1+ j-1]==DIAG) v1 += GAP_OPEN;
			v2 = V[cur_index2+ j];
			if(D[cur_index2+ j]==DIAG) v2 += GAP_OPEN;
			v3 = V[cur_index2+ j-1] + score[cur_index3+ j-1];
			//double dist = distFunc(firstSeq[i-1], secondSeq[j-1]);
			//Point::squaredDistance(atom1, atom2);
			V[cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (V[cur_index1+ j] == v3) D[cur_index1+ j] = DIAG;
			else if (V[cur_index1+ j] == v2) D[cur_index1+ j] = VERT;
			else D[cur_index1+ j] = HORI;
		}
	}
	//build(ali);
	i = m - 1;
	j = n - 1;
	int count = 0;
	int matches = 0;
	int cur_case;
	while ( (i > 0) || (j > 0) ) 
	{
		cur_case=D[i*DP_maximal+ j];
		switch (cur_case)
		{
			case DIAG:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case VERT:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case HORI:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> normal_global_II: invalid direction D(" << i << ", " << j << ") = " 
				<< D[i*DP_maximal+ j] << endl;
				exit(-1);
		}
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=V[(m-1)*DP_maximal+ n-1];
	return matches;
}

//======================= the following is three-layer dynamic programming ====================//__110203__//
int Advance_Align_Dyna_Prog(int n1,int n2,double *score,double GAP_OPEN,double GAP_EXT,
							vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;
	//init()
	double IN_MIN=-1000000;
	DP_D3[_S_][0*DP_maximal+ 0] = -1;
	DP_D3[_H_][0*DP_maximal+ 0] = -1;
	DP_D3[_V_][0*DP_maximal+ 0] = -1;
	DP_M3[_S_][0*DP_maximal+ 0] = 0;
	DP_M3[_H_][0*DP_maximal+ 0] = IN_MIN;
	DP_M3[_V_][0*DP_maximal+ 0] = IN_MIN;
	for (i = 1; i < m; i++) 
	{
		DP_D3[_S_][i*DP_maximal+ 0] = _V_;
		DP_D3[_H_][i*DP_maximal+ 0] = _V_;
		DP_D3[_V_][i*DP_maximal+ 0] = _V_;
		DP_M3[_S_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_H_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_V_][i*DP_maximal+ 0] = 0;//-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		DP_D3[_S_][0*DP_maximal+ j] = _H_;
		DP_D3[_H_][0*DP_maximal+ j] = _H_;
		DP_D3[_V_][0*DP_maximal+ j] = _H_;
		DP_M3[_S_][0*DP_maximal+ j] = IN_MIN;
		DP_M3[_H_][0*DP_maximal+ j] = 0;//-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		DP_M3[_V_][0*DP_maximal+ j] = IN_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = 1; j < n; j++) 
		{
			//condition upper
			gap_open = (j == n-1) ? 0 : GAP_OPEN;
			gap_ext  = (j == n-1) ? 0 : GAP_EXT;
			v1 = DP_M3[_V_][cur_index2+ j] + gap_ext;
			v2 = DP_M3[_S_][cur_index2+ j] + gap_open;
			v3 = DP_M3[_H_][cur_index2+ j] + gap_open;
			DP_M3[_V_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_V_][cur_index1+ j] == v1) DP_D3[_V_][cur_index1+ j] = _V_;
			else if(DP_M3[_V_][cur_index1+ j] == v2) DP_D3[_V_][cur_index1+ j] = _S_;
			else DP_D3[_V_][cur_index1+ j] = _H_;
			//condition left
			gap_open = (i == m-1) ? 0 : GAP_OPEN;
			gap_ext  = (i == m-1) ? 0 : GAP_EXT;
			v1 = DP_M3[_H_][cur_index1+ j-1] + gap_ext;
			v2 = DP_M3[_S_][cur_index1+ j-1] + gap_open;
			v3 = DP_M3[_V_][cur_index1+ j-1] + gap_open;
			DP_M3[_H_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_H_][cur_index1+ j] == v1) DP_D3[_H_][cur_index1+ j] = _H_;
			else if(DP_M3[_H_][cur_index1+ j] == v2) DP_D3[_H_][cur_index1+ j] = _S_;
			else DP_D3[_H_][cur_index1+ j] = _V_;
			//condition diag
			dist = score[cur_index3+ j-1];  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = DP_M3[_V_][cur_index2+ j-1] + dist;
			v2 = DP_M3[_H_][cur_index2+ j-1] + dist;
			v3 = DP_M3[_S_][cur_index2+ j-1] + dist;
			DP_M3[_S_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_S_][cur_index1+ j] == v3) DP_D3[_S_][cur_index1+ j] = _S_;
			else if (DP_M3[_S_][cur_index1+ j] == v1) DP_D3[_S_][cur_index1+ j] = _V_;
			else DP_D3[_S_][cur_index1+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=DP_M3[_V_][i*DP_maximal+ j];
	v2=DP_M3[_H_][i*DP_maximal+ j];
	v3=DP_M3[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	while ( (i > 0) || (j > 0) ) 
	{
		pre_case=DP_D3[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case _V_:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case _H_:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< DP_D3[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=maximal;
	return matches;
}
int Advance_Align_Dyna_Prog_Double(int n1,int n2,double *score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;
	//init()
	double IN_MIN=-1000000;
	DP_D3[_S_][0*DP_maximal+ 0] = -1;
	DP_D3[_H_][0*DP_maximal+ 0] = -1;
	DP_D3[_V_][0*DP_maximal+ 0] = -1;
	DP_M3[_S_][0*DP_maximal+ 0] = 0;
	DP_M3[_H_][0*DP_maximal+ 0] = IN_MIN;
	DP_M3[_V_][0*DP_maximal+ 0] = IN_MIN;
	for (i = 1; i < m; i++) 
	{
		DP_D3[_S_][i*DP_maximal+ 0] = _V_;
		DP_D3[_H_][i*DP_maximal+ 0] = _V_;
		DP_D3[_V_][i*DP_maximal+ 0] = _V_;
		DP_M3[_S_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_H_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		DP_D3[_S_][0*DP_maximal+ j] = _H_;
		DP_D3[_H_][0*DP_maximal+ j] = _H_;
		DP_D3[_V_][0*DP_maximal+ j] = _H_;
		DP_M3[_S_][0*DP_maximal+ j] = IN_MIN;
		DP_M3[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		DP_M3[_V_][0*DP_maximal+ j] = IN_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = DP_M3[_V_][cur_index2+ j] + gap_ext;
			v2 = DP_M3[_S_][cur_index2+ j] + gap_open;
			v3 = DP_M3[_H_][cur_index2+ j] + gap_open;
			DP_M3[_V_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_V_][cur_index1+ j] == v1) DP_D3[_V_][cur_index1+ j] = _V_;
			else if(DP_M3[_V_][cur_index1+ j] == v2) DP_D3[_V_][cur_index1+ j] = _S_;
			else DP_D3[_V_][cur_index1+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = DP_M3[_H_][cur_index1+ j-1] + gap_ext;
			v2 = DP_M3[_S_][cur_index1+ j-1] + gap_open;
			v3 = DP_M3[_V_][cur_index1+ j-1] + gap_open;
			DP_M3[_H_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_H_][cur_index1+ j] == v1) DP_D3[_H_][cur_index1+ j] = _H_;
			else if(DP_M3[_H_][cur_index1+ j] == v2) DP_D3[_H_][cur_index1+ j] = _S_;
			else DP_D3[_H_][cur_index1+ j] = _V_;
			//condition diag
			dist = score[cur_index3+ j-1];  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = DP_M3[_V_][cur_index2+ j-1] + dist;
			v2 = DP_M3[_H_][cur_index2+ j-1] + dist;
			v3 = DP_M3[_S_][cur_index2+ j-1] + dist;
			DP_M3[_S_][cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (DP_M3[_S_][cur_index1+ j] == v3) DP_D3[_S_][cur_index1+ j] = _S_;
			else if (DP_M3[_S_][cur_index1+ j] == v1) DP_D3[_S_][cur_index1+ j] = _V_;
			else DP_D3[_S_][cur_index1+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=DP_M3[_V_][i*DP_maximal+ j];
	v2=DP_M3[_H_][i*DP_maximal+ j];
	v3=DP_M3[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	while ( (i > 0) || (j > 0) ) 
	{
		pre_case=DP_D3[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case _V_:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case _H_:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< DP_D3[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=maximal;
	return matches;
}

//===================== Four_State DynaProg (for TM_align) =======================//__110230__//
int Advance_Align_Dyna_Prog_II(int n1,int n2,double *score,double GAP_OPEN,double GAP_EXT,
							   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;
	const int _D_  = 3;
	//init()
	double IN_MIN=-1000000;
	DP_D3[_S_][0*DP_maximal+ 0] = -1;
	DP_D3[_H_][0*DP_maximal+ 0] = -1;
	DP_D3[_V_][0*DP_maximal+ 0] = -1;
	DP_D3[_D_][0*DP_maximal+ 0] = -1;
	DP_M3[_S_][0*DP_maximal+ 0] = 0;
	DP_M3[_H_][0*DP_maximal+ 0] = IN_MIN;
	DP_M3[_V_][0*DP_maximal+ 0] = IN_MIN;
	DP_M3[_D_][0*DP_maximal+ 0] = 0;
	for (i = 1; i < m; i++) 
	{
		DP_D3[_S_][i*DP_maximal+ 0] = _V_;
		DP_D3[_H_][i*DP_maximal+ 0] = _V_;
		DP_D3[_V_][i*DP_maximal+ 0] = _V_;
		DP_D3[_D_][i*DP_maximal+ 0] = _V_;
		DP_M3[_S_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_H_][i*DP_maximal+ 0] = IN_MIN;
		DP_M3[_V_][i*DP_maximal+ 0] = 0;
		DP_M3[_D_][i*DP_maximal+ 0] = IN_MIN;
	}
	for (j = 1; j < n; j++) 
	{
		DP_D3[_S_][0*DP_maximal+ j] = _H_;
		DP_D3[_H_][0*DP_maximal+ j] = _H_;
		DP_D3[_V_][0*DP_maximal+ j] = _H_;
		DP_D3[_D_][0*DP_maximal+ j] = _H_;
		DP_M3[_S_][0*DP_maximal+ j] = IN_MIN;
		DP_M3[_H_][0*DP_maximal+ j] = 0;	
		DP_M3[_V_][0*DP_maximal+ j] = IN_MIN;
		DP_M3[_D_][0*DP_maximal+ j] = IN_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3,v4;
	double dist;
	for (i = 1; i < m; i++) 
	{
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = 1; j < n; j++) 
		{
			//condition upper
			gap_open = (j == n-1) ? 0 : GAP_OPEN;
			gap_ext  = (j == n-1) ? 0 : GAP_EXT;
			v1 = DP_M3[_V_][cur_index2+ j] + gap_ext;
			v2 = DP_M3[_S_][cur_index2+ j] + gap_open;
			v3 = DP_M3[_H_][cur_index2+ j] + gap_open;
			v4 = DP_M3[_D_][cur_index2+ j] + gap_open;
			DP_M3[_V_][cur_index1+ j] = std::max(v4,std::max(v1, std::max(v2, v3)));
			if (DP_M3[_V_][cur_index1+ j] == v1) DP_D3[_V_][cur_index1+ j] = _V_;
			else if(DP_M3[_V_][cur_index1+ j] == v2) DP_D3[_V_][cur_index1+ j] = _S_;
			else if(DP_M3[_V_][cur_index1+ j] == v4) DP_D3[_V_][cur_index1+ j] = _D_;
			else DP_D3[_V_][cur_index1+ j] = _H_;
			//condition left
			gap_open = (i == m-1) ? 0 : GAP_OPEN;
			gap_ext  = (i == m-1) ? 0 : GAP_EXT;
			v1 = DP_M3[_H_][cur_index1+ j-1] + gap_ext;
			v2 = DP_M3[_S_][cur_index1+ j-1] + gap_open;
			v3 = DP_M3[_V_][cur_index1+ j-1] + gap_open;
			v4 = DP_M3[_D_][cur_index1+ j-1] + gap_open;
			DP_M3[_H_][cur_index1+ j] = std::max(v4,std::max(v1, std::max(v2, v3)));
			if (DP_M3[_H_][cur_index1+ j] == v1) DP_D3[_H_][cur_index1+ j] = _H_;
			else if(DP_M3[_H_][cur_index1+ j] == v2) DP_D3[_H_][cur_index1+ j] = _S_;
			else if(DP_M3[_H_][cur_index1+ j] == v4) DP_D3[_H_][cur_index1+ j] = _D_;
			else DP_D3[_H_][cur_index1+ j] = _V_;
			//condition diag
			dist = score[cur_index3+ j-1];  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = DP_M3[_V_][cur_index2+ j-1] + dist;
			v2 = DP_M3[_H_][cur_index2+ j-1] + dist;
			v3 = DP_M3[_S_][cur_index2+ j-1] + dist;
			v4 = DP_M3[_D_][cur_index2+ j-1] + dist;
			DP_M3[_S_][cur_index1+ j] = std::max(v4,std::max(v1, std::max(v2, v3)));
			if (DP_M3[_S_][cur_index1+ j] == v3) DP_D3[_S_][cur_index1+ j] = _S_;
			else if (DP_M3[_S_][cur_index1+ j] == v4) DP_D3[_S_][cur_index1+ j] = _D_;
			else if (DP_M3[_S_][cur_index1+ j] == v1) DP_D3[_S_][cur_index1+ j] = _V_;
			else DP_D3[_S_][cur_index1+ j] = _H_;
			//condition diag0
			dist = 0.0;
			v1 = DP_M3[_V_][cur_index2+ j-1] + dist;
			v2 = DP_M3[_H_][cur_index2+ j-1] + dist;
			v3 = DP_M3[_S_][cur_index2+ j-1] + dist;
			v4 = DP_M3[_D_][cur_index2+ j-1] + dist;
			DP_M3[_D_][cur_index1+ j] = std::max(v4,std::max(v1, std::max(v2, v3)));
			if (DP_M3[_D_][cur_index1+ j] == v4) DP_D3[_D_][cur_index1+ j] = _D_;
			else if (DP_M3[_D_][cur_index1+ j] == v3) DP_D3[_D_][cur_index1+ j] = _S_;
			else if (DP_M3[_D_][cur_index1+ j] == v1) DP_D3[_D_][cur_index1+ j] = _V_;
			else DP_D3[_D_][cur_index1+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=DP_M3[_V_][i*DP_maximal+ j];
	v2=DP_M3[_H_][i*DP_maximal+ j];
	v3=DP_M3[_S_][i*DP_maximal+ j];
	v4=DP_M3[_D_][i*DP_maximal+ j];
	double maximal = std::max(v4,std::max(v1, std::max(v2, v3)));
	int k = -1;
	if(v4==maximal)k = _D_;
	else if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	while ( (i > 0) || (j > 0) ) 
	{
		pre_case=DP_D3[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case _D_:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case _V_:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case _H_:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< DP_D3[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=maximal;
	return matches;
}

//================ Fast_Dynamic_Programming ================//__110730__//
// [motivation]: given an initial alignment during structural alignment,
// we could only search for a very narrow adjacent area around the alignment.
// That is to say, the searching space for DP could be restricted in O(n),
// as well as the searching time.
// [structure]: we use "bound" to record the additional data structure,
// length is n1, first is n1's correspondence, second is end position
int Normal_Align_Dyna_Prog_Fast(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
								vector<pair<int,int> > &bound,vector<pair<int,int> > &alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//judge
	if(m*n>DYNA_PROG_MAXIMAL)
	{
		fprintf(stderr,"dynamic programming exceed!!\n");
		exit(-1);
	}
	//const value
	const int VERT  = 1;
	const int HORI  = 2;
	const int DIAG  = 3;
	//create V and D
	double *V=DP_M3[0];
	int *D=DP_D3[0];
	//init()
	V[0*DP_maximal+ 0] = 0;
	D[0*DP_maximal+ 0] = 0;
	//first line
	int wsstart,wsend;
	int prev_start,prev_end;
	double v1,v2,v3;
	double IN_MIN=-1000000;
	//real_start
	i=0;
	wsstart=bound[i].first;
	wsend=bound[i].second;
	for (j = 1; j <= wsend; j++) 
	{
		V[i*DP_maximal+ j] = V[i*DP_maximal+ j-1] + GAP_HEAD;
		D[i*DP_maximal+ j] = HORI;
	}
	prev_start=wsstart;
	prev_end=wsend;
	//othere line
	for(i=1;i<m;i++)
	{
		//current bound
		wsstart=bound[i].first;
		wsend=bound[i].second;
		//process
		int cur_index1=i*DP_maximal;
		int cur_index2=(i-1)*DP_maximal;
		int cur_index3=(i-1)*IN_maximal;
		for (j = wsstart; j <= wsend; j++) 
		{
			if(j==wsstart) //[start]
			{
				if(wsstart!=prev_start) //normal condition
				{
					v1=IN_MIN;
					v2 = V[cur_index2+ j];
					if(D[cur_index2+ j]==DIAG) v2 += GAP_OPEN;
					v3 = V[cur_index2+ j-1] + score[cur_index3+ j-1];
				}
				else        //vertical condition
				{
					v1=IN_MIN;
					if(wsstart==0) v2 = V[cur_index2+ j] + GAP_HEAD;
					else
					{
						v2 = V[cur_index2+ j];
						if(D[cur_index2+ j]==DIAG) v2 += GAP_OPEN;
					}
					v3=IN_MIN;
				}
			}
			else if(j>prev_end) //[end]
			{
				if(j==prev_end+1) //normal condition
				{
					v1 = V[cur_index1+ j-1];
					if(D[cur_index1+ j-1]==DIAG) v1 += GAP_OPEN;
					v2=IN_MIN;
					v3 = V[cur_index2+ j-1] + score[cur_index3+ j-1];
				}
				else          //horizontal condition
				{
					v1 = V[cur_index1+ j-1];
					if(D[cur_index1+ j-1]==DIAG) v1 += GAP_OPEN;
					v2=IN_MIN;
					v3=IN_MIN;
				}
			}
			else //[body]
			{
				v1 = V[cur_index1+ j-1];
				if(D[cur_index1+ j-1]==DIAG) v1 += GAP_OPEN;
				v2 = V[cur_index2+ j];
				if(D[cur_index2+ j]==DIAG) v2 += GAP_OPEN;
				v3 = V[cur_index2+ j-1] + score[cur_index3+ j-1];
			}
			//double dist = distFunc(firstSeq[i-1], secondSeq[j-1]);
			//Point::squaredDistance(atom1, atom2);
			V[cur_index1+ j] = std::max(v1, std::max(v2, v3));
			if (V[cur_index1+ j] == v3) D[cur_index1+ j] = DIAG;
			else if (V[cur_index1+ j] == v2) D[cur_index1+ j] = VERT;
			else D[cur_index1+ j] = HORI;
		}
		//next bound
		prev_start=wsstart;
		prev_end=wsend;
	}
	//build(ali);
	i = m - 1;
	j = n - 1;
	int count = 0;
	int matches = 0;
	int cur_case;
	while ( (i > 0) || (j > 0) ) 
	{
		cur_case=D[i*DP_maximal+ j];
		switch (cur_case)
		{
			case DIAG:
				DP_align1[count]=i;
				DP_align2[count]=j;
				i--;
				j--;
				++matches;
				break;
			case VERT:
				DP_align1[count]=i;
				DP_align2[count]=-j;
				i--;
				break;
			case HORI:
				DP_align1[count]=-i;
				DP_align2[count]=j;
				j--;
				break;
			default:
				cout << "ERROR!! -> normal_global_II: invalid direction D(" << i << ", " << j << ") = " 
				<< D[i*DP_maximal+ j] << endl;
				exit(-1);
		}
		count++;
	}
	while (j> 0) DP_align1[count]=-i,DP_align2[count]=j,j--,count++;
	while (i> 0) DP_align1[count]=i,DP_align2[count]=0, i--,count++;
	//reverse alignment
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[count-1-i];
		alignment[i].second=DP_align2[count-1-i];
	}
	ali_sco=V[(m-1)*DP_maximal+ n-1];
	return matches;
}

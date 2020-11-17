#include "Acc_Surface.h"


//====================== Acc_Surface DATA =======================//__130830__//
const int Acc_Surface_maxAcc[21]=
{101,242,171,161,127,190,195,73,193,174,182,200,192,214,127,122,144,259,235,148,-1};
const int Acc_Surface_AA2SUB[26]=
{0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};
const int Acc_Surface_AA1Coding[21]=
{0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18,20};


//--------------------- constructor -----------------//
Acc_Surface::Acc_Surface(int num,int ORDER)
{
	AC_MAXIMAL=num;
	AC_Init(AC_MAXIMAL);
	//icosahedron init
	AC_Init_Data();
	icosahedron_order=ORDER;
	AC_Init_Icosahedron(icosahedron_order);
	//three-code parameter
	ACC_BuryCut = 10.25;
	ACC_ExpCut = 42.9;
}
Acc_Surface::~Acc_Surface(void)
{
	AC_Dele(AC_MAXIMAL);
}

//-------------------- init -------------------//
void Acc_Surface::AC_Init(int AC_MAXIMAL)
{
	//input & output
	NewArray2D(&AC_mol,AC_MAXIMAL,4);
	AC_sidechain=new XYZ[AC_MAXIMAL*11];
	AC_side_rec=new int[AC_MAXIMAL];
	AC_side_num=new int[AC_MAXIMAL];
	AC_output=new int[AC_MAXIMAL];
	AC_normal=new int[AC_MAXIMAL];
	//process temporary
	//[minmaxbox & res_neibor]
	AC_mol_boxmin=new XYZ[AC_MAXIMAL];
	AC_mol_boxmax=new XYZ[AC_MAXIMAL];
	AC_res_neibor_p=new int[AC_MAXIMAL];
}
void Acc_Surface::AC_Dele(int AC_MAXIMAL)
{
	//input & output
	DeleteArray2D(&AC_mol,AC_MAXIMAL);
	delete [] AC_sidechain;
	delete [] AC_side_rec;
	delete [] AC_side_num;
	delete [] AC_output;
	delete [] AC_normal;
	//process temporary
	//[minmaxbox & res_neibor]
	delete [] AC_mol_boxmin;
	delete [] AC_mol_boxmax;
	delete [] AC_res_neibor_p;
}

//-------------------- input & init-------------------//
//[inner data]
void Acc_Surface::AC_Init_Data(void)
{
	AC_RN=1.65;       //atom_N radii   -> 1.65
	AC_RCA=1.87;      //atom_CA radii  -> 1.87
	AC_RC=1.76;       //atom_C radii   -> 1.76
	AC_RO=1.40;       //atom_O radii   -> 1.40
	AC_RSIDE=1.8;     //atom_sidechain -> 1.8
	AC_WATER=1.4;     //atom_water     -> 1.4
	AC_DATA[0]=AC_RN;
	AC_DATA[1]=AC_RCA;
	AC_DATA[2]=AC_RC;
	AC_DATA[3]=AC_RO;
	AC_DATA[4]=AC_RSIDE;
}
//[triangle for icosahedron]
void Acc_Surface::AC_Triangle(XYZ x1, XYZ x2, XYZ x3,long level)
{
	long level1;
	double xnorm;
	double vect[3];
	XYZ x4, x5, x6;
	if (level > 0) 
	{
		level1 = level - 1;
		x4=x1+x2;
		x5=x2+x3;
		x6=x1+x3;
		x4.xyz2double(vect);
		xnorm=Vector_Normalize(vect,3);
		x4.double2xyz(vect);
		x5.xyz2double(vect);
		xnorm=Vector_Normalize(vect,3);
		x5.double2xyz(vect);
		x6.xyz2double(vect);
		xnorm=Vector_Normalize(vect,3);
		x6.double2xyz(vect);
		AC_Triangle(x1, x4, x6, level1);
		AC_Triangle(x4, x5, x6, level1);
		AC_Triangle(x4, x2, x5, level1);
		AC_Triangle(x5, x3, x6, level1);
		return;
	}
	double x[3],y[3],z[3];
	x6=x1+x2+x3;
	x6.xyz2double(vect);
	xnorm=Vector_Normalize(vect,3);
	x6.double2xyz(vect);
	//-> push_back operation here for icosahedron_p
	icosahedron_p.push_back(x6);
	x5=x3-x1;
	x4=x2-x1;
	x5.xyz2double(x);
	x4.xyz2double(y);
	cross(z,x,y,3);
	x6.double2xyz(z);
	x6.xyz2double(vect);
	xnorm=Vector_Normalize(vect,3);
	x6.double2xyz(vect);
	//-> push_back operation here for icosahedron_a
	icosahedron_a.push_back(xnorm / 2.0);
	icosahedron_n++;
}
//[calculate icosahedron]
void Acc_Surface::AC_Init_Icosahedron(int order)
{
	XYZ v[12];
	double a, b;
	long i, j, k, level;
	k = 0;
	a = 0.8506508;   //YVERTEX
	b = 0.5257311;   //ZVERTEX
	for (i = 0; i < 2; i++) 
	{
		a = -a;
		for (j = 0; j < 2; j++,k++) 
		{
			b = -b;
			v[k].X = 0.0;
			v[k].Y = a;
			v[k].Z = b;
			k++;
			v[k].X = b;
			v[k].Y = 0.0;
			v[k].Z = a;
			k++;
			v[k].X = a;
			v[k].Y = b;
			v[k].Z = 0.0;
		}
	}
	icosahedron_n = 0;
	icosahedron_p.clear();
	icosahedron_a.clear();
	level = order;
	/* GET ALL 20 FACES OF ICOSAHEDRON */
	for (i = 0; i <= 9; i++) 
	{   /* FIND INTEGRATION POINTS */
		for (j = i + 1; j <= 10; j++) 
		{
			if(v[i].distance_square(v[j]) < 1.21)
			{
				for (k = j + 1; k <= 11; k++) 
				{
					if((v[i].distance_square(v[k])<1.21) && (v[j].distance_square(v[k])<1.21))
						AC_Triangle(v[i], v[j], v[k], level);
				}
			}
		}
	}
	a = 0.0;
	for (i = 0; i < icosahedron_n; i++)a += icosahedron_a[i];
	a = 12.56637 / a;   // #define FOURPI          12.56637
	for (i = 0; i < icosahedron_n; i++)icosahedron_a[i] *= a;
}
//[input molecular]
//input mol MUST has at least N,CA,C,O,CB five atoms!
//side_tot records the additional sidechain atoms!
void Acc_Surface::AC_Input_Mol(XYZ **mol,int moln,int *side_tot)
{
	int i,j,k;
	int tot;
	//process
	AC_side=0;
	for(i=0;i<moln;i++)
	{
		for(j=0;j<4;j++)AC_mol[i][j]=mol[i][j];
		tot=side_tot[i];
		AC_side_rec[i]=AC_side;
		AC_side_num[i]=tot;
		for(k=0;k<tot;k++)AC_sidechain[AC_side+k]=mol[i][4+k];
		AC_side+=tot;
	}
	AC_moln=moln;
}

//-------------------- process function -------------------//
//[residue level]
inline void Acc_Surface::AC_MinMax(XYZ v,double r,XYZ &vmin,XYZ &vmax)
{
	if (v.X - r < vmin.X)vmin.X = v.X - r;
	if (v.X + r > vmax.X)vmax.X = v.X + r;
	if (v.Y - r < vmin.Y)vmin.Y = v.Y - r;
	if (v.Y + r > vmax.Y)vmax.Y = v.Y + r;
	if (v.Z - r < vmin.Z)vmin.Z = v.Z - r;
	if (v.Z + r > vmax.Z)vmax.Z = v.Z + r;
}
void Acc_Surface::AC_Calc_ResBox(void)
{
	int i,j;
	int pos;
	int tot;
	double radii;
	for(i=0;i<AC_moln;i++)
	{
		//init minmaxbox
		AC_mol_boxmin[i]=INT_MAX_NUM;
		AC_mol_boxmax[i]=INT_MIN_NUM;
		for(j=0;j<4;j++) //Backbone
		{
			radii=AC_DATA[j]+AC_WATER;
			AC_MinMax(AC_mol[i][j],radii,AC_mol_boxmin[i],AC_mol_boxmax[i]);
		}
		if(AC_side>0)
		{
			pos=AC_side_rec[i];
			tot=AC_side_num[i];
			radii=AC_RSIDE+AC_WATER;
			for(j=0;j<tot;j++)AC_MinMax(AC_sidechain[pos+j],radii,AC_mol_boxmin[i],AC_mol_boxmax[i]);
		}
	}
}
int Acc_Surface::AC_Calc_ResNeib(int pos)
{
	int i;
	int neibor=0;
	AC_res_neibor_n=0;
	for(i=0;i<AC_moln;i++)
	{
		if(AC_mol_boxmin[pos].X<AC_mol_boxmax[i].X && 
		   AC_mol_boxmin[pos].Y<AC_mol_boxmax[i].Y && 
		   AC_mol_boxmin[pos].Z<AC_mol_boxmax[i].Z &&
		   AC_mol_boxmax[pos].X>AC_mol_boxmin[i].X &&
		   AC_mol_boxmax[pos].Y>AC_mol_boxmin[i].Y &&
		   AC_mol_boxmax[pos].Z>AC_mol_boxmin[i].Z)
		{
			AC_res_neibor_p[neibor]=i;
			neibor++;
		}
	}
	AC_res_neibor_n=neibor;
	return neibor;
}
//[atomic level]
inline int Acc_Surface::AC_InBox(XYZ v,double r,XYZ vmin,XYZ vmax)
{
	if(v.X-r<vmax.X && v.Y-r<vmax.Y && v.Z-r<vmax.Z &&
	   v.X+r>vmin.X && v.Y+r>vmin.Y && v.Z+r>vmin.Z)return 1;
	else return 0;
}
void Acc_Surface::AC_Calc_AtomNeib_Single(XYZ v,double r,double dist)
{
	//init_judge
	if (dist <= 0.00001) //#define EPS             0.00001
	return;
	//add_neibor
	AC_neibor_p.push_back( v-AC_atom_center );
	AC_neibor_a.push_back( (AC_atom_radii*AC_atom_radii-r*r+dist)/(2*AC_atom_radii) );
	//check_neibor
	if(AC_neibor_a[AC_neibor_n]<AC_neibor_a[0])
	{
		XYZ tmp;
		double tmpa;
		tmp=AC_neibor_p[AC_neibor_n];
		AC_neibor_p[AC_neibor_n]=AC_neibor_p[0];
		AC_neibor_p[0]=tmp;
		tmpa=AC_neibor_a[AC_neibor_n];
		AC_neibor_a[AC_neibor_n]=AC_neibor_a[0];
		AC_neibor_a[0]=tmpa;
	}
	AC_neibor_n++;
}
int Acc_Surface::AC_Calc_AtomNeib(XYZ v,double r)
{
	int i,j;
	int pos;  //res_neibor pos
	int spos; //res_neibor sidechain's pos
	int stot; //res_neibor sidechain's tot
	double radii;
	double dist;
	//init
	AC_neibor_n=0;
	AC_neibor_p.clear();
	AC_neibor_a.clear();
	AC_atom_center=v;
	AC_atom_radii=r;
	//collect nearby residue
	for(i=0;i<AC_res_neibor_n;i++)
	{
		pos=AC_res_neibor_p[i];
		if(AC_InBox(v,r,AC_mol_boxmin[pos],AC_mol_boxmax[pos]))
		{
			for(j=0;j<4;j++) //Backbone
			{
				radii=AC_DATA[j]+AC_WATER;
				dist=v.distance_square(AC_mol[pos][j]);
				if(dist<(r+radii)*(r+radii))AC_Calc_AtomNeib_Single(AC_mol[pos][j],radii,dist);
			}
			if(AC_side>0)  //sidechain
			{
				spos=AC_side_rec[pos];
				stot=AC_side_num[pos];
				radii=AC_RSIDE+AC_WATER;
				for(j=0;j<stot;j++)
				{
					dist=v.distance_square(AC_sidechain[spos+j]);
					if(dist<(r+radii)*(r+radii))AC_Calc_AtomNeib_Single(AC_sidechain[spos+j],radii,dist);
				}
			}
		}
	}
	return AC_neibor_n;
}
double Acc_Surface::AC_Calc_AtomAcc(void)
{
	int i,k;
	int lastk;
	double dist;
	double x[3],y[3];
	double f;
	lastk=0;
	f=0.0;
	for(i=0;i<icosahedron_n;i++) //total icosahedron number
	{
		icosahedron_p[i].xyz2double(x);
		AC_neibor_p[lastk].xyz2double(y);
		dist=dot(x,y,3);
		if(dist<=AC_neibor_a[lastk]) //the lastk CANNOT cover this water!!
		{
			for(k=0;k<AC_neibor_n;k++)
			{
				AC_neibor_p[k].xyz2double(y);
				dist=dot(x,y,3);
				if(dist>AC_neibor_a[k])break; //find a new k which CAN cover!!
			}
			if(k<AC_neibor_n)lastk=k;    //reset this k as the lastk
			else f+=icosahedron_a[i];    //no neibor can cover!! then ADD this water!!
		}
	}
	return AC_atom_radii*AC_atom_radii*f; //scale it with the radius of the current atom
}

//-------------------- major function -------------------//
//input mol SHOULD have at least N,CA,C,O,CB five atoms
void Acc_Surface::AC_Calc_SolvAcc(XYZ **mol,char *ami,int moln,char *acc,int *side_tot)
{
	int i,j;
	int spos,stot;
	double radii;
	double totacc;
	double maxacc;
	double nomacc;
	//input mol
	AC_Input_Mol(mol,moln,side_tot);
	//process
	AC_Calc_ResBox();  //calculate minmaxbox for residue
	for(i=0;i<AC_moln;i++)
	{
		AC_Calc_ResNeib(i); //find residue's neibors
		totacc=0.0;
		for(j=0;j<4;j++)    //Backbone atom
		{
			radii=AC_DATA[j]+AC_WATER;
			AC_Calc_AtomNeib(mol[i][j],radii);
			totacc+=AC_Calc_AtomAcc();
		}
		if(AC_side>0)
		{
			spos=AC_side_rec[i];
			stot=AC_side_num[i];
			radii=AC_RSIDE+AC_WATER;
			for(j=0;j<stot;j++) //sidechain atom
			{
				AC_Calc_AtomNeib(AC_sidechain[spos+j],radii);
				totacc+=AC_Calc_AtomAcc();
			}
		}
		//assign original ACC value
		totacc=(int)floor(totacc+0.5);
		AC_output[i]=(int)totacc;
		//assign normalized ACC value
		maxacc=Acc_Surface_maxAcc[Acc_Surface_AA1Coding[Acc_Surface_AA2SUB[ami[i]-'A']]];
		if(maxacc<0)maxacc=242;
		nomacc=1.0*totacc*100.0/maxacc;
		AC_normal[i]=(int)(nomacc+0.5);
		if(AC_normal[i]>100) AC_normal[i]=100;
		if(AC_normal[i]<0) AC_normal[i]=0;
		//assign ACC code
		if(nomacc<ACC_BuryCut)acc[i]='B';
		else if(nomacc>ACC_ExpCut)acc[i]='E';
		else acc[i]='M';
	}
	acc[AC_moln]='\0';
}

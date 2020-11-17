#include "PDB_Residue.h"
#include "Utility.h"
#include "PDB_Chain.h"
#include "PDB_Utility.h"
#include <iomanip>
using namespace std;


//-------- PDB_Chain_Define ---------// 
PDB_Chain::PDB_Chain()
{
	chain_id='A';
	length=0;
	SSE="";
	CLE="";
	SEQ="";
}
PDB_Chain::~PDB_Chain()
{
}
PDB_Chain::PDB_Chain(const PDB_Chain & chain)
{
	this->chain_id = chain.chain_id;
	this->length = chain.length;
	this->residues = chain.residues;
	this->phi_psi_omegas = chain.phi_psi_omegas;
	this->theta_taus = chain.theta_taus;
	this->SSE.assign(chain.SSE);
	this->CLE.assign(chain.CLE);
	this->SEQ.assign(chain.SEQ);
}
PDB_Chain &  PDB_Chain::operator =(const PDB_Chain & chain)
{
	if(this==&chain)return *this;
	this->chain_id = chain.chain_id;
	this->length = chain.length;
	this->residues = chain.residues;
	this->phi_psi_omegas = chain.phi_psi_omegas;
	this->theta_taus = chain.theta_taus;
	this->SSE.assign(chain.SSE);
	this->CLE.assign(chain.CLE);
	this->SEQ.assign(chain.SEQ);
	return *this;
}

//=================== main process part ====================//
//[init & calc]
//init
int PDB_Chain::initialize(int len, char id)
{
	length = len;
	allocate(len);
	chain_id = id;
	//init
	for( int i = 0; i < length; i++ )
	{
		char buff[7] = "A     ";
		buff[6] = '\0';
		sprintf( buff, "%c%4d ", chain_id, i+1 );
		residues.at(i).set_PDB_residue_number( string( buff ) );
		string s;
		residues.at(i).get_PDB_residue_number( s );
		debug_info( DEBUG_DEFAULT, "residue i: %d  '%s'\n", i, s.c_str() );
	}
	return 0;
}
int PDB_Chain::initialize_simple(int len, char id)
{
	length = len;
	allocate(len);
	chain_id = id;
	return 0;
}
void PDB_Chain::allocate(void)
{
	residues.resize( length );
	phi_psi_omegas.resize( length );
	theta_taus.resize( length );
}
void PDB_Chain::allocate( int n )
{
	length = n;
	allocate();
}

//============ virtual functions =============//
//phi_psi <-> XYZ
int PDB_Chain::refold_from_phi_psi_omega(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::refold_from_theta_tau(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_phi_psi_omega(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_theta_tau(void)
{
	return UNSUPPORT_ERROR;
}
//calculate CLE,SSE and SEQ
int PDB_Chain::calculate_CLE(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_SSE(void)
{
	return UNSUPPORT_ERROR;
}

/*------- Learning -------
//The difference between "at" and operator[] is that, 
//"at" does range checks. Otherwise, they are the same.
-----------------------------*/
//================================= get and set ============================//
//[main_part]
int PDB_Chain::get_length() const
{
	return length;
}
//residue
int PDB_Chain::get_residue(int i, PDB_Residue & residue) const
{
	if(i < 0 || i >= (int)residues.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::get_residue, index = %d, excced %d\n",
				i, residues.size());
		return -1;
	}
	residue = residues.at(i);
	return 0;
}
int PDB_Chain::set_residue(int i, const PDB_Residue & residue)
{
	if(i < 0 || i >= (int)residues.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_residue, index = %d, excced %d\n",
				i, residues.size());
		return -1;
	}
	residues.at(i) = residue;
	return 0;
}
//phi_psi_omega
int PDB_Chain::get_phi_psi_omega(int i, Phi_Psi_Omega & phi_psi_omega ) const
{
	if(i < 0 || i >= (int)phi_psi_omegas.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::get_phi_psi_omega, index = %d, excced %d\n",
				i, phi_psi_omegas.size());
		return -1;
	}
	phi_psi_omega = phi_psi_omegas.at(i);
	return 0;
}
int PDB_Chain::set_phi_psi_omega(int i, const Phi_Psi_Omega & phi_psi_omega)
{
	if(i < 0 || i >= (int)phi_psi_omegas.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_phi_psi_omega, index = %d, excced %d\n",
				i, phi_psi_omegas.size());
		return -1;
	}
	phi_psi_omegas.at(i) = phi_psi_omega;
	return 0;
}
//theta_tau
int PDB_Chain::get_theta_tau(int i, Theta_Tau & theta_tau) const
{
	if(i < 0 || i >= (int)theta_taus.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::get_theta_tau, index = %d, excced %d\n",
				i, theta_taus.size());
		return -1;
	}
	theta_tau = theta_taus.at(i);
	return 0;
}

int PDB_Chain::set_theta_tau(int i, const Theta_Tau & theta_tau)
{
	if(i < 0 || i >= (int)theta_taus.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_theta_tau, index = %d, excced %d\n",
				i, theta_taus.size());
		return -1;
	}
	theta_taus.at(i) = theta_tau;
	return 0;
}

//[vice_part]
string PDB_Chain::get_CLE() const
{
	return CLE;
}

int PDB_Chain::set_CLE(string & cle,int WARNING)
{
	if(WARNING==1)debug_info(DEBUG_WARNING, "@PDB_Chain::set_CLE: You should know what you are doing\n");
	if((int)cle.size() != (int)get_length())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_CLE, cle size %d is not equal with the size of the chain %d", cle.size(), get_length());
		return -1;
	}
	CLE.clear();
	CLE=cle;
	return 0;
}
string PDB_Chain::get_SSE() const
{
	return SSE;
}
int PDB_Chain::set_SSE(string & sse,int WARNING)
{
	if(WARNING==1)debug_info(DEBUG_WARNING, "@PDB_Chain::set_SSE: You should know what you are doing\n");
	if((int)sse.size() != (int)get_length())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_SSE, sse size %d is not equal with the size of the chain %d", sse.size(), get_length());
		return -1;
	}
	SSE.clear();
	SSE=sse;
	return 0;
}
string PDB_Chain::get_SEQ() const
{
	return SEQ;
}
int PDB_Chain::set_SEQ(string & seq,int WARNING)
{
	if(WARNING==1)debug_info(DEBUG_WARNING, "@PDB_Chain::set_SEQ: You should know what you are doing\n");
	if((int)seq.size() != (int)get_length())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_SEQ, seq size %d is not equal with the size of the chain %d", seq.size(), get_length());
		return -1;
	}
	SEQ.clear();
	SEQ=seq;
	return 0;
}
char PDB_Chain::get_chain_id() const
{
	return chain_id;
}
int PDB_Chain::set_chain_id(char id)
{
	chain_id = id;
	return 0;
}
//sequence set & get (they are different with get_SEQ!!!)
string PDB_Chain::get_sequence() const
{
	int size;
	size=(int)residues.size();
	string s;
	s.resize(size);
	for(int i=0;i<size;i++)s.at(i) = residues.at(i).get_AA();
	return s;
}
int PDB_Chain::set_sequence(string & seq)
{
	debug_info(DEBUG_WARNING, "@PDB_Chain::set_sequence: You should know what you are doing\n");
	if(seq.size() != residues.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_sequence, seq size %d is not equal with the size of residues %d", seq.size(), residues.size());
		return -1;
	}
	for(int i = 0; i < (int)seq.size(); i++)
	{
		residues.at(i).set_AA(seq.at(i));
		residues.at(i).PDB_residue_backbone_initialize( seq.at(i) );
		residues.at(i).PDB_residue_sidechain_initialize( seq.at(i) );
	}
	return 0;
}

//================================= printf related ============================//
int PDB_Chain::print_phi_psi_omega(ostream & os)
{
	for(int i = 0; i < (int)this->phi_psi_omegas.size(); i++)
	{
		os<<setw(4)<<i + 1<<" ";
		os<<setw(10)<<phi_psi_omegas.at(i).phi / M_PI * 180<<" ";
		os<<setw(10)<<phi_psi_omegas.at(i).psi / M_PI * 180<<" ";
		os<<setw(10)<<phi_psi_omegas.at(i).omega / M_PI * 180<<" ";
		os<<endl;
	}
	os<<endl;
	return 0;
}
int PDB_Chain::print_phi_psi_omega(vector <vector <double> > & out)
{
	out.clear();
	for(int i = 0; i < (int)this->phi_psi_omegas.size(); i++)
	{
		//print
		vector <double> cur;
		cur.push_back(phi_psi_omegas.at(i).phi);
		cur.push_back(phi_psi_omegas.at(i).psi);
		cur.push_back(phi_psi_omegas.at(i).omega);
		//push
		out.push_back(cur);
	}
	return 0;
}

int PDB_Chain::print_theta_tau(ostream & os)
{
	for(int i = 0; i < (int)this->theta_taus.size(); i++)
	{
		os<<setw(4)<<i + 1<<" ";
		os<<setw(10)<<theta_taus.at(i).theta / M_PI * 180<<" ";
		os<<setw(10)<<theta_taus.at(i).tau / M_PI * 180<<" ";
		os<<endl;
	}
	os<<endl;
	return 0;
}
int PDB_Chain::print_theta_tau(vector <vector <double> > & out)
{
	out.clear();
	for(int i = 0; i < (int)this->theta_taus.size(); i++)
	{
		//print
		vector <double> cur;
		cur.push_back(theta_taus.at(i).theta);
		cur.push_back(theta_taus.at(i).tau);
		//push
		out.push_back(cur);
	}
	return 0;
}

int PDB_Chain::print_chain(ostream & os)
{
	//init
	int j,k;
	int length;
	string TER="TER                                                                             ";
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	XYZ xyz;
	char fp[1024];
	int atom_serial_number;
	//process
	length=get_length();
	for(j=0;j<length;j++)
	{
		 //init
		PDB_Residue & PDB = residues.at(j);
		//chains[i].get_residue(j, PDB);
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind);
		//backbone_out
		for(k=0;k<4;k++)
		{
			//check
			if(PDB.get_backbone_part_index(k)==0)
			{
				debug_info(DEBUG_WARNING, "@PDB_Chain::print_chain, get_backbone_part_index[%d]==0\n", k);
				continue;
			}
			//get
			atomname=backbone_atom_name_decode(k);
			PDB.get_backbone_atom(k, xyz, atom_serial_number, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
				atom_serial_number,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			os<<fp;
		}			
		//sidechain_out
		number=PDB.get_sidechain_totnum();
		for(k=0;k<number;k++)
		{
			//check
			if(PDB.get_sidechain_part_index(k)==0)
			{
				debug_info(DEBUG_WARNING, "@PDB_Chain::print_chain, get_sidechain_part_index[%d]==0\n", k);
				continue;
			}
			//get
			atomname=sidechain_atom_name_decode(k,amino);
			PDB.get_sidechain_atom(k, xyz, atom_serial_number, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
				atom_serial_number,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			os<<fp;
		}
	}
	sprintf(fp, "%s\n",TER.c_str());
	os<<fp;
	return 0;
}

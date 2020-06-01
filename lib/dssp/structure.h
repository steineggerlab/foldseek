// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <vector>
#include <limits.h>
#include <string>


#include "primitives3d.h"

struct MAtom;
class MResidue;
class MChain;
class MProtein;



// forward declaration of buffer
template<typename T, uint32_t N> class buffer;
typedef buffer<MResidue*,100>	MResidueQueue;

const uint32_t kHistogramSize = 30;

// a limited set of known atoms. This is an obvious candidate for improvement of DSSP.
enum MAtomType
{
    kUnknownAtom,
    kHydrogen,
    // ...
            kCarbon,
    kNitrogen,
    kOxygen,
    kFluorine,
    // ...
            kPhosphorus,
    kSulfur,
    kChlorine,
    kMagnesium,
    kPotassium,
    kCalcium,
    kZinc,
    kSelenium,

    kAtomTypeCount
};

MAtomType MapElement(std::string inElement);

// for now, MAtom contains exactly what the ATOM line contains in a PDB file
struct MAtom
{
    uint32_t		mSerial;
    std::string	mName;
    char		mAltLoc;
    std::string	mResName;
    std::string	mChainID;
    int16_t		mResSeq;
    std::string	mICode;
    MAtomType	mType;
    MPoint		mLoc;
    double		mOccupancy;
    double		mTempFactor;
    std::string	mElement;
    int			mCharge;

    void		SetChainID(const std::string& inChainID){ mChainID = inChainID;}
    std::string	GetName() const							{ return mName; }
    void		Translate(const MPoint& inTranslation)	{ mLoc += inTranslation; }
    void		Rotate(const MQuaternion& inRotation)	{ mLoc.Rotate(inRotation); }
    void		WritePDB(std::ostream& os) const;

    operator const MPoint&() const			{ return mLoc; }
    operator MPoint&()						{ return mLoc; }
};

enum MResidueType
{
    kUnknownResidue,

    //
            kAlanine,				// A	ala
    kArginine,				// R	arg
    kAsparagine,			// N	asn
    kAsparticAcid,			// D	asp
    kCysteine,				// C	cys
    kGlutamicAcid,			// E	glu
    kGlutamine,				// Q	gln
    kGlycine,				// G	gly
    kHistidine,				// H	his
    kIsoleucine,			// I	ile
    kLeucine,				// L	leu
    kLysine,				// K	lys
    kMethionine,			// M	met
    kPhenylalanine,			// F	phe
    kProline,				// P	pro
    kSerine,				// S	ser
    kThreonine,				// T	thr
    kTryptophan,			// W	trp
    kTyrosine,				// Y	tyr
    kValine,				// V	val

    kResidueTypeCount
};

struct MResidueInfo
{
    MResidueType		type;
    char				code;
    char				name[4];
};

// a residue number to info mapping
extern const MResidueInfo kResidueInfo[];

MResidueType MapResidue(std::string inName);

struct HBond
{
    MResidue*		residue;
    double			energy;
};

enum MBridgeType
{
    btNoBridge, btParallel, btAntiParallel
};

struct MBridgeParner
{
    MResidue*		residue;
    uint32_t			ladder;
    bool			parallel;
};

enum MHelixFlag
{
    helixNone, helixStart, helixEnd, helixStartAndEnd, helixMiddle
};

enum MSecondaryStructure
{
    loop,		//' '
    alphahelix,	// H
    betabridge, // B
    strand,		// E
    helix_3,	// G
    helix_5,	// I
    turn,		// T
    bend		// S
};

class MResidue
{
public:
    MResidue(const MResidue& residue);
    MResidue(uint32_t inNumber, char inTypeCode, MResidue* inPrevious);
    MResidue(uint32_t inNumber,
             MResidue* inPrevious, const std::vector<MAtom>& inAtoms);

    void				SetChainID(const std::string& inChainID);
    std::string			GetChainID() const				{ return mChainID; }

    MResidueType		GetType() const					{ return mType; }

    const MAtom&		GetCAlpha() const				{ return mCA; }
    const MAtom&		GetC() const					{ return mC; }
    const MAtom&		GetN() const					{ return mN; }
    const MAtom&		GetO() const					{ return mO; }
    const MAtom&		GetH() const					{ return mH; }

    double				Phi() const;
    double				Psi() const;
    std::pair<double, char> Alpha() const;
    double				Kappa() const;
    double				TCO() const;

    double				Accessibility() const			{ return mAccessibility; }

    void				SetSecondaryStructure(MSecondaryStructure inSS)
    { mSecondaryStructure = inSS; }
    MSecondaryStructure	GetSecondaryStructure() const	{ return mSecondaryStructure; }

    const MResidue*		Next() const					{ return mNext; }
    const MResidue*		Prev() const					{ return mPrev; }

    void				SetPrev(MResidue* inResidue);

    void				SetBetaPartner(uint32_t n, MResidue* inResidue, uint32_t inLadder,
                                       bool inParallel);
    MBridgeParner		GetBetaPartner(uint32_t n) const;

    void				SetSheet(uint32_t inSheet)	{ mSheet = inSheet; }
    uint32_t				GetSheet() const			{ return mSheet; }

    bool				IsBend() const				{ return mBend; }
    void				SetBend(bool inBend)		{ mBend = inBend; }

    MHelixFlag			GetHelixFlag(uint32_t inHelixStride) const;
    bool				IsHelixStart(uint32_t inHelixStride) const;
    void				SetHelixFlag(uint32_t inHelixStride, MHelixFlag inHelixFlag);

    void				SetSSBridgeNr(uint8_t inBridgeNr);
    uint8_t				GetSSBridgeNr() const;

    void				AddAtom(MAtom& inAtom);

    HBond*				Donor()						{ return mHBondDonor; }
    HBond*				Acceptor()					{ return mHBondAcceptor; }

    const HBond*		Donor() const				{ return mHBondDonor; }
    const HBond*		Acceptor() const			{ return mHBondAcceptor; }

    bool				ValidDistance(const MResidue& inNext) const;

    static bool			TestBond(const MResidue* a, const MResidue* b)
    {
        return a->TestBond(b);
    }

    // bridge functions
    MBridgeType			TestBridge(MResidue* inResidue) const;

    uint16_t				GetSeqNumber() const		{ return mSeqNumber; }
    std::string			GetInsertionCode() const	{ return mInsertionCode; }

    void				SetNumber(uint16_t inNumber)	{ mNumber = inNumber; }
    uint16_t				GetNumber() const			{ return mNumber; }

    void				Translate(const MPoint& inTranslation);
    void				Rotate(const MQuaternion& inRotation);

    void				WritePDB(std::ostream& os);

    static double		CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor);

    std::vector<MAtom>&	GetSideChain()				{ return mSideChain; }
    const std::vector<MAtom>&
    GetSideChain() const		{ return mSideChain; }

    void				GetPoints(std::vector<MPoint>& outPoints) const;

    void				CalculateSurface(const std::vector<MResidue*>& inResidues);

    void				GetCenterAndRadius(MPoint& outCenter, double& outRadius) const
    { outCenter = mCenter; outRadius = mRadius; }

    static bool			NoChainBreak(const MResidue* from, const MResidue* to);

protected:

    double				CalculateSurface(
            const MAtom& inAtom, double inRadius,
            const std::vector<MResidue*>& inResidues);

    bool				TestBond(const MResidue* other) const;

    void				ExtendBox(const MAtom& atom, double inRadius);
    bool				AtomIntersectsBox(const MAtom& atom, double inRadius) const;

    std::string			mChainID;
    MResidue*			mPrev;
    MResidue*			mNext;
    int32_t				mSeqNumber, mNumber;
    std::string			mInsertionCode;
    MResidueType		mType;
    uint8_t				mSSBridgeNr;
    double				mAccessibility;
    MSecondaryStructure	mSecondaryStructure;
    MAtom				mC, mN, mCA, mO, mH;
    HBond				mHBondDonor[2], mHBondAcceptor[2];
    std::vector<MAtom>	mSideChain;
    MBridgeParner		mBetaPartner[2];
    uint32_t				mSheet;
    MHelixFlag			mHelixFlags[3];	//
    bool				mBend;
    MPoint				mBox[2];		// The 3D box containing all atoms
    MPoint				mCenter;		// and the 3d Sphere containing all atoms
    double				mRadius;

private:
    MResidue&			operator=(const MResidue& residue);
};

class MChain
{
public:

    MChain(const MChain& chain);
    MChain(const std::string& inChainID) : mChainID(inChainID) {}
    ~MChain();

    MChain&				operator=(const MChain& chain);

    std::string			GetChainID() const					{ return mChainID; }
    void				SetChainID(const std::string& inChainID);

    MResidue*			GetResidueBySeqNumber(uint16_t inSeqNumber, const std::string& inInsertionCode);

    void				GetSequence(std::string& outSequence) const;

    void				Translate(const MPoint& inTranslation);
    void				Rotate(const MQuaternion& inRotation);

    void				WritePDB(std::ostream& os);

    std::vector<MResidue*>&
    GetResidues()						{ return mResidues; }
    const std::vector<MResidue*>&
    GetResidues() const					{ return mResidues; }

    bool				Empty() const						{ return mResidues.empty(); }

private:
    std::string			mChainID;
    std::vector<MResidue*>
            mResidues;
};

class MProtein
{
public:
    MProtein();
    MProtein(const std::string& inID, MChain* inChain);
    ~MProtein();

//						MProtein(std::istream& is, bool inCAlphaOnly = false);

    void				ReadPDB(std::istream& is, bool inCAlphaOnly = false);
    void				ReadmmCIF(std::istream& is, bool inCAlphaOnly = false);

    const std::string&	GetID() const					{ return mID; }
    const std::string&	GetHeader() const				{ return mHeader; }
    std::string			GetCompound() const;
    std::string			GetSource() const;
    std::string			GetAuthor() const;
    const std::vector<std::string>&
    GetDbRef() const				{ return mDbRef; }

    void				CalculateSecondaryStructure(bool inPreferPiHelices = true);

    void				GetStatistics(uint32_t& outNrOfResidues, uint32_t& outNrOfChains,
                                      uint32_t& outNrOfSSBridges, uint32_t& outNrOfIntraChainSSBridges,
                                      uint32_t& outNrOfHBonds, uint32_t outNrOfHBondsPerDistance[11]) const;

    void				GetCAlphaLocations(const std::string& inChainID, std::vector<MPoint>& outPoints) const;
    MPoint				GetCAlphaPosition(const std::string& inChainID, int16_t inPDBResSeq) const;

    void				Center();
    void				Translate(const MPoint& inTranslation);
    void				Rotate(const MQuaternion& inRotation);

    void				WritePDB(std::ostream& os);

    void				GetPoints(std::vector<MPoint>& outPoints) const;

    std::string			GetFirstChainID() const								{ return mChains.front()->GetChainID(); }

    void				SetChain(const std::string& inChainID, const MChain& inChain);

    MChain&				GetChain(const std::string& inChainID);
    const MChain&		GetChain(const std::string& inChainID) const;

    const std::vector<MChain*>&
    GetChains() const									{ return mChains; }

    template<class OutputIterator>
    void				GetSequences(OutputIterator outSequences) const;

    MResidue*			GetResidue(const std::string& inChainID, uint16_t inSeqNumber, const std::string& inInsertionCode);

    // statistics
    uint32_t				GetNrOfHBondsInParallelBridges() const				{ return mNrOfHBondsInParallelBridges; }
    uint32_t				GetNrOfHBondsInAntiparallelBridges() const			{ return mNrOfHBondsInAntiparallelBridges; }

    void				GetResiduesPerAlphaHelixHistogram(uint32_t outHistogram[30]) const;
    void				GetParallelBridgesPerLadderHistogram(uint32_t outHistogram[30]) const;
    void				GetAntiparallelBridgesPerLadderHistogram(uint32_t outHistogram[30]) const;
    void				GetLaddersPerSheetHistogram(uint32_t outHistogram[30]) const;

private:

    void				AddResidue(const std::vector<MAtom>& inAtoms);

    void				CalculateHBondEnergies(const std::vector<MResidue*>& inResidues);
    void				CalculateAlphaHelices(const std::vector<MResidue*>& inResidues, bool inPreferPiHelices);
    void				CalculateBetaSheets(const std::vector<MResidue*>& inResidues);
    void				CalculateAccessibilities(const std::vector<MResidue*>& inResidues);

    // a thread entry point
    void				CalculateAccessibility(MResidueQueue& inQueue,
                                               const std::vector<MResidue*>& inResidues);

    std::string			mID, mHeader;

    std::vector<std::string>
            mDbRef;
    std::string			mCompound, mSource, mAuthor;
    std::vector<MChain*>mChains;
    uint32_t				mResidueCount, mChainBreaks;

    std::vector<std::pair<MResidue*,MResidue*> >
            mSSBonds;
    uint32_t				mIgnoredWaterMolecules;

    // statistics
    uint32_t				mNrOfHBondsInParallelBridges, mNrOfHBondsInAntiparallelBridges;
    uint32_t				mParallelBridgesPerLadderHistogram[kHistogramSize];
    uint32_t				mAntiparallelBridgesPerLadderHistogram[kHistogramSize];
    uint32_t				mLaddersPerSheetHistogram[kHistogramSize];
};

// inlines

// GetSequences can be used to quickly get all sequences in a vector<string> e.g.
template<class OutputIterator>
void MProtein::GetSequences(OutputIterator outSequences) const
{
    for (std::vector<MChain*>::const_iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
    {
        std::string seq;
        (*chain)->GetSequence(seq);
        *outSequences++ = seq;
    }
}

#endif
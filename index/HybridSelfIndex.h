#ifndef HYBRIDSELFINDEX_H_
#define HYBRIDSELFINDEX_H_

#include <sdsl/suffix_arrays.hpp>
#include <RMQRMM64.h>
#include "wt_fbb.hpp"

using namespace sdsl;
using namespace std;
using namespace rmqrmm;

#ifndef SIGMA
#define SIGMA 256
#endif
#define INI_SIZE 1000 			// Initial size of the array occ[ ] to locate symbols (string of length 1).

// Sampling values for the FMI on the filtered text.
#define S_SA 32
#define S_ISA 1073741824		// If you not desire to display symbols, then set S_ISA with a larger power of 2

#define SAMP_GC 16		// sampling size for Tables with Gap-Encoding for phrases in the filtered Text. We recommend to use a small value. 16 looks very well in the most of cases.
#define N_LOOKUP 128	// length of the Lookup Table for X. (the final length will be in N_LOOKUP/2 to N_LOOKUP)
#define N_SCAN 256		// when a range of X contains more than N_SCAN items we perform a binary search in that segment; otherwise we only SCAN it to find the predecessor

class HybridSelfIndex {
	// typedef structure used at query time, exclusively when m > 2*M, for searches of primary occurrences of p_{0..m-1}.
	typedef struct PriOcc{
		ulong *L;
		ulong n;
		ulong *Phr;
		ulong *Pos;
		uint *Dx;
		bool *B;
	}PriOcc;
private:
	static uint SAMPGC;
	static uint POT_GC;
	static uchar MIN_CHAR;
	static uchar REPLACE_CHAR;	// WE REPLACE ANY SIMBOL WHOSE ASCII CODE <= MIN_CHAR BY 'REPLACE_CHAR'

	//========================================================================
	// Original Parser (these are not stored only and are USED ONLY IN CONDTRUCTION TIME)
	ulong *ARR_POS;		// with lgPOS bits per cell
	uint lgPOS;
	ulong *ARR_LEN;		// with lgLEN bits per cell
	uint lgLEN;

	//========================================================================
	// Structures to mapping from/to the phrase boundaries between Txt and fTxt (Tables stored with Gap-Encoding):

	bit_vector_il<512> BL_il;
	bit_vector_il<512>::rank_1_type rankBL;

	uint nSamP;			// numbers of sampled positions of PhraT and PhraFT

	ulong *PhraT;		// Table of length z that stores the positions of Txt of each phrase boundary
	uint lgPT;			// bits for each gap in PhraT
	ulong *SGCPT;		// Table of length nSamP for the sampling of PhraT

	ulong *PhraFT;		// Table of length z that stores the positions of fTxt of each phrase boundary
	uint lgPFT;			// bits for each gap in PhraFT
	ulong *SGCPFT;		// Table of length nSamP for the sampling of PhraFT

	//========================================================================
	uint nG;			// number of points in the new Grid
	uint lgNG;			// log(nG) bits per cell in the Array Src[]

	ulong *X;			// increasing sequences for position of sources
	uint b_X;			// power of 2 of the wide block for lookup table for X
	uint *LOOKUP_X;		// lookup Table for X
	uint *COUNT_X;		// counters associated to LOOKUP_X
	uint nSampX;		// this is the final length of the tables LOOKUP_X and COUNT_X (N_LOOKUP/2 <= nSampX <= N_LOOKUP)
	RMQRMM64 *rmqY;		// Range Maximal Query structure on Y coordinate of sources

	ulong *Pr;			// Table (of length nG) of pointers to the phrase identifier, stored explicitly in nG*lgZ bits

	ulong *Src;			// Array of length nG that stores the position in G of the source for each phrase, stored explicitly in nG*lgNG bits

	//========================================================================
	csa_wt<wt_fbb<hyb_vector<> >, S_SA, S_ISA> FMI;  			// FMI of the filtered text

	//========================================================================

	void swimInMinQX(ulong *Q, ulong newX, ulong i, ulong k, ulong *Source);
	void setTopMinQX(ulong *Q, uint k, ulong *Source);

	void buildBasicStructures(ulong **Source);
	void buildPtrRMQonG(ulong *Source);

	void buildFMI();
	void readParsing(char *inputFile);
	void checkSymbols(uchar *seq);

	void testPredecessor();
	bool findPredecesor(ulong x, ulong *pred);
	bool isPrimaryLeftLarge(ulong x, uint len, ulong *phra, uint *dx);
	bool isPrimaryLeft_b(ulong x, uint len, ulong *pIni, uint *dIni);
	bool isPrimaryMidle(ulong x, uint len, ulong *phra, uint *dx);
	bool isPrimaryRight(ulong x, uint len, ulong *phra, uint *dx);

	ulong getPosPhraT(ulong phra);
	ulong getPosPhraFT(ulong phra);
	bool isPrimary(ulong x, uint len, ulong *pIni, uint *dIni);
	void locatePryOcc(uchar *pat, uint m, ulong *nOcc, ulong **L, ulong* currN);
	void locateSecOcc(ulong l, ulong r, ulong posX, uint m, ulong *nOcc, ulong **occ, ulong *currN);
	ulong searchPhraFilTxt(ulong x, uint *dx);
	ulong searchPhraTxt(ulong x, ulong *pos, uint *len);

	void locateAChar(uchar *pat, ulong *nOcc, ulong **occ);
	void locateUptoM(uchar *pat, uint m, ulong *nOcc, ulong **occ);
	void setTopMinQ(ulong *Q, ulong pos);
	bool isOccInArr(ulong u, ulong len, ulong *A);
	void locateUpto2M(uchar *pat, uint m, ulong *nOcc, ulong **occ);

	bool isOccInArr(ulong u, ulong len, ulong *A, ulong *pos);
	void locateLargeM_b(uchar *pat, uint m, ulong *nOcc, ulong **occ);

	void extract(ulong sp, ulong len, uchar *A);

public:
	static bool TRACE;		// true: prints all details for stdout
	static bool SHOW_SIZE;	// true: prints the breakdown size of the index
	static bool CREATE_FMI_TEST;	// true: create the TMI_TEST of the original test

	char dirStore[300];	// directory to save/load the data structure

	ulong sizeDS;		// size in bytes for the complete data structure

	uchar fSymb;		// mismatch symbol
	uchar *Txt;			// original text of length n (used in construction time)
	ulong n;			// length of the concatenated text T = d_1+d_2+ .. + d_D

	uchar *fTxt;		// filtered text of length nFT
	ulong nFT;			// length of the filtered text

	uint lgN;
	uint M;				// optimal pattern length
	uint lgM;

	ulong orig_z;		// original number of LZ77 phrases
	uint z;				// number of LZ77 phrases after reduction
	uint lgZ;

	// FMI of the original text. It is created only if it want to test the index, that is when CREATE_FMI_TEST=true
	csa_wt<wt_fbb<hyb_vector<> >, 32, 32> FMI_TEST;

	// parserFile: filename of the LZ77 parser
	// maxM: M value for the structure
	// dirSaveLoad: folder where save/load the structure
	HybridSelfIndex(char *parserFile, uint optM, char dirSaveLoad[300]);

	// load the structure from dirSaveLoad
	HybridSelfIndex(char dirSaveLoad[300]);

	// writes in *nOcc the number of occurrences of the pattern *pat[0..m-1] allocating these in **occ.
	void locate(uchar *pat, uint m, ulong *nOcc, ulong **occ);

	// extracts the 'len' characters that represent the original segment T[sp.. sp+len-1] and allocates these into A[0...len-1]
	void extract(ulong sp, ulong len, uchar **A);

	// used for test !
	bool isPrimaryT(ulong x, uint m, ulong *phra);

	// save load structure to/from the base-name 'pathFile', typically a directory with a prefix.
	void saveStructure();
	void loadStructure();

	virtual ~HybridSelfIndex();
};

#endif /* HYBRIDSELFINDEX_H_ */

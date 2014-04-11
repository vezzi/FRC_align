#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <sstream>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"


#ifdef INLINE_DISABLED
#define INLINE
#else
#define INLINE inline
#endif


using namespace BamTools;
using namespace std;

#define DEFAULT_CHANNEL std::cout
#define ERROR_CHANNEL std::cerr

#define VERBOSE_CHANNEL std::cerr
#define DEBUG_CHANNEL std::cerr
#define DEFAULT_CHANNEL std::cout

static inline std::string package_description() {
	std::string line("FRC");
	line.append(" version ");
	line.append("1.2.0");
	return line;
}

enum readStatus {unmapped, lowQualty, singleton, pair_wrongChrs,
	pair_proper, pair_wrongDistance, pair_wrongOrientation};


enum Feature {LOW_COVERAGE_AREA, HIGH_COVERAGE_AREA, LOW_NORMAL_AREA, HIGH_NORMAL_AREA, HIGH_SINGLE_AREA, HIGH_SPANNING_AREA, HIGH_OUTIE_AREA, COMPRESSION_AREA, STRECH_AREA, TOTAL};


enum FeatureTypes {
	FRC_TOTAL,
	LOW_COV_PE,
	HIGH_COV_PE,
	LOW_NORM_COV_PE,
	HIGH_NORM_COV_PE,
	HIGH_SINGLE_PE,
	HIGH_SPAN_PE,
	HIGH_OUTIE_PE,
	COMPR_PE,
	STRECH_PE,
	HIGH_SINGLE_MP,
	HIGH_OUTIE_MP,
	HIGH_SPAN_MP,
	COMPR_MP,
	STRECH_MP,
};



static int StringToNumber ( string Text ) {
	stringstream ss(Text);
	int result;
	return ss >> result ? result : 0;
}


struct LibraryStatistics{
	float C_A;
	float S_A;
	float C_D;
	float C_M;
	float C_S;
	float C_W;
	float insertMean;
	float insertStd;
};






static readStatus computeReadType(BamAlignment al, uint32_t max_insert, bool is_mp) {
	if (!al.IsMapped()) {
		return unmapped;
	}
	if((al.IsDuplicate()) || (al.IsFailedQC())) {
		return lowQualty;
	}
	if(!(al.IsMateMapped())) {
		return singleton;
	}
	if (al.IsMateMapped() && al.RefID != al.MateRefID) {
		return pair_wrongChrs;
	}
	//If I am here the read must be aligned, with the pair/mate aligned on the same contig/scaffold
	uint32_t startRead   = al.Position; // start position on the contig
	uint32_t startPaired = al.MatePosition;
	int iSize = al.InsertSize;
	if (iSize < 0) { iSize = -1 * iSize;}
	//Now check if reads belong to a proper pair: both reads aligned on the same contig at the expected distance and orientation
	if (iSize > max_insert) {
		return pair_wrongDistance;
	}
	if (! is_mp) { // I have a paired end
		if(startRead < startPaired) { //
			if(!(al.IsReverseStrand()) && (al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		} else {
			if((al.IsReverseStrand()) && !(al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		}
	} else {
		if(startRead < startPaired) { //
			if((al.IsReverseStrand()) && !(al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		} else {
			if(!(al.IsReverseStrand()) && (al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		}
	}
}









#endif /*TYPES_H_*/

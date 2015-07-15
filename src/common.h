#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <sstream>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include <boost/filesystem.hpp>


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
	line.append("1.3.0");
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

static string returnFeatureName(FeatureTypes type) {
	switch (type) {
	case FRC_TOTAL:
		return "TOTAL";
		break;
	case LOW_COV_PE:
		return "LOW_COV_PE";
		break;
	case HIGH_COV_PE:
		return "HIGH_COV_PE";
		break;
	case LOW_NORM_COV_PE:
		return "LOW_NORM_COV_PE";
		break;
	case HIGH_NORM_COV_PE:
		return "HIGH_NORM_COV_PE";
		break;
	case HIGH_SINGLE_PE:
		return "HIGH_SINGLE_PE";
		break;
	case HIGH_SPAN_PE:
		return "HIGH_SPAN_PE";
		break;
	case HIGH_OUTIE_PE:
		return "HIGH_OUTIE_PE";
		break;
	case COMPR_PE:
		return "COMPR_PE";
		break;
	case STRECH_PE:
		return "STRECH_PE";
		break;
	case HIGH_SINGLE_MP:
		return "HIGH_SINGLE_MP";
		break;
	case HIGH_OUTIE_MP:
		return "HIGH_OUTIE_MP";
		break;
	case HIGH_SPAN_MP:
		return "HIGH_SPAN_MP";
		break;
	case COMPR_MP:
		return "COMPR_MP";
		break;
	case STRECH_MP:
		return "STRECH_MP";
		break;
	default:
		cout << "THis whould never happen\n";
	}

}

static int StringToNumber ( string Text ) {
	stringstream ss(Text);
	int result;
	return ss >> result ? result : 0;
}


struct LibraryStatistics{
	uint32_t reads;
	uint32_t mappedReads;
	uint32_t unmappedReads;
	uint32_t matedReads;
	uint32_t wrongDistanceReads;
	uint32_t lowQualityReads;
	uint32_t wronglyOrientedReads;
	uint32_t matedDifferentContig;
	uint32_t singletonReads;

	float C_A;
	float S_A;
	float C_D;
	float C_M;
	float C_S;
	float C_W;
	float insertMean;
	float insertStd;
	string library_name;
};






static readStatus computeReadType(BamAlignment al, uint32_t max_insert, bool is_mp) {
	if (!al.IsMapped()) {
		return unmapped;
	}
	if((al.IsDuplicate()) || (al.IsFailedQC()) || (!al.IsPrimaryAlignment())) {
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



static LibraryStatistics computeLibraryStats(string bamFileName, uint64_t genomeLength, uint32_t max_insert, bool is_mp) {
	BamReader bamFile;
	bamFile.Open(bamFileName);
	LibraryStatistics library;
	string library_name = boost::filesystem::path(bamFileName).stem().string();
	library.library_name = library_name;

	//All var declarations
	uint32_t reads 		 	  = 0;
	uint32_t unmappedReads 	  = 0;
	uint32_t lowQualityReads  = 0;
	uint32_t mappedReads 	  = 0;
	uint64_t mappedReadsLength= 0;

	uint64_t insertsLength = 0; // total inserts length
	float insertMean;
	float insertStd;
	// mated reads (not necessary correctly mated)
	uint32_t matedReads 	  = 0;        // reads that align on a contig with the mate
	uint64_t matedReadsLength = 0;  // total length of mated reads
	// wrongly distance
	uint32_t wrongDistanceReads 		= 0;  // number of paired reads too far away
	uint64_t wrongDistanceReadsLength   = 0; // length  of paired reads too far away
	// wrongly oriented reads
	uint32_t wronglyOrientedReads 		= 0;       // number of wrongly oriented reads
	uint64_t wronglyOrientedReadsLength = 0; // length of wrongly oriented reads
	// singletons
	uint32_t singletonReads 	  = 0; // number of singleton reads
	uint64_t singletonReadsLength = 0;     // total length of singleton reads
	// mates on different contigs
	uint32_t matedDifferentContig 		= 0; // number of contig placed in a different contig
	uint64_t matedDifferentContigLength = 0; // total number of reads placed in different contigs

	float C_A = 0; // total read coverage
	float S_A = 0; // total span coverage
	float C_M = 0; // coverage induced by proper pairs (same contig and correct orientation)
	float C_W = 0; // coverage induced by wrongly mated pairs
	float C_S = 0; // coverage induced by singletons
	float C_D = 0; // coverage induced by reads with mate on a different contigs

	// compute mean and std on the fly
	float Mk = 0;
	float Qk = 0;
	uint32_t counterK = 1;
	//Keep header for further reference
	int32_t currentTid = -1;
	int32_t iSize;

	BamAlignment al;
	while ( bamFile.GetNextAlignmentCore(al) ) {
		reads ++;
		readStatus read_status = computeReadType(al, max_insert, is_mp);
		if (read_status != unmapped and read_status != lowQualty) {
			mappedReads ++;
			mappedReadsLength += al.Length;
		}

		if (al.IsFirstMate() && read_status == pair_proper) {
			iSize = abs(al.InsertSize);
			if(counterK == 1) {
				Mk = iSize;
				Qk = 0;
				counterK++;
			} else {
				float oldMk = Mk;
				float oldQk = Qk;
				Mk = oldMk + (iSize - oldMk)/counterK;
				Qk = oldQk + (counterK-1)*(iSize - oldMk)*(iSize - oldMk)/(float)counterK;
				counterK++;
			}
			insertsLength += iSize;
		}

		switch (read_status) {
		  case unmapped:
			  unmappedReads ++;
		     break;
		  case lowQualty:
			  lowQualityReads ++;
		     break;
		  case singleton:
			  singletonReads ++;
			  singletonReadsLength += al.Length ;
			  break;
		  case pair_wrongChrs:
			  matedDifferentContig ++;
			  matedDifferentContigLength += al.Length ;
			  break;
		  case pair_proper:
			  matedReads ++;
			  matedReadsLength += al.Length ;
			  break;
		  case pair_wrongDistance:
			  wrongDistanceReads ++;
			  wrongDistanceReadsLength += al.Length ;
			  break;
		  case pair_wrongOrientation:
			  wronglyOrientedReads ++;
			  wronglyOrientedReadsLength += al.Length ;
			  break;
		  default:
		     cout << "This should never be printed\n";
		     break;
		}

	}

	library.reads                 =  reads;
	library.mappedReads           =  mappedReads;
	library.unmappedReads         = unmappedReads;
	library.matedReads            = matedReads ;
	library.wrongDistanceReads    = wrongDistanceReads;
	library.lowQualityReads       = lowQualityReads ;
	library.wronglyOrientedReads  = wronglyOrientedReads ;
	library.matedDifferentContig  = matedDifferentContig ;
	library.singletonReads        =  singletonReads ;

	uint32_t total = matedReads + wrongDistanceReads +  wronglyOrientedReads +  matedDifferentContig + singletonReads  ;

	library.C_A = C_A = mappedReadsLength/(float)genomeLength;
	library.S_A = S_A = insertsLength/(float)genomeLength;
	library.C_M = C_M = matedReadsLength/(float)genomeLength;
	library.C_W = C_W = wronglyOrientedReadsLength/(float)genomeLength;
	library.C_S = C_S = singletonReadsLength/(float)genomeLength;
	library.C_D = C_D = matedDifferentContigLength/(float)genomeLength;
	library.insertMean = insertMean = Mk;
	Qk = sqrt(Qk/counterK);
	library.insertStd = insertStd = Qk;

	bamFile.Close();
	return library;
}


static void print_contigMetricsFileHeader(ofstream &ContigMetricsFile) {
	ContigMetricsFile << "contigID" << ","; //contigID
	ContigMetricsFile << "READ_COVERAGE" << ",";//read coverage
	ContigMetricsFile << "SPAN_COVERAGE" << ",";// span coverage
	ContigMetricsFile << "MEAN_INSERT_SIZE" << ",";//mean insert size
	ContigMetricsFile << "CORRECTLY_MATED_COV" << ",";//correctly mated coverage
	ContigMetricsFile << "WRONGLY_ORIENTED_COV" << ",";//wrongly oriented coverage
	ContigMetricsFile << "SINGLETON_COV" << ",";//singleton coverage
	ContigMetricsFile << "MATED_DIFFERENT_CTG_COV" << "\n";//Mated Different Contigs coverage

}

static void print_AssemblyMetrics(LibraryStatistics library, string type , ofstream &AssemblyMetricsFile) {
	AssemblyMetricsFile << "###LIBRARY STATISTICS\n";
	AssemblyMetricsFile << "BAM,LIB_TYPE,InsertSizeMean,InsertSizeStd,READS,MAPPED,UNMAPPED,PROPER,WRONG_DIST,ZERO_QUAL,WRONG_ORIENTATION,WRONG_CONTIG,";
	AssemblyMetricsFile << "SINGLETON,MEAN_COVERAGE,SPANNING_COVERAGE,PROPER_PAIRS_COVERAGE,WRONG_MATE_COVERAGE,SINGLETON_MATE_COV,DIFFERENT_CONTIG_COV\n";

	AssemblyMetricsFile <<  library.library_name <<  ",";
	AssemblyMetricsFile <<  type                 <<  ",";
	AssemblyMetricsFile <<  library.insertMean   <<  ",";
	AssemblyMetricsFile <<  library.insertStd    <<  ",";


	AssemblyMetricsFile << library.reads  				<<  ",";
	AssemblyMetricsFile << library.mappedReads 			<<  ",";
	AssemblyMetricsFile << library.unmappedReads		<<  ",";
	AssemblyMetricsFile << library.matedReads			<<  ",";
	AssemblyMetricsFile << 	library.wrongDistanceReads 	<<  ",";
	AssemblyMetricsFile << library.lowQualityReads 		<<  ",";
	AssemblyMetricsFile << library.wronglyOrientedReads <<  ",";
	AssemblyMetricsFile << library.matedDifferentContig	<<  ",";
	AssemblyMetricsFile << library.singletonReads		<<  ",";
	AssemblyMetricsFile << library.C_A 					<<  ",";
	AssemblyMetricsFile << library.S_A 					<<  ",";
	AssemblyMetricsFile << library.C_M 					<<  ",";
	AssemblyMetricsFile << library.C_W 					<<  ",";
	AssemblyMetricsFile << library.C_S 					<<  ",";
	AssemblyMetricsFile << library.C_D 					<<  "";

	AssemblyMetricsFile 								<<  "\n";


}







#endif /*TYPES_H_*/

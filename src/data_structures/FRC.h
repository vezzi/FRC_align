/*
 * FRC.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef FRC_H_
#define FRC_H_

using namespace std;

#include "options/Options.h"
#include "common.h"
#include "data_structures/Contig.h"
#include <iostream>
#include <algorithm>
#include <vector>




class contigFeatures {

	unsigned long int contigLength;
	unsigned int LOW_COVERAGE_AREA;
	unsigned int HIGH_COVERAGE_AREA;
	unsigned int LOW_NORMAL_AREA;
	unsigned int HIGH_NORMAL_AREA;
	unsigned int HIGH_SINGLE_AREA;
	unsigned int HIGH_SPANNING_AREA;
	unsigned int HIGH_OUTIE_AREA;
	unsigned int COMPRESSION_AREA;
	unsigned int STRECH_AREA;
	unsigned int TOTAL;

public:

	contigFeatures();
	~contigFeatures();

	unsigned long int getContigLength();
	void setContigLength(unsigned int contigLength);
	void setLOW_COVERAGE_AREA(unsigned int numFeat);
	void setHIGH_COVERAGE_AREA(unsigned int numFeat);
	void setLOW_NORMAL_AREA(unsigned int numFeat);
	void setHIGH_NORMAL_AREA(unsigned int numFeat);
	void setHIGH_SINGLE_AREA(unsigned int numFeat);
	void setHIGH_SPANNING_AREA(unsigned int numFeat);
	void setHIGH_OUTIE_AREA(unsigned int numFeat);
	void setCOMPRESSION_AREA(unsigned int numFeat);
	void setSTRECH_AREA(unsigned int numFeat);
	void setTOTAL(unsigned int numFeat);

	void computeTOTAL();

	void updateLOW_COVERAGE_AREA(unsigned int numFeat);
	void updateHIGH_COVERAGE_AREA(unsigned int numFeat);
	void updateLOW_NORMAL_AREA(unsigned int numFeat);
	void updateHIGH_NORMAL_AREA(unsigned int numFeat);
	void updateHIGH_SINGLE_AREA(unsigned int numFeat);
	void updateHIGH_SPANNING_AREA(unsigned int numFeat);
	void updateHIGH_OUTIE_AREA(unsigned int numFeat);
	void updateCOMPRESSION_AREA(unsigned int numFeat);
	void updateSTRECH_AREA(unsigned int numFeat);
	//void updateFeatures(windowStatistics* window)


	unsigned int getLOW_COVERAGE_AREA();
	unsigned int getHIGH_COVERAGE_AREA();
	unsigned int getLOW_NORMAL_AREA();
	unsigned int getHIGH_NORMAL_AREA();
	unsigned int getHIGH_SINGLE_AREA();
	unsigned int getHIGH_SPANNING_AREA();
	unsigned int getHIGH_OUTIE_AREA();
	unsigned int getCOMPRESSION_AREA();
	unsigned int getSTRECH_AREA();
	unsigned int getTOTAL();

	vector<pair<unsigned int, unsigned int> > LOW_COVERAGE_AREAS;
	vector<pair<unsigned int, unsigned int> > HIGH_COVERAGE_AREAS;
	vector<pair<unsigned int, unsigned int> > LOW_NORMAL_AREAS;
	vector<pair<unsigned int, unsigned int> > HIGH_NORMAL_AREAS;
	vector<pair<unsigned int, unsigned int> > HIGH_SINGLE_AREAS;
	vector<pair<unsigned int, unsigned int> > HIGH_SPANNING_AREAS;
	vector<pair<unsigned int, unsigned int> > HIGH_OUTIE_AREAS;
	vector<pair<unsigned int, unsigned int> > COMPRESSION_AREAS;
	vector<pair<unsigned int, unsigned int> > STRECH_AREAS;
	vector<pair<unsigned int, unsigned int> > TOTAL_AREAS;

	void print();

};



class FRC {

    vector<contigFeatures> CONTIG;
    unsigned int contigs;

    float C_A; // total read coverage
    float S_A; // total span coverage

    float C_M; // coverage induced by correctly aligned pairs
    float C_W; // coverage induced by wrongly mated pairs
    float C_S; // coverage induced by singletons
    float C_C; // coverage induced by reads with mate on a diferent contif

    float insertMean;
    float insertStd;

public:

	FRC();
	FRC(unsigned int contigs);
	~FRC();

	unsigned int getFeature(unsigned int ctg, Feature f);
	void setFeature(unsigned int ctg, Feature f, unsigned int value);
	void setContigLength(unsigned int ctg, unsigned int contigLength);
	unsigned int getContigLength(unsigned int ctg);

	void setC_A(float C_A);
	void setS_A(float S_A);
	void setC_M(float C_M);
	void setC_W(float C_W);
	void setC_S(float C_S);
	void setC_C(float C_C);
	void setInsertMean(float insertMean);
	void setInsertStd(float insertStd);

	void sortFRC();
	void computeLowCoverageArea(unsigned int ctg, Contig *contig);
	void computeHighCoverageArea(unsigned int ctg, Contig *contig);
	void computeLowNormalArea(unsigned int ctg, Contig *contig);
	void computeHighNormalArea(unsigned int ctg, Contig *contig);
	void computeHighSingleArea(unsigned int ctg, Contig *contig);
	void computeHighSpanningArea(unsigned int ctg, Contig *contig);
	void computeHighOutieArea(unsigned int ctg, Contig *contig);
	void computeCompressionArea(unsigned int ctg, Contig *contig);
	void computeStrechArea(unsigned int ctg, Contig *contig);

	void computeTOTAL(unsigned int ctg);

	void printContig(unsigned int ctg);


};




#endif /* FRC_H_ */

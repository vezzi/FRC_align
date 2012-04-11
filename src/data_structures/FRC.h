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

#include "Features.h"


class contigFeatures {
	string contigID;
	unsigned long int contigLength;

	unsigned int TOTAL;

public:


	Features PE;
	Features MP;

	contigFeatures();
	~contigFeatures();

	void setID(string ID);
	string getID();

	void setContigLength(unsigned int contigLength);
	unsigned long int getContigLength();


	unsigned int getTotal();

	vector<ternary> SUSPICIOUS_AREAS;

	void printFeatures(ofstream &file);

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
    float Expansion;
    float Compression;

    float insertMean;
    float insertStd;

public:

	FRC();
	FRC(unsigned int contigs);
	~FRC();

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
	void setID(unsigned int i, string ID);
	string  getID(unsigned int i);

	void sortFRC();
	void obtainCoverage(unsigned int ctg, Contig *contig);
	void computeLowCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int WindowSize, unsigned int WindowStep);
	void computeHighCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeLowNormalArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighNormalArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighSingleArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighSpanningArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighOutieArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeCompressionArea(string type, unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep);
	void computeStrechArea(string type, unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep);

	unsigned int getTotal(unsigned int ctg);


	void printFeatures(unsigned int ctg, ofstream &f);


};




#endif /* FRC_H_ */

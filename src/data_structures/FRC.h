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
	unsigned int HIGH_SPANING_AREA;
	unsigned int HIGH_OUTIE;
	unsigned int COMPRESSION_AREA;
	unsigned int STRECH_AREA;
	unsigned int TOTAL;

public:

	contigFeatures();
	~contigFeatures();

	unsigned long int getContigLength();
	void setContigLength(unsigned int contigLength);

	void updateLOW_COVERAGE_AREA();
	void updateHIGH_COVERAGE_AREA();
	void updateLOW_NORMAL_AREA();
	void updateHIGH_NORMAL_AREA();
	void updateHIGH_SINGLE_AREA();
	void updateHIGH_SPANING_AREA();
	void updateHIGH_OUTIE();
	void updateCOMPRESSION_AREA();
	void updateSTRECH_AREA();
	//void updateFeatures(windowStatistics* window)

	unsigned int getLOW_COVERAGE_AREA();
	unsigned int getHIGH_COVERAGE_AREA();
	unsigned int getLOW_NORMAL_AREA();
	unsigned int getHIGH_NORMAL_AREA();
	unsigned int getHIGH_SINGLE_AREA();
	unsigned int getHIGH_SPANING_AREA();
	unsigned int getHIGH_OUTIE();
	unsigned int getCOMPRESSION_AREA();
	unsigned int getSTRECH_AREA();
	unsigned int getTOTAL();

	void print();

};



class FRC {
    vector<contigFeatures> CONTIG;
    unsigned int contigs;

public:

	FRC();
	FRC(unsigned int contigs);
	~FRC();

	void update(unsigned int ctg, Feature f);
	unsigned int getFeature(unsigned int ctg, Feature f);
	void setContigLength(unsigned int ctg, unsigned int contigLength);
	unsigned int getContigLength(unsigned int ctg);


	void sortFRC();

};




#endif /* FRC_H_ */

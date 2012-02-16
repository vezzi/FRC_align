/*
 * Contig.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef CONTIG_H_
#define CONTIG_H_

using namespace std;


#include "data_structures/Features.h"
#include "samtools/sam.h"
#include "options/Options.h"
#include "common.h"


enum data {readCov, insertCov, cmCov, woCov, wdCov, singCov, mdcCov};

class Position {
public:
	unsigned int ReadCoverage;
	unsigned int StratingInserts;
	unsigned int InsertCoverage;
	unsigned int CorrectlyMated;
	unsigned int WronglyOriented;
	unsigned int WronglyDistance;
	unsigned int Singleton;
	unsigned int MatedDifferentContig;
	unsigned long int insertsLength;

	Position();

	~Position();

};


#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)

class Contig{
	unsigned int contigLength;
	unsigned int minInsert;
	unsigned int maxInsert;
	unsigned int windowSize;
	unsigned int windowStep;

	float lowCoverageFeat;
	float highCoverageFeat;
	float lowNormalFeat;
	float highNormalFeat;
	float highSingleFeat;
	float highSpanningFeat;
	float highOutieFeat;
	float CE_statistics;

	Position *CONTIG;

	void updateCov(unsigned int strat, unsigned int end, data type);

public:
	Contig();
	Contig(unsigned int contigLength, unsigned int minInsert, unsigned int maxInsert);
	~Contig();

	void updateContig(bam1_t* b); // given an alignment it updates the contig situation


	unsigned int getLowCoverageAreas(float C_A);
	unsigned int getHighCoverageAreas(float C_A);
	unsigned int getLowNormalAreas(float C_M);
	unsigned int getHighNormalAreas(float C_M);
	unsigned int getHighSingleAreas();
	unsigned int getHighSpanningAreas();
	unsigned int getHighOutieAreas();
	unsigned int getCompressionAreas(float insertionMean, float insertionStd);
	unsigned int getExpansionAreas(float insertionMean, float insertionStd);

	void print();

};





#endif /* CONTIG_H_ */

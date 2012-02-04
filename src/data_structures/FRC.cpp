/*
 * FRC.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#include "FRC.h"


bool sortContigs(contigFeatures i, contigFeatures j) {
	return (i.getContigLength() > j.getContigLength());
}

FRC::FRC() {

}

FRC::~FRC() {
	if (CONTIG.empty() != 0) {
		CONTIG.~vector();
	}
}

FRC::FRC(unsigned int contigs) {
	this->CONTIG.resize(contigs);
}


void FRC::sortFRC() {
	sort(CONTIG.begin(), CONTIG.end(), sortContigs);
}

void FRC::update(unsigned int ctg, Feature f) {
	switch(f) {
	case LOW_COVERAGE_AREA: this->CONTIG[ctg].updateLOW_COVERAGE_AREA(); break;
	case HIGH_COVERAGE_AREA: this->CONTIG[ctg].updateHIGH_COVERAGE_AREA() ; break;
	case LOW_NORMAL_AREA: this->CONTIG[ctg].updateLOW_NORMAL_AREA(); break;
	case HIGH_NORMAL_AREA: this->CONTIG[ctg].updateHIGH_NORMAL_AREA(); break;
	case HIGH_SINGLE_AREA: this->CONTIG[ctg].updateHIGH_SINGLE_AREA() ; break;
	case HIGH_SPANNING_AREA: this->CONTIG[ctg].updateHIGH_SPANNING_AREA() ; break;
	case HIGH_OUTIE_AREA: this->CONTIG[ctg].updateHIGH_OUTIE_AREA() ; break;
	case COMPRESSION_AREA: this->CONTIG[ctg].updateCOMPRESSION_AREA() ; break;
	case STRECH_AREA: this->CONTIG[ctg].updateSTRECH_AREA() ; break;
	default: cout << f << "undefined feature, exit with error\n";

	}
}

void FRC::setContigLength(unsigned int ctg, unsigned int contigLength) {
	CONTIG[ctg].setContigLength(contigLength);
}

unsigned int FRC::getContigLength(unsigned int ctg) {
	return this->CONTIG[ctg].getContigLength();
}

unsigned int FRC::getFeature(unsigned int ctg, Feature f) {
	switch(f) {
	case LOW_COVERAGE_AREA: return  this->CONTIG[ctg].getLOW_COVERAGE_AREA(); break;
	case HIGH_COVERAGE_AREA: return  this->CONTIG[ctg].getHIGH_COVERAGE_AREA() ; break;
	case LOW_NORMAL_AREA: return this->CONTIG[ctg].getLOW_NORMAL_AREA(); break;
	case HIGH_NORMAL_AREA:return  this->CONTIG[ctg].getHIGH_NORMAL_AREA(); break;
	case HIGH_SINGLE_AREA: return this->CONTIG[ctg].getHIGH_SINGLE_AREA() ; break;
	case HIGH_SPANNING_AREA: return this->CONTIG[ctg].getHIGH_SPANNING_AREA() ; break;
	case HIGH_OUTIE_AREA: return this->CONTIG[ctg].getHIGH_OUTIE_AREA() ; break;
	case COMPRESSION_AREA: return this->CONTIG[ctg].getCOMPRESSION_AREA() ; break;
	case STRECH_AREA: return this->CONTIG[ctg].getSTRECH_AREA() ; break;
	case TOTAL: return this->CONTIG[ctg].getTOTAL(); break;
	default: cout << f << "undefined feature, exit with error\n"; return 100;
	}
	return 100;
}




void FRC::setFeature(unsigned int ctg, Feature f, unsigned int value) {
	switch(f) {
	case LOW_COVERAGE_AREA:   this->CONTIG[ctg].setLOW_COVERAGE_AREA(value); break;
	case HIGH_COVERAGE_AREA: this->CONTIG[ctg].setHIGH_COVERAGE_AREA(value) ; break;
	case LOW_NORMAL_AREA: this->CONTIG[ctg].setLOW_NORMAL_AREA(value); break;
	case HIGH_NORMAL_AREA:  this->CONTIG[ctg].setHIGH_NORMAL_AREA(value); break;
	case HIGH_SINGLE_AREA: this->CONTIG[ctg].setHIGH_SINGLE_AREA(value) ; break;
	case HIGH_SPANNING_AREA: this->CONTIG[ctg].setHIGH_SPANNING_AREA(value) ; break;
	case HIGH_OUTIE_AREA:  this->CONTIG[ctg].setHIGH_OUTIE_AREA(value) ; break;
	case COMPRESSION_AREA: this->CONTIG[ctg].setCOMPRESSION_AREA(value) ; break;
	case STRECH_AREA:  this->CONTIG[ctg].setSTRECH_AREA(value) ; break;
	case TOTAL: this->CONTIG[ctg].setTOTAL(value); break;
	default: cout << f << "undefined feature, exit with error\n";;
	}
}

void FRC::computeLowCoverageArea(unsigned int ctg, Contig *contig) {
	unsigned int lowCoverageFeatures = contig->getLowCoverageAreas(this->C_A);
	this->CONTIG[ctg].setLOW_COVERAGE_AREA(lowCoverageFeatures);
}

void FRC::computeHighCoverageArea(unsigned int ctg, Contig *contig) {
	unsigned int highCoverageFeatures = contig->getHighCoverageAreas(this->C_A);
	this->CONTIG[ctg].setHIGH_COVERAGE_AREA(highCoverageFeatures);
}

void FRC::computeLowNormalArea(unsigned int ctg, Contig *contig) {
	unsigned int lowNormalFeatures = contig->getLowNormalAreas(this->C_M);
	this->CONTIG[ctg].setLOW_NORMAL_AREA(lowNormalFeatures);
}

void FRC::computeHighNormalArea(unsigned int ctg, Contig *contig) {
	unsigned int highNormalFeatures = contig->getHighNormalAreas(this->C_M);
	this->CONTIG[ctg].setHIGH_NORMAL_AREA(highNormalFeatures);

}

void FRC::computeHighSingleArea(unsigned int ctg, Contig *contig) {
	unsigned int highSingleFeatures = contig->getHighSingleAreas();
	this->CONTIG[ctg].setHIGH_SINGLE_AREA(highSingleFeatures);
}

void FRC::computeHighSpanningArea(unsigned int ctg, Contig *contig) {
	unsigned int highSpanningFeatures = contig->getHighSpanningAreas();
	this->CONTIG[ctg].setHIGH_SPANNING_AREA(highSpanningFeatures);
}

void FRC::computeHighOutieArea(unsigned int ctg, Contig *contig) {
	unsigned int highOutieFeatures = contig->getHighOutieAreas();
	this->CONTIG[ctg].setHIGH_OUTIE_AREA(highOutieFeatures);

}

void FRC::computeCompressionArea(unsigned int ctg, Contig *contig) {
	unsigned int compressionFeatures = contig->getCompressionAreas(this->insertMean, this->insertStd);
	this->CONTIG[ctg].setCOMPRESSION_AREA(compressionFeatures);

}

void FRC::computeStrechArea(unsigned int ctg, Contig *contig) {
	unsigned int expansionFeatures = contig->getExpansionAreas(this->insertMean, this->insertStd);
	this->CONTIG[ctg].setSTRECH_AREA(expansionFeatures);
}

void FRC::computeTOTAL(unsigned int ctg) {
	this->CONTIG[ctg].computeTOTAL();
}



void FRC::setC_A(float C_A) {
	this->C_A = C_A;
}
void FRC::setS_A(float S_A) {
	this->S_A = S_A;
}
void FRC::setC_M(float C_M) {
	this->C_M = C_M;
}
void FRC::setC_W(float C_W) {
	this->C_W = C_W;
}
void FRC::setC_S(float C_S) {
	this->C_S = C_S;
}
void FRC::setC_C(float C_C) {
	this->C_C = C_C;
}
void FRC::setInsertMean(float insertMean) {
	this->insertMean = insertMean;
}
void FRC::setInsertStd(float insertStd) {
	this->insertStd = insertStd;
}



void FRC::printContig(unsigned int ctg) {
	CONTIG[ctg].print();
}



contigFeatures::contigFeatures() {
	contigLength = 0;
	LOW_COVERAGE_AREA = 0;
	HIGH_COVERAGE_AREA = 0;
	LOW_NORMAL_AREA = 0;
	HIGH_NORMAL_AREA = 0;
	HIGH_SINGLE_AREA = 0;
	HIGH_SPANNING_AREA = 0;
	HIGH_OUTIE_AREA = 0;
	COMPRESSION_AREA = 0;
	STRECH_AREA = 0;
	TOTAL = 0;
}

contigFeatures::~contigFeatures() {

}


unsigned long int contigFeatures::getContigLength() {
	return this->contigLength;
}

void contigFeatures::setContigLength(unsigned int contigLength) {
	this->contigLength = contigLength;
}


void contigFeatures::setLOW_COVERAGE_AREA(unsigned int numFeat) {
	this->LOW_COVERAGE_AREA = numFeat;
}
void contigFeatures::setHIGH_COVERAGE_AREA(unsigned int numFeat) {
	this->HIGH_COVERAGE_AREA = numFeat;
}
void contigFeatures::setLOW_NORMAL_AREA(unsigned int numFeat) {
	this->LOW_NORMAL_AREA = numFeat;
}
void contigFeatures::setHIGH_NORMAL_AREA(unsigned int numFeat) {
	this->HIGH_NORMAL_AREA = numFeat;
}
void contigFeatures::setHIGH_SINGLE_AREA(unsigned int numFeat) {
	this->HIGH_SINGLE_AREA = numFeat;
}
void contigFeatures::setHIGH_SPANNING_AREA(unsigned int numFeat) {
	this->HIGH_SPANNING_AREA = numFeat;
}
void contigFeatures::setHIGH_OUTIE_AREA(unsigned int numFeat) {
	this->HIGH_OUTIE_AREA = numFeat;
}
void contigFeatures::setCOMPRESSION_AREA(unsigned int numFeat) {
	this->COMPRESSION_AREA = numFeat;
}
void contigFeatures::setSTRECH_AREA(unsigned int numFeat) {
	this->STRECH_AREA = numFeat;
}

void contigFeatures::computeTOTAL() {
	this->TOTAL = this->COMPRESSION_AREA + this->HIGH_COVERAGE_AREA + this->HIGH_NORMAL_AREA + this->HIGH_OUTIE_AREA +
			this->HIGH_SINGLE_AREA + this->HIGH_SINGLE_AREA + this->HIGH_SPANNING_AREA + this->LOW_COVERAGE_AREA +
			this->LOW_NORMAL_AREA + this->STRECH_AREA ;
}


void contigFeatures::setTOTAL(unsigned int numFeat) {
	this->STRECH_AREA = numFeat;
}


void contigFeatures::updateLOW_COVERAGE_AREA() {
	LOW_COVERAGE_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_COVERAGE_AREA() {
	HIGH_COVERAGE_AREA++;
	TOTAL++;
}

void contigFeatures::updateLOW_NORMAL_AREA() {
	LOW_NORMAL_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_NORMAL_AREA() {
	HIGH_NORMAL_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_SINGLE_AREA() {
	HIGH_SINGLE_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_SPANNING_AREA() {
	HIGH_SPANNING_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_OUTIE_AREA() {
	HIGH_OUTIE_AREA++;
	TOTAL++;
}

void contigFeatures::updateCOMPRESSION_AREA() {
	COMPRESSION_AREA++;
	TOTAL++;
}

void contigFeatures::updateSTRECH_AREA() {
	STRECH_AREA++;
	TOTAL++;
}



/*	void updateFeatures(windowStatistics* window) {
		float C_A_i = window->readsLength_win/(float)window->insertsLength_win;
		float S_A_i = window->insertsLength_win/(float)window->insertsLength_win;
		float C_M_i = window->correctlyMatedReadsLength_win/(float)window->insertsLength_win;
		float C_W_i = (window->wronglyDistanceReadsLength_win + window->wronglyOrientedReadsLength_win)/(float)window->insertsLength_win;
		float C_S_i = window->singletonReadsLength_win/(float)window->insertsLength_win;
		float C_C_i = window->matedDifferentContigLength_win/(float)window->insertsLength_win;
	}
 */


unsigned int contigFeatures::getLOW_COVERAGE_AREA() {return LOW_COVERAGE_AREA;}
unsigned int contigFeatures::getHIGH_COVERAGE_AREA() {return HIGH_COVERAGE_AREA;}
unsigned int contigFeatures::getLOW_NORMAL_AREA() {return LOW_NORMAL_AREA;}
unsigned int contigFeatures::getHIGH_NORMAL_AREA() {return HIGH_NORMAL_AREA;}
unsigned int contigFeatures::getHIGH_SINGLE_AREA() {return HIGH_SINGLE_AREA;}
unsigned int contigFeatures::getHIGH_SPANNING_AREA() {return HIGH_SPANNING_AREA;}
unsigned int contigFeatures::getHIGH_OUTIE_AREA() {return HIGH_OUTIE_AREA;}
unsigned int contigFeatures::getCOMPRESSION_AREA() {return COMPRESSION_AREA;}
unsigned int contigFeatures::getSTRECH_AREA() {return STRECH_AREA;}
unsigned int contigFeatures::getTOTAL() {return TOTAL;}




void contigFeatures::print() {
	cout << "contigLength " << contigLength << "\n";
	cout << "LOW_COVERAGE_AREA " << LOW_COVERAGE_AREA << "\n";
	cout << "HIGH_COVERAGE_AREA " << HIGH_COVERAGE_AREA << "\n";
	cout << "LOW_NORMAL_AREA " << LOW_NORMAL_AREA << "\n";
	cout << "HIGH_NORMAL_AREA " << HIGH_NORMAL_AREA << "\n";
	cout << "HIGH_SINGLE_AREA " << HIGH_SINGLE_AREA << "\n";
	cout << "HIGH_SPANNING_AREA " << HIGH_SPANNING_AREA << "\n";
	cout << "HIGH_OUTIE_AREA " << HIGH_OUTIE_AREA << "\n";
	cout << "COMPRESSION_AREA " << COMPRESSION_AREA << "\n";
	cout << "STRECH_AREA " << STRECH_AREA << "\n";
	cout << "TOTAL " << TOTAL << "\n-----\n";
}





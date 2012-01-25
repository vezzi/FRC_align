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
	CONTIG.~vector();
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
	case HIGH_SPANING_AREA: this->CONTIG[ctg].updateHIGH_SPANING_AREA() ; break;
	case HIGH_OUTIE: this->CONTIG[ctg].updateHIGH_OUTIE() ; break;
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
	case HIGH_SPANING_AREA: return this->CONTIG[ctg].getHIGH_SPANING_AREA() ; break;
	case HIGH_OUTIE: return this->CONTIG[ctg].getHIGH_OUTIE() ; break;
	case COMPRESSION_AREA: return this->CONTIG[ctg].getCOMPRESSION_AREA() ; break;
	case STRECH_AREA: return this->CONTIG[ctg].getSTRECH_AREA() ; break;
	case TOTAL: return this->CONTIG[ctg].getTOTAL(); break;
	default: cout << f << "undefined feature, exit with error\n"; return 100;
	}
	return 100;
}




contigFeatures::contigFeatures() {
	contigLength = 0;
	LOW_COVERAGE_AREA = 0;
	HIGH_COVERAGE_AREA = 0;
	LOW_NORMAL_AREA = 0;
	HIGH_NORMAL_AREA = 0;
	HIGH_SINGLE_AREA = 0;
	HIGH_SPANING_AREA = 0;
	HIGH_OUTIE = 0;
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

void contigFeatures::updateHIGH_SPANING_AREA() {
	HIGH_SPANING_AREA++;
	TOTAL++;
}

void contigFeatures::updateHIGH_OUTIE() {
	HIGH_OUTIE++;
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
unsigned int contigFeatures::getHIGH_SPANING_AREA() {return HIGH_SPANING_AREA;}
unsigned int contigFeatures::getHIGH_OUTIE() {return HIGH_OUTIE;}
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
	cout << "HIGH_SPANING_AREA " << HIGH_SPANING_AREA << "\n";
	cout << "HIGH_OUTIE " << HIGH_OUTIE << "\n";
	cout << "COMPRESSION_AREA " << COMPRESSION_AREA << "\n";
	cout << "STRECH_AREA " << STRECH_AREA << "\n";
	cout << "TOTAL " << TOTAL << "\n-----\n";
}





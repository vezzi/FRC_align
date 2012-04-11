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
	this->contigs = contigs;
	this->CONTIG.resize(contigs);
}


void FRC::sortFRC() {
	sort(CONTIG.begin(), CONTIG.end(), sortContigs);
}


void FRC::setContigLength(unsigned int ctg, unsigned int contigLength) {
	CONTIG[ctg].setContigLength(contigLength);
}

unsigned int FRC::getContigLength(unsigned int ctg) {
	return this->CONTIG[ctg].getContigLength();
}


void FRC::obtainCoverage(unsigned int ctg, Contig *contig) {


	float coverage = contig->getCoverage();
	cout << "contig " << ctg << " has mean cov " << coverage << "\n";

}



void FRC::computeLowCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getLowCoverageAreas(C_A,windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateLOW_COVERAGE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateLOW_COVERAGE_AREA(feat);
	}
	for(unsigned int i=0; i< contig->lowCoverageAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "LOW_COV_"+type;
		tmp.start = contig->lowCoverageAreas.at(i).first;
		tmp.end = contig->lowCoverageAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeHighCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighCoverageAreas(this->C_A, windowSize, windowStep);
	cout << "contig " << ctg << " has " << feat << " HIGH COVERAGE FEATS\n";
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_COVERAGE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_COVERAGE_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highCoverageAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_COV_"+type;
		tmp.start = contig->highCoverageAreas.at(i).first;
		tmp.end = contig->highCoverageAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}

void FRC::computeLowNormalArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getLowNormalAreas(this->C_M, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateLOW_NORMAL_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateLOW_NORMAL_AREA(feat);
	}
	for(unsigned int i=0; i < contig->lowNormalAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "LOW_NORM_COV_"+type;
		tmp.start = contig->lowNormalAreas.at(i).first;
		tmp.end = contig->lowNormalAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}

void FRC::computeHighNormalArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighNormalAreas(this->C_M, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_NORMAL_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_NORMAL_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highNormalAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_NORM_COV_"+type;
		tmp.start = contig->highNormalAreas.at(i).first;
		tmp.end = contig->highNormalAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}

void FRC::computeHighSingleArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighSingleAreas( windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_SINGLE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_SINGLE_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highSingleAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_SINGLE_"+type;
		tmp.start = contig->highSingleAreas.at(i).first;
		tmp.end = contig->highSingleAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);

	}

}

void FRC::computeHighSpanningArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighSpanningAreas( windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_SPANNING_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_SPANNING_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highSpanningAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_SPAN_"+type;
		tmp.start = contig->highSpanningAreas.at(i).first;
		tmp.end = contig->highSpanningAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);

	}
}

void FRC::computeHighOutieArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighOutieAreas( windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_OUTIE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_OUTIE_AREA(feat);
	}

	for(unsigned int i=0; i < contig->highOutieAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_OUTIE_"+type;
		tmp.start = contig->highOutieAreas.at(i).first;
		tmp.end = contig->highOutieAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeCompressionArea(string type,unsigned int ctg, Contig *contig, float Zscore , unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getCompressionAreas(this->insertMean, this->insertStd, Zscore, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateCOMPRESSION_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateCOMPRESSION_AREA(feat);
	}

	for(unsigned int i=0; i< contig->compressionAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "COMPR_"+type;
		tmp.start = contig->compressionAreas.at(i).first;
		tmp.end = contig->compressionAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeStrechArea(string type,unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getExpansionAreas(this->insertMean, this->insertStd, Zscore, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateSTRECH_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateSTRECH_AREA(feat);
	}

	for(unsigned int i=0; i < contig->expansionAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "STRECH_"+type;
		tmp.start = contig->expansionAreas.at(i).first;
		tmp.end = contig->expansionAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
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

void FRC::setID(unsigned int i, string ID) {
	CONTIG[i].setID(ID);
}

string  FRC::getID(unsigned int i) {
	return CONTIG[i].getID();
}





void FRC::printFeatures(unsigned int ctg, ofstream &file) {
	this->CONTIG[ctg].printFeatures(file);

}


unsigned int FRC::getTotal(unsigned int ctg) {
	return this->CONTIG[ctg].getTotal();
}



contigFeatures::contigFeatures() {
	contigLength = 0;
	TOTAL = 0;
	SUSPICIOUS_AREAS.clear();
}

contigFeatures::~contigFeatures() {

}

void contigFeatures::setContigLength(unsigned int contigLength) {
	this->contigLength = contigLength;
}

unsigned long int contigFeatures::getContigLength() {
	return this->contigLength;
}




void contigFeatures::setID(string ID) {
	this->contigID = ID;
}

string contigFeatures::getID() {
	return this->contigID;
}

unsigned int contigFeatures::getTotal() {
	TOTAL = PE.returnTotal() + MP.returnTotal();
	return TOTAL;
}

bool sortTernary(ternary t1, ternary t2) {return (t1.start < t2.start);}

void contigFeatures::printFeatures(ofstream &file) {
	unsigned int total = 0;
	sort(SUSPICIOUS_AREAS.begin(), SUSPICIOUS_AREAS.end(), sortTernary);
	for(unsigned int i=0; i < SUSPICIOUS_AREAS.size(); i++) {
		file << this->contigID << " " << SUSPICIOUS_AREAS[i].feature << " " << SUSPICIOUS_AREAS[i].start << " " << SUSPICIOUS_AREAS[i].end << "\n";
		total += floor((SUSPICIOUS_AREAS[i].end - SUSPICIOUS_AREAS[i].start + 1)/(float)1000 + 0.5);
	}
	//this->TOTAL = total;


}







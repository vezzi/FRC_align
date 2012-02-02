/*
 * Contig.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#include "Contig.h"

Position::Position() {
	ReadCoverage = 0;
	StratingInserts = 0;
	InsertCoverage = 0;
	CorrectlyMated = 0;
	WronglyOriented = 0;
	WronglyDistance = 0;
	Singleton = 0;
	MatedDifferentContig = 0;
}

Position::~Position() {

}


Contig::Contig() {
	contigLength = 0;
	peMinInsert = 0;
	peMaxInsert = 0;
	windowSize = 1000;
	windowStep = 200;

	lowCoverageFeat = 1/(float)3;
	highCoverageFeat = 3;
	lowNormalFeat = 1/(float)3;
	highNormalFeat = 3;
	highSingleFeat = 0.6;
	highSpanningFeat = 0.6;
	highOutieFeat = 0.6;
	CE_statistics = 3;
}

Contig::Contig(unsigned int contigLength, unsigned int peMinInsert, unsigned int peMaxInsert) {
	this->contigLength = contigLength;
	this->peMinInsert = peMinInsert;
	this->peMaxInsert = peMaxInsert;
	this->CONTIG =  new Position[contigLength];
	windowSize = 1000;
	windowStep = 200;

	lowCoverageFeat = 1/(float)3;
	highCoverageFeat = 3;
	lowNormalFeat = 1/(float)3;
	highNormalFeat = 3;
	highSingleFeat = 0.6;
	highSpanningFeat = 0.6;
	highOutieFeat = 0.6;
	CE_statistics = 3;
}

Contig::~Contig() {
	if(CONTIG != NULL) {
		delete [] CONTIG;
	}
}


void Contig::updateCov(unsigned int start, unsigned int end, data type) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
//		cout << "hoops, end longer than contig length when updating CONTIG " << type << "\n";
//		cout << "\tcontig length " << this->contigLength << " starting point " << start << " ending point " << end << "\n";
		end = this->contigLength;
	}
	// now update
	if(type == insertCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].InsertCoverage++;
	} else if(type == readCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].ReadCoverage++;
	} else if(type ==  cmCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].CorrectlyMated++;
	} else if(type == woCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].WronglyOriented++;
	}else if(type == wdCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].WronglyDistance++;
	} else if(type == singCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].Singleton++;
	} else if(type == mdcCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].MatedDifferentContig++;
	} else {
		cout << "hoops, unknown type " << type << " there must be something wrong!!!\n";
	}
}






void Contig::updateContig(bam1_t* b) {
	const bam1_core_t* core =  &b->core;
	uint32_t* cigar = bam1_cigar(b);
	int32_t alignmentLength = 0;
	int32_t startRead=0;
	int32_t endRead=0;
	int32_t startPaired=0;
	int32_t startInsert=0;
	int32_t endInsert=0;
	uint32_t iSize=0;
	if(!(core->flag&BAM_FUNMAP) && !(core->flag&BAM_FDUP) && !(core->flag&BAM_FSECONDARY) && !(core->flag&BAM_FQCFAIL)) { // if read has been mapped and it is not a DUPLICATE or a SECONDARY alignment
		alignmentLength = bam_cigar2qlen(core,cigar);
		startRead = core->pos; // start position on the contig
		endRead = startRead + alignmentLength ; // position where reads ends
		updateCov(startRead, endRead, readCov); // update coverage

		iSize = abs(core->isize);

		if ((core->flag&BAM_FREAD1) //First in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			startPaired = core->mpos;
			if(startRead < startPaired) {
				iSize = (startPaired + core->l_qseq -1) - startRead; // insert size, I consider both reads of the same length
				startInsert = startRead;
				endInsert = startRead + iSize;
				if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					//cout << "correctly oriented\n";
					if (peMinInsert <= iSize && iSize <= peMaxInsert) { //this is a right insert
						updateCov(startRead, endRead, cmCov); // update good read coverage
						CONTIG[startRead].StratingInserts++; // another insert starts here
						updateCov(startInsert,endInsert, insertCov); // update spanning coverage
					} else {
						updateCov(startRead, endRead, wdCov);
					}
				} else {
					//pair is wrongly oriented
					//cout << "wrongly aligned\n";
					updateCov(startRead, endRead, woCov);
				}
			} else {
				iSize = (startRead + alignmentLength - 1) - startPaired;
				startInsert = startPaired;
				endInsert = startInsert + iSize;
				if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					//cout << "correctly oriented\n";
					if (peMinInsert <= iSize && iSize <= peMaxInsert) { //this is a right insert
						updateCov(startRead, endRead, cmCov); // update good read coverage
						CONTIG[startRead].StratingInserts++; // another insert starts here
						updateCov(startInsert,endInsert, insertCov); // update spanning coverage
					} else {
						updateCov(startRead, endRead, wdCov);
					}
				} else {
					updateCov(startRead, endRead, wdCov);
				}
			}
		} else  if ((core->flag&BAM_FREAD2) //Second in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			startPaired = core->mpos;
			if(startRead > startPaired) {
				iSize = (startRead + alignmentLength -1) - startPaired;
				if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					//cout << "correctly oriented\n";
					if (peMinInsert <= iSize && iSize <= peMaxInsert) { //this is a right insert, no need to update insert coverage
						updateCov(startRead, endRead, cmCov); // update good read coverage
					} else {
						updateCov(startRead, endRead, wdCov);
					}
				} else {
					//pair is wrongly oriented
					//cout << "wrongly aligned\n";
					updateCov(startRead, endRead, woCov);
				}
			} else {
				iSize = (startPaired + core->l_qseq -1) - startRead;
				if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					//cout << "correctly oriented\n";
					if (peMinInsert <= iSize && iSize <= peMaxInsert) { //this is a right insert, no need to update insert coverage
						updateCov(startRead, endRead, cmCov); // update good read coverage
					} else {
						updateCov(startRead, endRead, wdCov);
					}
				} else {
					//pair is wrongly oriented
					//cout << "wrongly aligned\n";
					updateCov(startRead, endRead, woCov);
				}
			}
		} else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
			//Count inter-chrom pairs
			updateCov(startRead, endRead, mdcCov);
//			actualWindow->matedDifferentContigLength_win += bam_cigar2qlen(core,cigar);
		} else if(core->flag&BAM_FMUNMAP) {
			// if mate read is unmapped
			updateCov(startRead, endRead, singCov);
//			actualWindow->singletonReadsLength_win =+ bam_cigar2qlen(core,cigar);
		}

	}
	//cout << "readStart " << startRead << " readEnd " << endRead << " iSize " << iSize << " insStart " << startInsert << " insEnd " << endInsert << "\n";

}



void Contig::print() {
	cout << "Contig size " << this->contigLength << "\n";
	for(unsigned int i= 0; i < this->contigLength; i++) {
		if(i % 6 == 0 && i > 0) {
			cout << "\n";
		}
		cout << "(" << i << ":" << this->CONTIG[i].ReadCoverage  << "," << this->CONTIG[i].InsertCoverage << "," << this->CONTIG[i].CorrectlyMated << "," <<
				this->CONTIG[i].MatedDifferentContig  << "," << "," <<   this->CONTIG[i].Singleton << "," <<   this->CONTIG[i].WronglyOriented
				<< "," << this->CONTIG[i].StratingInserts << ") " ;
	}
	cout << "\n\n";
}


unsigned int Contig::getLowCoverageAreas(float C_A) {
	//compute length of low coverage areas in the contig
	//use a 1K sliding window
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < this->windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov < lowCoverageFeat*C_A ) { // this is a feature
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov < lowCoverageFeat*C_A ) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}

		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov < lowCoverageFeat*C_A ) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
//TODO compute statistics for the eventual last overlapping window
		if(feat) { // a feature are reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
		}
	}


	return features;
}



unsigned int Contig::getHighCoverageAreas(float C_A) {
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < this->windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov > highCoverageFeat*C_A ) { // this is a feature
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov > highCoverageFeat*C_A ) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}

		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov > highCoverageFeat*C_A ) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
//TODO compute statistics for the eventual last overlapping window
		if(feat) { // a feature are reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
		}
	}


	return features;
}




unsigned int Contig::getLowNormalAreas(float C_M) {
	//compute length of low coverage areas in the contig
	//use a 1K sliding window
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < this->windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].CorrectlyMated ;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov < lowNormalFeat*C_M ) { // this is a feature
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov < lowNormalFeat*C_M ) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov < lowNormalFeat*C_M ) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
//TODO compute statistics for the eventual last overlapping window
		if(feat) { // a feature are reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
		}
	}
	return features;
}



unsigned int Contig::getHighNormalAreas(float C_M) {
	//compute length of low coverage areas in the contig
	//use a 1K sliding window
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < this->windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].CorrectlyMated ;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov > highNormalFeat*C_M ) { // this is a feature
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov > highNormalFeat*C_M ) {  // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov > highNormalFeat*C_M ) {  // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
//TODO compute statistics for the eventual last overlapping window
		if(feat) { // a feature are reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
		}
	}
	return features;
}

unsigned int Contig::getHighSingleAreas( ) {
	unsigned int totalCoverage = 0;
	unsigned int singleReadCoverage = 0;
	unsigned int features = 0;
	float meanTotalCov;
	float meanSingleCov;
	if(this->contigLength < this->windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			singleReadCoverage += CONTIG[i].Singleton;
		}
		meanTotalCov = totalCoverage/(float)this->contigLength; // this is the "window" total coverage
		meanSingleCov = singleReadCoverage/(float)this->contigLength; // this is the "window" single read coverage
		if( meanSingleCov > highSingleFeat*meanTotalCov ) { // this is a feature
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			singleReadCoverage += CONTIG[i].Singleton;
		}
		meanTotalCov = totalCoverage/(float)this->contigLength; // this is the "window" total coverage
		meanSingleCov = singleReadCoverage/(float)this->contigLength; // this is the "window" single read coverage
		if( meanSingleCov > highSingleFeat*meanTotalCov ) { //first window's covrage
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			singleReadCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage ;
				singleReadCoverage += CONTIG[i].Singleton;
			}
			meanTotalCov = totalCoverage/(float)(endWindow - startWindow); // compute window total coverage
			meanSingleCov = singleReadCoverage/(float)(endWindow - startWindow); // compute window single read coverage
			if( meanSingleCov > highSingleFeat*meanTotalCov ) {
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
//TODO compute statistics for the eventual last overlapping window
		if(feat) { // a feature are reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
		}
	}
	return features;
}











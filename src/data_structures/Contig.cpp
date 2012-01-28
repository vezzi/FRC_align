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


Contig::Contig() {
	contigLength = 0;
	peMinInsert = 0;
	peMaxInsert = 0;
}

Contig::Contig(unsigned int contigLength, unsigned int peMinInsert, unsigned int peMaxInsert) {
	this->contigLength = contigLength;
	this->peMinInsert = peMinInsert;
	this->peMaxInsert = peMaxInsert;
	CONTIG = new Position[contigLength];
}

Contig::~Contig() {
	delete CONTIG;
}


void Contig::updateCov(unsigned int start, unsigned int end, data type) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
		cout << "hoops, end longer than contig length when updating CONTIG\n";
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
	int32_t start;
	int32_t end;
	uint32_t iSize;
	if (core->qual > 0) {
		uint32_t* cigar = bam1_cigar(b);

		if ((core->flag&BAM_FREAD1) //First in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			//pair is mapped on the same contig and I'm working on the first pair
			start = MIN(core->pos,core->mpos);
			end = start+abs(core->isize);
			iSize = end-start; // compute insert size

			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				CONTIG[start].StratingInserts++; // another insert starts here
				updateCov(start,end, insertCov);
			}

			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				if(core->pos < core->mpos) { // read I'm processing is the first
					if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, cmCov);
					} else {
						// wrong orientation
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, woCov);
					}
				} else {
					if(!(core->flag&BAM_FMREVERSE) && (core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, cmCov);   //TODO: chack this I`m not really sure about ranges!!!
//						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
					} else {
						// wrong orientation
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, woCov);
//						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				}
			} else {
				//wrong distance
				start = core->pos; // start position on the contig
				end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
				updateCov(start, end, wdCov);
//				actualWindow->wronglyDistanceReadsLength_win += bam_cigar2qlen(core,cigar);
			}
		} else  if ((core->flag&BAM_FREAD2) //Second in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		)
			// if I'm considering the second read in a pair I must check it is is a correctly mated read and if this is tha case update the right variables
		{
			start = MIN(core->pos,core->mpos);
			end = start+abs(core->isize);
			iSize = end-start; // compute insert size
// NOW I CAN SKIP THIS, BECAUSE I'M WORKING ON THE ALL CONTIG
//			if (peMinInsert <= iSize && iSize <= peMaxInsert) { // I have to check if the mate is outside window boundaries
//				if(start <= actualWindow->windowStart || end >= actualWindow->windowEnd) {
//					actualWindow->insertsLength_win += iSize; // update number of inserts and total length
//					actualWindow->inserts++;
//				}
//			}
			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				if(core->pos > core->mpos) { // read I'm processing is the first
					if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, cmCov);
//						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); //  update number of correctly mapped and their length
					} else {
						// wrong orientation
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, woCov);
//						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				} else {
					if((core->flag&BAM_FMREVERSE) && !(core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, cmCov);
//						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
					} else {
						// wrong orientation
						start = core->pos; // start position on the contig
						end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
						updateCov(start, end, woCov);
//						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				}
			} else {
				//wrong distance
				start = core->pos; // start position on the contig
				end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
				updateCov(start, end, wdCov);
//				actualWindow->wronglyDistanceReadsLength_win += bam_cigar2qlen(core,cigar);

			}
		} else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
			//Count inter-chrom pairs
			start = core->pos; // start position on the contig
			end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
			updateCov(start, end, mdcCov);
//			actualWindow->matedDifferentContigLength_win += bam_cigar2qlen(core,cigar);
		} else if(core->flag&BAM_FMUNMAP) {
			// if mate read is unmapped
			start = core->pos; // start position on the contig
			end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
			updateCov(start, end, singCov);
//			actualWindow->singletonReadsLength_win =+ bam_cigar2qlen(core,cigar);
		}

		if (core->flag&BAM_FDUP) {   //This is a duplicate. Don't count it!.
		} else { // otherwise always increase the coverage
			start = core->pos; // start position on the contig
			end = start + bam_cigar2qlen(core,cigar) - 1; // position where reads ends
			updateCov(start, end, readCov);
//			actualWindow->readsLength_win += bam_cigar2qlen(core,cigar);
		}
	}
}



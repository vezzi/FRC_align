/*
 * Translocations.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi
 */

#include "Translocation.h"



bool sortPairs(pair<uint32_t, uint32_t> i, pair<uint32_t, uint32_t>  j) {
	return (i.first < j.first);
}

bool sortLinksChr2(Link i, Link  j) {
	return (i.chr2_start < j.chr2_start);
}

bool sortLinksChr1(Link i, Link  j) {
	return (i.chr1_start < j.chr1_start);
}






Window::Window(int windowSize, int windowStep, int max_insert, uint16_t minimum_mapping_quality,
		bool outtie, float mean_insert, float std_insert, int minimumPairs,
		float meanCoverage, string outputFileHeader) {
	this->windowSize 		 = windowSize;
	this->windowStep 		 = windowStep;
	this->max_insert		 = max_insert;
	this->minimum_mapping_quality = minimum_mapping_quality;
	this->outtie			 = outtie;
	this->mean_insert		 = mean_insert;
	this ->std_insert		 = std_insert;
	this->minimumPairs		 = minimumPairs;
	this->meanCoverage		 = meanCoverage;

	this->outputFileHeader   = outputFileHeader;
	string inter_chr_events = outputFileHeader + "_inter_chr_events.tab";
	this->interChrVariations.open(inter_chr_events.c_str());
	this->interChrVariations << "chrA\tstartOnA\tendOnA\tchrB\tstartOnB\tendOnB\tLinksFromWindow\tLinksToChrB\tLinksInCurrentEvent\tCoverageOnChrA\tExpectedLinks<<RatioEL_LiCE\n";
	string intra_chr_events = outputFileHeader + "_intra_chr_events.tab";
	this->intraChrVariations.open(intra_chr_events.c_str());
	this->intraChrVariations << "chrA\tstartOnA\tendOnA\tchrB\tstartOnB\tendOnB\tLinksFromWindow\tLinksToChrB\tLinksInCurrentEvent\tCoverageOnChrA\tExpectedLinks<<RatioEL_LiCE\n";

	this->coverage          	= 0;
	this->currentWindowStart	= 0;
	this->currentWindowEnd		= 0;
	this->windowOpen			= false;
	this->chr					=-1;

}




void Window::initTrans(SamHeader head) {
	uint32_t contigsNumber = 0;
	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		this->contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
		this->position2contig[contigsNumber]  = sequence->Name;
		contigsNumber++;
	}

}

void Window::insertRead(BamAlignment alignment) {

	readStatus alignmentStatus = computeReadType(alignment, this->max_insert, this->outtie);
	//if(alignmentStatus == unmapped or alignmentStatus == lowQualty ) {
	//	return; // in case the alignment is of no use discard it
	//}
	//TODO: maybe do not count singletons

	if(this->chr == -1) { // first read being inserted I need to initialize the window object
		cout << "working on sequence " << position2contig[alignment.RefID] << "\n";
		this->resetWindow(alignment.Position, alignment.RefID); // this is executed only when the first read in inserted
	}

	if(alignment.RefID != this->chr) { // I am moving to a new chromosomes, need to check if the current window can be used or not
		if(this->windowOpen) {
			this->computeVariations();
		}
		cout << "working on sequence " << position2contig[alignment.RefID] << "\n";
		this->resetWindow(alignment.Position, alignment.RefID); // this is executed only when the first read in inserted
	}

	if(this->windowOpen) { //is window is open I need to check that I am not going out of boundaries
		if(alignment.Position > this->currentWindowEnd) { //I am out of the limit, I need to check if my current buffer contains variations
			bool varFound = this->computeVariations();
			if(varFound) {
				this->resetWindow(alignment.Position, alignment.RefID);
			} else {
				this->goToNextWindow(alignment.Position);
				if (this->currentWindowEnd < alignment.Position) {
					cout << "attenzione!!! "<<  this->currentWindowEnd  << " " << alignment.Position << "\n";
				}
			}
		}
	}

	if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
		if(alignment.RefID < alignment.MateRefID or (alignment.RefID == alignment.MateRefID and alignment.Position < alignment.MatePosition)) {  // insert only "forward" variations
			if(! this->windowOpen) { //I found a candidate but I need to open the window for it
				this->resetWindow(alignment.Position, alignment.RefID);
				this->windowOpen = true; //and set the window to open
			}
			this->TranslocationEvents.push_back(alignment);
		}
	}
	// it is a read that can be used to compute coverage
	if(this->windowOpen) { //If currently I am building a window
		if(alignmentStatus != unmapped and alignmentStatus != lowQualty ) {
			this->alignmentsOnWindow.push_back(alignment);
		}
	}


}


bool Window::computeVariations() {
	//by construction I have only forward links, i.e., from chr_i to chr_j with i<j
	this->computeCoverage();
	//if( this->coverage > 5*this->meanCoverage ) {
	//	return false;
	//}

	bool found = false;
	Translocations *Trans;
	Trans = new Translocations();
	int linksFromWindow = 0;
	for(list<BamAlignment>::iterator alignment = TranslocationEvents.begin(); alignment != TranslocationEvents.end(); ++alignment) {
		uint32_t startRead_1 			= alignment->Position;
		uint32_t chromosomeRead_1		= alignment->RefID;
		uint16_t qualityAligRead_1		= alignment->MapQuality;
		uint32_t startRead_2 			= alignment->MatePosition;
		uint32_t chromosomeRead_2		= alignment->MateRefID;
		//un-fortunatly I do not have mapping quality of read_2
		if(qualityAligRead_1 >= minimum_mapping_quality) {
			linksFromWindow ++;
			Trans->insertConnection(startRead_1, chromosomeRead_2,startRead_2);
		}
	}

	for (map<uint32_t, vector<Link> > ::iterator it1=Trans->Connections.begin(); it1!= Trans->Connections.end(); ++it1) {
		int chr2 = it1->first;
		vector<Link>  LinksToChr2 = it1->second;
		int numLinksToChr2 = LinksToChr2.size(); // number of links between chr1 and chr2 in this window (of WindowLength)
		sort(LinksToChr2.begin(), LinksToChr2.end(), sortLinksChr2); // sort all links by chr2 position
		uint32_t currentPair = 0;
		while(currentPair < numLinksToChr2) {
			vector<Link>  LinksFormingCurrentWindow; // this vector stores only current links. Mailnly used to compute real window size
			LinksFormingCurrentWindow.push_back(LinksToChr2[currentPair]); // memorize current link
			uint32_t startAt = LinksToChr2[currentPair].chr2_start;
			uint32_t stopAt  = LinksToChr2[currentPair].chr2_start + (windowSize + std_insert*10); // tollarance window
			uint32_t pairsFormingLink = 1; // number of links forming a bridge between two windows
			uint32_t nextPair = currentPair + 1;
			while(nextPair < LinksToChr2.size() and LinksToChr2[nextPair].chr2_start < stopAt ) {
				LinksFormingCurrentWindow.push_back(LinksToChr2[nextPair]);
				pairsFormingLink ++;
				nextPair ++;
			}
			uint32_t secondWindowLength = LinksToChr2[nextPair-1].chr2_start - startAt + 1 ; // windoSize (real) on second chr
			sort(LinksFormingCurrentWindow.begin(), LinksFormingCurrentWindow.end(), sortLinksChr1); // sort all links by chr1 position
			uint32_t realFirstWindowStart = LinksFormingCurrentWindow[0].chr1_start;
			uint32_t realFirstWindowEnd   = LinksFormingCurrentWindow[LinksFormingCurrentWindow.size() -1].chr1_start;
			uint32_t firstWindowLength  =  realFirstWindowEnd - realFirstWindowStart +1; // real window size on chr1
			//TODO: need to compute coverage only on window effective size
			float coverageRealFirstWindow = this->computeCoverage(realFirstWindowStart, realFirstWindowEnd);
			//ExpectedLinks(uint32_t sizeA, uint32_t sizeB, uint32_t gap, float insert_mean, float insert_stddev, float coverage, uint32_t readLength)
			//float expectedLinksInWindow = ExpectedLinks(this->windowSize, secondWindowLength, 0, mean_insert, std_insert, this->coverage, 100);
			float expectedLinksInWindow = ExpectedLinks(firstWindowLength, secondWindowLength, 0, mean_insert, std_insert, coverageRealFirstWindow, 100);
			//and pairsFormingLink/(float)linksFromWindow >= 0.2
			if( pairsFormingLink >= minimumPairs  and coverageRealFirstWindow < 5*this->meanCoverage) { //ration between coverage
				found = true;
				currentPair = nextPair + 1;
				if(this->chr == chr2) {
					intraChrVariations << position2contig[this->chr]  << "\t" <<     realFirstWindowStart   << "\t" <<       realFirstWindowEnd               << "\t"  ;
					intraChrVariations << position2contig[chr2]       << "\t" <<         startAt            << "\t" <<    LinksToChr2[nextPair-1].chr2_start  << "\t"  ;
					intraChrVariations <<      linksFromWindow        << "\t" <<        numLinksToChr2      << "\t" <<          pairsFormingLink              << "\t";
					intraChrVariations <<     coverageRealFirstWindow << "\t" <<      expectedLinksInWindow << "\t" << pairsFormingLink/expectedLinksInWindow << "\n";

				} else {
					interChrVariations << position2contig[this->chr]  << "\t" <<     realFirstWindowStart   << "\t" <<       realFirstWindowEnd               << "\t"  ;
					interChrVariations << position2contig[chr2]       << "\t" <<         startAt            << "\t" <<    LinksToChr2[nextPair-1].chr2_start  << "\t"  ;
					interChrVariations <<     linksFromWindow        << "\t" <<        numLinksToChr2      << "\t" <<          pairsFormingLink              << "\t";
					interChrVariations <<     coverageRealFirstWindow << "\t" <<      expectedLinksInWindow << "\t" << pairsFormingLink/expectedLinksInWindow << "\n";
				}

			} else {
				currentPair ++;
			}
		}
	}

	return found;

}

void Window::goToNextWindow(int nextAlignmentPosition) {
	//the window is open if I am here!!!
	int nextMinimumWindowStart = currentWindowStart + windowStep;
	while(nextMinimumWindowStart + this->windowSize < nextAlignmentPosition) {
			nextMinimumWindowStart += windowStep; // this is needed in order to avoid cases in which my next read is far away from my window end
	}

	while(TranslocationEvents.size() > 0 && TranslocationEvents.begin()->Position <  nextMinimumWindowStart) {
			TranslocationEvents.erase(TranslocationEvents.begin());
	}

	if(TranslocationEvents.size() == 0) { //if I have emptied my buffer
		resetWindow(nextMinimumWindowStart,this->chr); //I simply close the window, nothing else to do
	}  else{ //otherwise I need to remove entries from currentWindow and change the boundaries accordingly
		int nextWindowStart = TranslocationEvents.begin()->Position; //This is the new starting window point, it coincides with the first candidate read to witness a variation
		int i= 0;
		while(i < alignmentsOnWindow.size()  && alignmentsOnWindow.begin()->Position <  nextWindowStart) {
			i++; //take the pointer to the rightmost alignment
			alignmentsOnWindow.erase(alignmentsOnWindow.begin());
		}
		if(alignmentsOnWindow.size() == 0) {
			cout << "error: alignmentsOnWindow cannot be empty at this point\n";
			return;
		}

		currentWindowStart = nextWindowStart;
		currentWindowEnd   = currentWindowStart + windowSize;
		//window is still open!!!
	}
}

void Window::resetWindow(int position, uint32_t chr) {
	TranslocationEvents.clear();
	alignmentsOnWindow.clear();
	currentWindowStart 			 = position;
	currentWindowEnd   			 = position + windowSize;
	this->chr 		 = chr;
	this->windowOpen = false;

}


float Window::computeCoverage(uint32_t start, uint32_t end) {
	float totalReadSizeOnWindow = 0;
	for(list<BamAlignment>::iterator it = alignmentsOnWindow.begin(); it != alignmentsOnWindow.end(); ++it) {
		if(it->Position >= start and it->Position <= end) {
			totalReadSizeOnWindow += it->Length;
		}
	}
	return totalReadSizeOnWindow/(float)(end -start + 1);
}


float Window::computeCoverage() {
	float totalReadSizeOnWindow = 0;
	for(list<BamAlignment>::iterator it = alignmentsOnWindow.begin(); it != alignmentsOnWindow.end(); ++it) {
		totalReadSizeOnWindow += it->Length;
	}
	this->coverage = (totalReadSizeOnWindow/(float)windowSize);
	return this->coverage;
}

Translocations::Translocations() { }

void Translocations::insertConnection(uint32_t chr2, uint32_t pos2) {
	Link connection;
	connection.chr2_start = pos2;
	connection.chr2_end = pos2 + 100;
	connection.supportingPairs = 1;
	Connections[chr2].push_back(connection);
}


void Translocations::insertConnection(uint32_t pos1, uint32_t chr2, uint32_t pos2) {
	Link connection;
	connection.chr1_start = pos1;
	connection.chr2_start = pos2;
	connection.chr2_end = pos2 + 100;
	connection.supportingPairs = 1;
	Connections[chr2].push_back(connection);
}







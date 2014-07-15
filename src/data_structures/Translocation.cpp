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

bool sortLinks(Link i, Link  j) {
	return (i.chr2_start < j.chr2_start);
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
	string intra_chr_events = outputFileHeader + "_intra_chr_events.tab";
	this->intraChrVariations.open(intra_chr_events.c_str());

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
	if(alignmentStatus == unmapped or alignmentStatus == lowQualty ) {
		return; // in case the alignment is of no use discard it
	}
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
				this->goToNextWindow();
			}
		}
	}

	if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
		if(alignment.RefID < alignment.MateRefID or alignment.Position < alignment.MatePosition) {  // insert only "forward" variations
			if(! this->windowOpen) { //I found a candidate but I need to open the window for it
				this->resetWindow(alignment.Position, alignment.RefID);
				this->windowOpen = true; //and set the window to open
			}
			this->TranslocationEvents.push_back(alignment);
		}
	}
	// it is a read that can be used to compute coverage
	if(this->windowOpen) { //If currently I am building a window
		this->alignmentsOnWindow.push_back(alignment);
	}


}


bool Window::computeVariations() {
	//by construction I have only forward links, i.e., from chr_i to chr_j with i<j
	this->computeCoverage();
	//if(this->coverage > this->meanCoverage*4) {
	//	return true; // in this way I can reset it!!!! NOT PROPER WAY TO DO IT
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
			Trans->insertConnection(chromosomeRead_2,startRead_2);
		}
	}

	for (map<uint32_t, vector<Link> > ::iterator it1=Trans->Connections.begin(); it1!= Trans->Connections.end(); ++it1) {
		int chr2 = it1->first;
		vector<Link>  LinksToChr2 = it1->second;
		int supportingPairs = LinksToChr2.size(); // number of links between chr1 and chr2
		sort(LinksToChr2.begin(), LinksToChr2.end(), sortLinks);
		uint32_t currentPair = 0;
		while(currentPair < LinksToChr2.size()) {
			uint32_t startAt = LinksToChr2[currentPair].chr2_start;
			uint32_t stopAt  = LinksToChr2[currentPair].chr2_start + (windowSize + std_insert*10);
			uint32_t pairsInWindow = 1; // number of links forming a bridge between two windows
			uint32_t nextPair = currentPair + 1;
			while(nextPair < LinksToChr2.size() and LinksToChr2[nextPair].chr2_start < stopAt ) {
				pairsInWindow ++;
				nextPair ++;
			}
			uint32_t secondWindowLength = LinksToChr2[nextPair-1].chr2_start - startAt;
			//ExpectedLinks(uint32_t sizeA, uint32_t sizeB, uint32_t gap, float insert_mean, float insert_stddev, float coverage, uint32_t readLength)
			float expectedLinksInWindow = ExpectedLinks(this->windowSize, secondWindowLength, 0, mean_insert, std_insert, this->coverage, 100);

			float coverage1 =  (float)(supportingPairs*100)/(float)(windowSize);
			float coverage2 =  (float)(pairsInWindow*100)/(float)(secondWindowLength);
			//&& (this->coverage/meanCoverage < 3) && (coverage1/meanCoverage > 0.2)
			if( pairsInWindow >= minimumPairs and pairsInWindow/(float)linksFromWindow >= 0.2) { //ration between coverage
				found = true;
				currentPair = nextPair + 1;
				if(this->chr == chr2) {
					intraChrVariations << position2contig[this->chr] << "\t" << this->currentWindowStart  << "\t" << this->currentWindowEnd              << "\t"   << supportingPairs  << "\t";
					intraChrVariations << position2contig[chr2]      << "\t" <<            startAt        << "\t" << LinksToChr2[nextPair-1].chr2_start  << "\t"   <<  pairsInWindow   <<"\t";
					intraChrVariations << expectedLinksInWindow      << "\t" << this->coverage            << "\t" << linksFromWindow << "\n";
				} else {
					interChrVariations << position2contig[this->chr] << "\t" << this->currentWindowStart  << "\t" << this->currentWindowEnd              << "\t"   << supportingPairs  << "\t";
					interChrVariations << position2contig[chr2]      << "\t" <<            startAt        << "\t" << LinksToChr2[nextPair-1].chr2_start  << "\t"   <<  pairsInWindow   <<"\t";
					interChrVariations << expectedLinksInWindow      << "\t" << this->coverage            << "\t" << linksFromWindow << "\n";
				}

			} else {
				currentPair ++;
			}
		}
	}

	return found;

}

void Window::goToNextWindow() {
	//the window is open if I am here!!!
	int nextMinimumWindowStart = currentWindowStart + windowStep;
//TODO try the same trick

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








/*
    FRC: computes the FRC curve starting from alignments
    Copyright (C) 2011  F. Vezzi(vezi84@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include "radix.h"
#include "samtools/sam.h"
#include "options/Options.h"




class windowStatistics {
public:
	windowStatistics() {
		windowLength = 0;
		readsLength_win=0;
		insertsLength_win=0;
		correctlyMatedReadsLength_win=0;
		wronglyOrientedReadsLength_win=0;
		singletonReadsLength_win=0;
		matedDifferentContigLength_win=0;
	}

	void reset() {
		windowLength = 0;
		readsLength_win=0;
		inserts=0;
		insertsLength_win=0;
		correctlyMatedReadsLength_win=0;
		wronglyDistanceReadsLength_win = 0;
		wronglyOrientedReadsLength_win=0;
		singletonReadsLength_win=0;
		matedDifferentContigLength_win=0;

	}

	void print() {
		cout << "from " << windowStart << " to " << windowEnd << "\n";
		cout << "windowLength " << windowLength << "\n";
		cout << "readsLength " << readsLength_win << "\n";
		cout << "inserts " << inserts << "\n";
		cout << "insersLength " << insertsLength_win << "\n";
		cout << "correctlyMatedReadsLength " << correctlyMatedReadsLength_win << "\n";
		cout << "wronglyOrientedReadsLength " << wronglyOrientedReadsLength_win << "\n";
		cout << "singletonReadsLength " << singletonReadsLength_win << "\n";
		cout << "matedDifferentContigLength " << matedDifferentContigLength_win  << "\n";
		cout << "-----\n";
	}

	uint32_t windowStart;
	uint32_t windowEnd;
	uint32_t windowLength;

	// reads aligned
	uint64_t readsLength_win; // length of reads placed in window
	// insert length
	uint32_t inserts; // number of inserts inside window
	uint64_t insertsLength_win; // total length of inserts inside window
	// correctly aligned mates
	uint64_t correctlyMatedReadsLength_win; // length of correctly mated reads inside window
	// wrongly oriented reads
	uint64_t wronglyOrientedReadsLength_win; // length of wrongly oriented reads inside window
	// wrongly distance reads
	uint64_t wronglyDistanceReadsLength_win; // total length of reads placed in different contigs  inside window
	// singletons
	uint64_t singletonReadsLength_win; // total length of singleton reads  inside window
	// mates on different contigs
	uint64_t matedDifferentContigLength_win;// total number of reads placed in different contigs  inside window


/*
	float C_A; // total read coverage
	float S_A; // total span coverage
	float C_M; // coverage induced by correctly aligned pairs
	float C_W; // coverage induced by wrongly mated pairs
	float C_S; // coverage induced by singletons
	float C_C; // coverage induced by reads with mate on a different contig
	float CEstatistics;
	*/
};


class ContigFeature {
public:
	ContigFeature() {
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
	}

	void updateLOW_COVERAGE_AREA() {
		LOW_COVERAGE_AREA++;
		TOTAL++;
	}
	void updateHIGH_COVERAGE_AREA() {
		HIGH_COVERAGE_AREA++;
		TOTAL++;
	}
	void updateLOW_NORMAL_AREA() {
		LOW_NORMAL_AREA++;
		TOTAL++;
	}
	void updateHIGH_NORMAL_AREA() {
		HIGH_NORMAL_AREA++;
		TOTAL++;
	}
	void updateHIGH_SINGLE_AREA() {
		HIGH_SINGLE_AREA++;
		TOTAL++;
	}
	void updateHIGH_SPANING_AREA() {
		HIGH_SPANING_AREA++;
		TOTAL++;
	}
	void updateHIGH_OUTIE() {
		HIGH_OUTIE++;
		TOTAL++;
	}
	void updateCOMPRESSION_AREA() {
		COMPRESSION_AREA++;
		TOTAL++;
	}
	void updateSTRECH_AREA() {
		STRECH_AREA++;
		TOTAL++;
	}


/*
	void updateFeatures(windowStatistics* window) {
		float C_A_i = window->readsLength_win/(float)window->insertsLength_win;
		float S_A_i = window->insertsLength_win/(float)window->insertsLength_win;
		float C_M_i = window->correctlyMatedReadsLength_win/(float)window->insertsLength_win;
		float C_W_i = (window->wronglyDistanceReadsLength_win + window->wronglyOrientedReadsLength_win)/(float)window->insertsLength_win;
		float C_S_i = window->singletonReadsLength_win/(float)window->insertsLength_win;
		float C_C_i = window->matedDifferentContigLength_win/(float)window->insertsLength_win;


	}
*/

	void print() {
		cout << "contigLength " << contigLength << "\n";
		cout << "LOW_COVERAGE_AREA " << LOW_COVERAGE_AREA << "\n";
		cout << "HIGH_COVERAGE_AREA " << HIGH_COVERAGE_AREA << "\n";
		cout << "LOW_NORMAL_AREA " << LOW_NORMAL_AREA << "\n";
		cout << "HIGH_NORMAL_AREA " << HIGH_NORMAL_AREA << "\n";
		cout << "HIGH_SINGLE_AREA " << HIGH_SINGLE_AREA << "\n";
		cout << "HIGH_SPANING_AREA " << HIGH_SPANING_AREA << "\n";
		cout << "HIGH_OUTIE " << HIGH_OUTIE << "\n";
		cout << "COMPRESSION_AREA " << COMPRESSION_AREA << "\n";
		cout << "STRECH_AREA " << STRECH_AREA << "\n-----\n";
	}

	uint64_t contigLength;
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
};


float lowCoverageFeat = 1/(float)3;
float highCoverageFeat = 3;
float lowNormalFeat = 1/(float)3;
float highNormalFeat = 3;
float highSingleFeat = 0.6;
float highSpanningFeat = 0.6;
float highOutieFeat = 0.6;
float CE_statistics = 3;

bool sortContigs(ContigFeature i, ContigFeature j) {
	return (i.contigLength > j.contigLength);
}

#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)

#define EXIT_IF_NULL(P) \
  if (P == NULL) \
    return 1;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core)
{

  if (core->flag&BAM_FUNMAP) {
    return false;
  }

  return true;
}


/**
 * Open a .sam/.bam file.
 * @returns NULL is open failed.
 */
samfile_t * open_alignment_file(std::string path)
{
  samfile_t * fp = NULL;
  std::string flag = "r";
  if (path.substr(path.size()-3).compare("bam") == 0) {
    //BAM file!
    flag += "b";
  }
  if ((fp = samopen(path.c_str(), flag.c_str() , 0)) == 0) {
    fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
  }
  return fp;
}


void updateWindow(bam1_t* b, windowStatistics* actualWindow, unsigned int peMinInsert, unsigned int peMaxInsert ) {
	const bam1_core_t* core =  &b->core;
	if (core->qual >= 0) {
		uint32_t* cigar = bam1_cigar(b);

		if ((core->flag&BAM_FREAD1) //First in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			//pair is mapped on the same contig and I'm looking the first pair
			int32_t start = MIN(core->pos,core->mpos);
			int32_t end = start+abs(core->isize);
			int32_t iSize = end-start; // compute insert size



			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				actualWindow->insertsLength_win += iSize; // update number of inserts and total length
				actualWindow->inserts++;
			}

			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				if(core->pos < core->mpos) { // read I'm processing is the first
					if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
					} else {
						// wrong orientation
						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				} else {
					if(!(core->flag&BAM_FMREVERSE) && (core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
					} else {
						// wrong orientation
						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				}
			} else {
				//wrong distance
				actualWindow->wronglyDistanceReadsLength_win += bam_cigar2qlen(core,cigar);
			}
		} else  if ((core->flag&BAM_FREAD2) //Second in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		)
			// if I'm considering the second read in a pair I must check it is is a correctly mated read and if this is tha case update the right variables
		{
			int32_t start = MIN(core->pos,core->mpos);
			int32_t end = start+abs(core->isize);
			int32_t iSize = end-start; // compute insert size

			if (peMinInsert <= iSize && iSize <= peMaxInsert) { // I have to check if the mate is outside window boundaries
				if(start <= actualWindow->windowStart || end >= actualWindow->windowEnd) {
					actualWindow->insertsLength_win += iSize; // update number of inserts and total length
					actualWindow->inserts++;
				}
			}

			if (peMinInsert <= iSize && iSize <= peMaxInsert) {
				if(core->pos > core->mpos) { // read I'm processing is the first
					if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); //  update number of correctly mapped and their length
					} else {
						// wrong orientation
						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				} else {
					if((core->flag&BAM_FMREVERSE) && !(core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
						actualWindow->correctlyMatedReadsLength_win +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
					} else {
						// wrong orientation
						actualWindow->wronglyOrientedReadsLength_win += bam_cigar2qlen(core,cigar);
					}
				}
			} else {
				//wrong distance
				actualWindow->wronglyDistanceReadsLength_win += bam_cigar2qlen(core,cigar);

			}
		} else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
			//Count inter-chrom pairs
			actualWindow->matedDifferentContigLength_win += bam_cigar2qlen(core,cigar);
		} else if(core->flag&BAM_FMUNMAP) {
			// if mate read is unmapped
			actualWindow->singletonReadsLength_win =+ bam_cigar2qlen(core,cigar);
		}


		if (core->flag&BAM_FDUP) {   //This is a duplicate. Don't count it!.
		} else {
			actualWindow->readsLength_win += bam_cigar2qlen(core,cigar);
		}
	}
}




/**
 * Main of app
 */


int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string alignmentFile = "";
	string assemblyFile = "";
	vector<string> peFiles;
	uint32_t peMinInsert = 100;
	uint32_t peMaxInsert = 1000000;
	unsigned int WINDOW = 1000;
	unsigned long int estimatedGenomeSize;

	samfile_t *fp;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
			("assembly", po::value<string>(), "assembly file name in fasta format")
			("sam", po::value<string>(), "alignment file")
			("pe", po::value< vector < string > >(), "paired reads, one pair after the other")
			("pe-min-insert",  po::value<int>(), "minimum allowed insert size")
			("pe-max-insert",  po::value<int>(), "maximum allowed insert size")
			("window",  po::value<unsigned int>(), "window size for features computation")
			("genome-size", po::value<unsigned long int>(), "estimated genome size (if not supplied genome size is belived to be assembly length")
			;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(2);
	}
	if (vm.count("help")) {
		DEFAULT_CHANNEL << desc << endl;
		exit(0);
	}
	if (!vm.count("assembly") || !vm.count("sam") || !vm.count("pe")) {
		DEFAULT_CHANNEL << "--assembly, --sam and,  --pe are mandatory" << endl;
		exit(0);
	}

	if (vm.count("assembly")) {
		assemblyFile = vm["assembly"].as<string>();
	}
	if (vm.count("sam")) {
		alignmentFile = vm["sam"].as<string>();
	}
	if (vm.count("pe")) {
		peFiles = vm["pe"].as<vector<string> >();
	}

	if (vm.count("pe-min-insert")) {
		peMinInsert = vm["pe-min-insert"].as<int>();
	}

	if (vm.count("pe-max-insert")) {
		peMaxInsert = vm["pe-max-insert"].as<int>();
	}

	if (vm.count("window")) {
		WINDOW = vm["window"].as<unsigned int>();
	}

	if (vm.count("genome-size")) {
		estimatedGenomeSize = vm["genome-size"].as<unsigned long int>();
	} else {
		estimatedGenomeSize = 0;
	}

	cout << "assembly file name is " << assemblyFile << endl;
	cout << "sam file name is " << alignmentFile << endl;
	for(int i=0; i< peFiles.size() ; i++ ) {
		cout << peFiles.at(i) << endl;
	}

	fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);

	//Initialize bam entity
	bam1_t *b = bam_init1();

	//All var declarations
	uint64_t genomeLength = 0; // total genome length
	uint32_t contigs = 0; // number of contigs/scaffolds
// total reads
	uint32_t reads = 0;
	uint64_t readsLength = 0;   // total length of reads
	uint32_t unmappedReads = 0;
	uint32_t mappedReads = 0;
	uint32_t zeroQualityReads = 0;

	uint32_t contigSize = 0;
	uint32_t duplicates = 0;

	uint32_t inserts = 0;       // total number of inserts
	uint64_t insertsLength = 0; // total inserts length
	float insertMean;
	float insertStd;

// mated reads (not necessary correctly mated)
	uint32_t matedReads = 0;        // length of reads that align on a contig with the mate
	uint64_t matedReadsLength = 0;  // total length of mated reads

 // correctly aligned mates
	uint32_t correctlyMatedReads = 0; // total number of correctly mated reads
	uint64_t correctlyMatedReadsLength = 0; // length of correctly mated reads

// wrongly oriented reads
	uint32_t wronglyOrientedReads = 0; // number of wrongly oriented reads
	uint64_t wronglyOrientedReadsLength = 0; // length of wrongly oriented reads

// wrongly distance reads
	uint32_t wronglyDistanceReads       = 0; // number of reads at the wrong distance
	uint64_t wronglyDistanceReadsLength = 0;  // total length of reads placed in different contigs


// singletons
	uint32_t singletonReads = 0; // number of singleton reads
	uint64_t singletonReadsLength = 0;     // total length of singleton reads

// mates on different contigs
	uint32_t matedDifferentContig = 0; // number of contig placed in a different contig
	uint64_t matedDifferentContigLength = 0; // total number of reads placed in different contigs

	float C_A = 0; // total read coverage
	float S_A = 0; // total span coverage

	float C_M = 0; // coverage induced by correctly aligned pairs
	float C_W = 0; // coverage induced by wrongly mated pairs
	float C_S = 0; // coverage induced by singletons
	float C_C = 0; // coverage induced by reads with mate on a diferent contif


// compute mean and std on the fly
	float Mk = 0;
	float Qk = 0;
	uint32_t counterK = 1;


    //Keep header for further reference
    bam_header_t* head = fp->header;
    int32_t currentTid = -1;


    while (samread(fp, b) >= 0) {
	      //Get bam core.
	      const bam1_core_t *core = &b->core;

	      if (core == NULL) {
	    	  //There is something wrong with the read/file
	    	  printf("Input file is corrupt!");
	    	  //Leak everything and exit!
	    	  return -1;
	      }
	      ++reads;

	      if (!is_mapped(core)) {
	    	  ++unmappedReads;
	      } else {
	    	  if (core->tid != currentTid) {
	    		  if (currentTid != -1) {
	    			  //			  	 cout << "working on new contig: " << head->target_name[core->tid] << endl;
	    		  }
	    		  //Get length of next section
	    		  contigSize = head->target_len[core->tid];
	    		  contigs++;
	    		  if (contigSize < 1) {//We can't have such sizes! this can't be right
	    			  fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid],contigSize);
	    		  }
	    		  genomeLength += contigSize;
	    		  currentTid = core->tid;
	    	  }

	    	  //If read has quality == 0, we won't count it as mapped
	    	  if (core->qual >= 0) {
	    		  uint32_t* cigar = bam1_cigar(b);
	    		  ++mappedReads;
	    		  //Now check if reads belong to a proper pair: both reads aligned on the same contig at the expected distance and orientation

	    		  if ((core->flag&BAM_FREAD1) //First in pair
	    				  && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
	    				  && (core->tid == core->mtid) /*Mate on the same chromosome*/
	    		  ) {
	    			  //pair is mapped on the same contig and I'm looking the first pair
	    			  int32_t start = MIN(core->pos,core->mpos);
	    			  int32_t end = start+abs(core->isize);
	    			  int32_t iSize = end-start; // compute insert size
	    			  inserts ++;
	    			  insertsLength += iSize; // update number of inserts and total length

	    			  if (peMinInsert <= iSize && iSize <= peMaxInsert) {
	    				  if(counterK == 1) {
	    					  Mk = iSize;
	    					  Qk = 0;
	    					  counterK++;
	    				  } else {
	    					  float oldMk = Mk;
	    					  float oldQk = Qk;
	    					  Mk = oldMk + (iSize - oldMk)/counterK;
	    					  Qk = oldQk + (counterK-1)*(iSize - oldMk)*(iSize - oldMk)/(float)counterK;
	    					  counterK++;
	    				  }
	    			  }

	    			  if (peMinInsert <= iSize && iSize <= peMaxInsert) {
	    				  if(core->pos < core->mpos) { // read I'm processing is the first
	    					  if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  // wrong orientation
	    						  wronglyOrientedReads++;
	    						  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  if(!(core->flag&BAM_FMREVERSE) && (core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  // wrong orientation
	    						  wronglyOrientedReads++;
	    						  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  }
	    			  } else {
	    				  //wrong distance
	    				  wronglyDistanceReads++;
	    				  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
	    			  }
	    		  } else  if ((core->flag&BAM_FREAD2) //Second in pair
	    				  && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
	    				  && (core->tid == core->mtid) /*Mate on the same chromosome*/
	    		  )
	    	// if I'm considering the second read in a pair I must check it is is a correctly mated read and if this is tha case update the right variables
	    			  {
	    			  int32_t start = MIN(core->pos,core->mpos);
	    			  int32_t end = start+abs(core->isize);
	    			  int32_t iSize = end-start; // compute insert size
	    			  if (peMinInsert <= iSize && iSize <= peMaxInsert) {
	    				  if(core->pos > core->mpos) { // read I'm processing is the first
	    					  if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { // pairs are one in front of the other
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); //  update number of correctly mapped and their length
	    					  } else {
	    						  // wrong orientation
	    						  wronglyOrientedReads++;
	    						  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  if((core->flag&BAM_FMREVERSE) && !(core->flag&BAM_FREVERSE)) { // pairs are one in front of the other
	    						  correctlyMatedReads ++;
  	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  // wrong orientation
	    						  wronglyOrientedReads++;
	    						  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  }
	    			  } else {
	    				  //wrong distance
	    				  wronglyDistanceReads++;
	    				  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);

	    			  }
	    		  } else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
	    			  //Count inter-chrom pairs
	    			  matedDifferentContig++;
	    			  matedDifferentContigLength += bam_cigar2qlen(core,cigar);
	    		  } else if(core->flag&BAM_FMUNMAP) {
		    		  // if mate read is unmapped
	    			  singletonReads++;
	    			  singletonReadsLength =+ bam_cigar2qlen(core,cigar);
	    		  }


	    		  if (core->flag&BAM_FPROPER_PAIR) {
	    			  //Is part of a proper pair
	    			  matedReads ++; // increment number of mated reads
	    			  matedReadsLength += bam_cigar2qlen(core,cigar); // add the length of the read aligne as proper mate (not necessary correctly mated)
	    		  }

	    		  if (core->flag&BAM_FDUP) {   //This is a duplicate. Don't count it!.
	    			  ++duplicates;
	    		  } else {
	    			  readsLength += bam_cigar2qlen(core,cigar);
	    		  }
	    	  } else {
	    		  ++zeroQualityReads;

	    	  }
	      }
    }

    genomeLength += contigSize;
    if (estimatedGenomeSize == 0) {
    	estimatedGenomeSize = genomeLength;
    }
    cout << "BAM file has been read and statistics have been computed\n";

    cout << "total number of contigs " << contigs << endl;
    cout << "total reads number " << reads << "\n";
    cout << "total mapped reads " << mappedReads << "\n";
    cout << "total unmapped reads " << unmappedReads << "\n";
    cout << "proper pairs " << matedReads << "\n";
    cout << "zero quality reads " << zeroQualityReads << "\n";
    cout << "correctly oriented " << correctlyMatedReads << "\n";
    cout << "wrongly oriented " << wronglyOrientedReads << "\n";
    cout << "wrongly distance " << wronglyDistanceReads << "\n";
    cout << "wrongly contig " <<  matedDifferentContig << "\n";
    cout << "singletons " << singletonReads << "\n";

    uint32_t total = correctlyMatedReads + wronglyOrientedReads + wronglyDistanceReads + matedDifferentContig + singletonReads;
    cout << "total " << total << "\n";
    cout << "\n-------\n";

    C_A = readsLength/(float)genomeLength;
    S_A = insertsLength/(float)genomeLength;
    C_M = correctlyMatedReadsLength/(float)genomeLength;
    C_W = (wronglyDistanceReadsLength + wronglyOrientedReadsLength)/(float)genomeLength;
    C_S = (singletonReadsLength)/(float)genomeLength;
    C_C = matedDifferentContigLength/(float)genomeLength;


    cout << "C_A = " << C_A << endl;
    cout << "S_A = " << S_A << endl;
    cout << "C_M = " << C_M << endl;
    cout << "C_W = " << C_W << endl;
    cout << "C_S = " << C_S << endl;
    cout << "C_C = " << C_C << endl;

    cout << "\n";
    cout << "Mean Insert length = " << Mk << endl;
    insertMean = Mk;
    Qk = sqrt(Qk/counterK);
    cout << "Std Insert length = " << Qk << endl;
    insertStd = Qk;




    //now close file and parse it again to compute FRC curve
    samclose(fp);

    cout << "\n----------\nNow computing FRC \n------\n";

    fp = open_alignment_file(alignmentFile);
    EXIT_IF_NULL(fp);
    //Keep header for further reference
    head = fp->header;
    currentTid = -2;
    reads = 0;

    vector<ContigFeature> CONTIG (contigs);
    uint32_t featuresTotal = 0;

//    ContigFeature *CONTIG = new ContigFeature[contigs];
    uint32_t contig=0;

	uint32_t windowStart = 0;
	uint32_t windowEnd   = windowStart + WINDOW;
	windowStatistics* actualWindow = new windowStatistics();
/*	actualWindow->C_A = C_A;
	actualWindow->S_A = S_A;
	actualWindow->C_M = C_M;
	actualWindow->C_W = C_W;
	actualWindow->C_S = C_S;
	actualWindow->C_C = C_C;
	actualWindow->CEstatistics = 3;
*/
	actualWindow->windowLength = (windowEnd -windowStart + 1);
	actualWindow->windowStart = windowStart;
	actualWindow->windowEnd = windowEnd;

   	float C_A_i = 0;  // read coverage of window
   	float S_A_i = 0; // span coverage of window
   	float C_M_i = 0; // coverage of correctly aligned reads of window
   	float C_W_i = 0; // coverage of wrongly aligned reads
   	float C_S_i = 0; // singleton coverage of window
   	float C_C_i = 0; // coverage of reads with mate on different contigs
   	float Z_i   = 0; // CE statistics

    while (samread(fp, b) >= 0) {
    	//Get bam core.
    	const bam1_core_t *core = &b->core;
    	if (core == NULL) {
     		printf("Input file is corrupt!");
     		return -1;
    	}
    	++reads;

    	// new contig
    	if (!is_mapped(core)) {
    		++unmappedReads;
    	} else {
    		if (core->tid != currentTid) { // another contig or simply the first one
    			if(currentTid == -2) {
    				cout << "first contig read\n";
    			}
    			if(currentTid > -1) {
    				//compute statistics and update contig feature for the last segment of the contig
    				//actualWindow->print();
        			C_A_i = actualWindow->readsLength_win/(float)actualWindow->windowLength;  // read coverage of window
        		   	S_A_i = actualWindow->insertsLength_win/(float)actualWindow->windowLength; // span coverage of window
        		   	C_M_i = actualWindow->correctlyMatedReadsLength_win/(float)actualWindow->windowLength; // coverage of correctly aligned reads of window
        		   	C_W_i = (actualWindow->wronglyDistanceReadsLength_win+actualWindow->wronglyOrientedReadsLength_win)/(float)actualWindow->windowLength; // coverage of wrongly aligned reads
        		   	C_S_i = actualWindow->singletonReadsLength_win/(float)actualWindow->windowLength; // singleton coverage of window
        		   	C_C_i = actualWindow->matedDifferentContigLength_win/(float)actualWindow->windowLength; // coverage of reads with mate on different contigs
        		   	if(actualWindow->inserts > 0) {
        		   		float localMean = actualWindow->insertsLength_win/(float)actualWindow->inserts;
        		   		Z_i   = (localMean - insertMean)/(float)(insertStd/sqrt(actualWindow->inserts)); // CE statistics
        		   	} else {
        		   		Z_i = -100;
        		   	}
        			//actualWindow->print();
        		   	// NOW UPDATE CONTIG'S FEATURES
        		   	if(C_A_i < lowCoverageFeat*C_A) {
        		   		CONTIG[contig].updateLOW_COVERAGE_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_A_i > highCoverageFeat*C_A) {
        		   		CONTIG[contig].updateHIGH_COVERAGE_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_M_i <  lowNormalFeat*C_M) {
        		   		CONTIG[contig].updateLOW_NORMAL_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_M_i > highNormalFeat*C_M) {
        		   		CONTIG[contig].updateHIGH_NORMAL_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_S_i > highSingleFeat*C_A_i) {
        		   		CONTIG[contig].updateHIGH_SINGLE_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_C_i > highSpanningFeat*C_A_i) {
        		   		CONTIG[contig].updateHIGH_SPANING_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(C_W_i > highOutieFeat*C_A_i) {
        		   		CONTIG[contig].updateHIGH_OUTIE();
        		   		featuresTotal++;
        		   	}
        		   	if(Z_i < -CE_statistics) {
        		   		CONTIG[contig].updateCOMPRESSION_AREA();
        		   		featuresTotal++;
        		   	}
        		   	if(Z_i > CE_statistics) {
        		   		CONTIG[contig].updateSTRECH_AREA();
        		   		featuresTotal++;
        		   	}
        		   	//        		   	CONTIG[contig].print();

        		   	actualWindow->reset();
        		   	contig++; // if I'm not reading the first contig update the contig counter
    			}
    			//Get length of next section
    			contigSize = head->target_len[core->tid];
    			CONTIG[contig].contigLength = contigSize; // new contig (0 at first iteration)

    			if (contigSize < 1) {//We can't have such sizes! this can't be right
    				fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid],contigSize);
    			}
    			currentTid = core->tid;
    			actualWindow->windowStart = 0; // reset window start
    			actualWindow->windowEnd  = actualWindow->windowStart + WINDOW; // reset window end
    			if(actualWindow->windowEnd > contigSize) {
    				actualWindow->windowEnd = contigSize;
    			}
    			actualWindow->windowLength = (actualWindow->windowEnd - actualWindow->windowStart + 1); // set window length

    		}

    		if (core->pos > actualWindow->windowEnd) {
    			//compute statistics and update contig feature
    			C_A_i = actualWindow->readsLength_win/(float)actualWindow->windowLength;  // read coverage of window
    		   	S_A_i = actualWindow->insertsLength_win/(float)actualWindow->windowLength; // span coverage of window
    		   	C_M_i = actualWindow->correctlyMatedReadsLength_win/(float)actualWindow->windowLength; // coverage of correctly aligned reads of window
    		   	C_W_i = (actualWindow->wronglyDistanceReadsLength_win+actualWindow->wronglyOrientedReadsLength_win)/(float)actualWindow->windowLength; // coverage of wrongly aligned reads
    		   	C_S_i = actualWindow->singletonReadsLength_win/(float)actualWindow->windowLength; // singleton coverage of window
    		   	C_C_i = actualWindow->matedDifferentContigLength_win/(float)actualWindow->windowLength; // coverage of reads with mate on different contigs
    		   	if(actualWindow->inserts > 0) {
    		   		float localMean = actualWindow->insertsLength_win/(float)actualWindow->inserts;
    		   		Z_i   = (localMean - insertMean)/(float)(insertStd/sqrt(actualWindow->inserts)); // CE statistics
    		   	} else {
    		   		Z_i = -100;
    		   	}
    			//actualWindow->print();
    		   	// NOW UPDATE CONTIG'S FEATURES
    		   	if(C_A_i < lowCoverageFeat*C_A) {
    		   		CONTIG[contig].updateLOW_COVERAGE_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_A_i > highCoverageFeat*C_A) {
    		   		CONTIG[contig].updateHIGH_COVERAGE_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_M_i <  lowNormalFeat*C_M) {
    		   		CONTIG[contig].updateLOW_NORMAL_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_M_i > highNormalFeat*C_M) {
    		   		CONTIG[contig].updateHIGH_NORMAL_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_S_i > highSingleFeat*C_A_i) {
    		   		CONTIG[contig].updateHIGH_SINGLE_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_C_i > highSpanningFeat*C_A_i) {
    		   		CONTIG[contig].updateHIGH_SPANING_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(C_W_i > highOutieFeat*C_A_i) {
    		   		CONTIG[contig].updateHIGH_OUTIE();
    		   		featuresTotal++;
    		   	}
    		   	if(Z_i < -CE_statistics) {
    		   		CONTIG[contig].updateCOMPRESSION_AREA();
    		   		featuresTotal++;
    		   	}
    		   	if(Z_i > CE_statistics) {
    		   		CONTIG[contig].updateSTRECH_AREA();
    		   		featuresTotal++;
    		   	}

    			actualWindow->reset();
    			actualWindow->windowStart = actualWindow->windowEnd + 1;
    			actualWindow->windowEnd += WINDOW;
    			if(actualWindow->windowEnd > contigSize) {
    				actualWindow->windowEnd = contigSize;
    			}
    			actualWindow->windowLength = (actualWindow->windowEnd - actualWindow->windowStart + 1); // set window length
    			//actualWindow->windowStart = windowStart;
    			//actualWindow->windowEnd = windowEnd;
    			updateWindow(b, actualWindow,  peMinInsert, peMaxInsert);
    		} else {
    			//this read contributes to the features of the current window
    			updateWindow(b, actualWindow,  peMinInsert, peMaxInsert);
    		}


    	}


    }
 // TODO: UPDATE LAST CONTIG

    cout << "number of reads " << reads << "\n";
    //SORT CONTIGS ACCORDING TO THEIR LENGTH
    sort(CONTIG.begin(),CONTIG.end(),sortContigs);
//    for(uint32_t i=0; i< contigs; i++) {
//   	CONTIG[i].print();
//    }

    cout << "Now computing FRC \n";
    ofstream myfile;
    myfile.open ("FRC.txt");
    myfile << "features coverage\n";
    float step = featuresTotal/(float)100;
    float partial=0;
    while(partial <= featuresTotal) {
    	uint32_t contigStep = 0;
    	uint64_t contigLengthStep = 0;
    	uint32_t featuresStep = 0;
    	while(featuresStep <= partial) {
    		contigLengthStep += CONTIG[contigStep].contigLength;
    		featuresStep += CONTIG[contigStep].TOTAL;
    		contigStep++;
    	}
    	float coveragePartial =  contigLengthStep/(float)estimatedGenomeSize;
    	myfile << partial << " " << coveragePartial << "\n";
    	partial += step;
    }
    myfile.close();


}




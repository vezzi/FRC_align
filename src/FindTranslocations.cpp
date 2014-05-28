/*
   Francesco Vezzi
 */

#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <string>

#include <sstream>
#include <iostream>
#include <fstream>

//#include "common.h"
#include "data_structures/Translocations.h"
#include  <boost/program_options.hpp>
namespace po = boost::program_options;


void findTranslocations(string file, int32_t min_insert,  int32_t max_insert, bool outtie, uint64_t genomeSize,
		uint16_t minimum_mapping_quality, uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs, float coverage,
		string outputFileHeader);

int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string alignmentFile		    = "";       // alignment file name
	bool outtie 				    = true;	 // library orientation
	uint32_t windowSize 		    = 1000;     // Window size
	uint32_t windowStep 		    = 100;     // Window step
	uint32_t minimumSupportingPairs = 10;
	int min_insert				    = 100;      // min insert size
	int max_insert				    = 1000000;  // max insert size
	string outputFileHeader         = "output"; // default output name
	int minimum_mapping_quality     = 20;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
						("bam", po::value<string>(), "alignment file in bam format, expected sorted by read name. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads")
						("min-insert",  po::value<int>(), "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read")
						("max-insert",  po::value<int>(), "paired reads maximum allowed insert size. Used in order to filter outliers.")
						("orientation", po::value<string>(), "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie")
						("output",  po::value<string>(), "Header output file names")
						("minimum-supporting-pairs",  po::value<unsigned int>(), "Minimum number of supporting pairs in order to call a variation event (default 10)")
						("minimum-mapping-quality",  po::value<int>(), "Minimum mapping quality to consider an alignment (default 20)")
						("window-size",  po::value<unsigned int>(), "Size of the sliding window (default 1000)")
						("window-step",  po::value<unsigned int>(), "size of the step in overlapping window (must be lower than window-size) (default 100)")
						;
//TODO: add minimum number of reads to support a translocation, window size etc.

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

	// PARSE SAM/BAM file
	if (!vm.count("bam")) {
		DEFAULT_CHANNEL << "Please specify --bam " << endl;
		exit(0);
	}
	if (vm.count("bam")) {
		alignmentFile = vm["bam"].as<string>();
	}
	if (!vm.count("min-insert") || !vm.count("max-insert")) {
		DEFAULT_CHANNEL << "Please specify min-insert and max-insert " << endl;
		exit(0);
	}
	if (vm.count("min-insert")) {
		min_insert = vm["min-insert"].as<int>();
		if(min_insert <= 0) {
			DEFAULT_CHANNEL << "minimum insert should be at least 1\n";
			DEFAULT_CHANNEL << desc << endl;
			exit(2);
		}
	}
	if (vm.count("max-insert")) {
		max_insert = vm["max-insert"].as<int>();
	}
	if (vm.count("orientation")) {
		string userOrientation = vm["orientation"].as<string>();
		if(userOrientation.compare("innie") == 0) {
			outtie = false;
		} else if (userOrientation.compare("outtie") == 0) {
			outtie = true;
		} else {
			DEFAULT_CHANNEL << "outtie (<- ->) or innie (-> <-) only allowed orientations\n";
			DEFAULT_CHANNEL << desc << endl;
			exit(2);
		}
	}

	if (vm.count("output")) {
		string header = vm["output"].as<string>();
		outputFileHeader = header ;
	}

	if (vm.count("minimum-mapping-quality")) {
		minimum_mapping_quality = vm["minimum-mapping-quality"].as<int>();
	}

	if (vm.count("window-size")) {
		windowSize = vm["window-size"].as<unsigned int>();
	}

	if (vm.count("window-step")) {
		windowStep = vm["window-step"].as<unsigned int>();
		if (windowStep > windowSize) {
			DEFAULT_CHANNEL << "window-step cannot be larger than window-size\n";
			exit(2);
		}
	}

	if (vm.count("minimum-supporting-pairs")) {
		minimumSupportingPairs = vm["minimum-supporting-pairs"].as<unsigned int>();
		}





	if(vm.count("sam")){
		cout << "sam file name is " << alignmentFile << endl;
		cout << "library min " << min_insert << "\n";
		cout << "library max " << max_insert << "\n";
		if(outtie) {
			cout << "library orientation <- ->\n";
		} else {
			cout << "library orientation -> <-\n";
		}
	}


	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;
	BamReader bamFile;
	bamFile.Open(alignmentFile);



	SamHeader head = bamFile.GetHeader();
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		genomeLength += StringToNumber(sequence->Length);
		contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
		position2contig[contigsNumber] = contig2position[sequence->Name];
		contigsNumber++;
	}
	bamFile.Close();

	cout << "total number of contigs " 	<< contigsNumber << endl;
	cout << "assembly length " 			<< genomeLength << "\n";


	LibraryStatistics library;
	library = computeLibraryStats(alignmentFile, genomeLength, max_insert, outtie);

	float coverage = library.C_A;
	findTranslocations(alignmentFile, min_insert, max_insert, outtie, genomeLength,
			minimum_mapping_quality, windowSize, windowStep, minimumSupportingPairs, coverage, outputFileHeader);
	// now find translocations
	/*
	 *
	 fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);
	head = fp->header; // sam header
	findTranslocations(outputFileDescriptor, fp,MinInsert, MaxInsert, estimatedGenomeSize);
	samclose(fp); // close the file
	*/

}






void findTranslocations(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint64_t genomeSize,
		uint16_t minimum_mapping_quality, uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs,
		float coverage, string outputFileHeader) {
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	uint32_t contigsNumber =  head.Sequences.Size();
	// now create Translocation DB
	Translocations *Trans;
	Trans = new Translocations(contigsNumber);
	Trans->initTrans(head);


	uint32_t cc = 0;
	SamSequenceDictionary sequences  = head.Sequences;
	map<unsigned int,string> position2contig;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		position2contig[cc]  = sequence->Name;
		cc++;
	}

	//Initialize bam entity
	vector<BamAlignment> currentReads;
	BamAlignment alignment;
	bamFile.GetNextAlignment(alignment);
	if(computeReadType(alignment, max_insert, outtie)!= lowQualty) {
		currentReads.push_back(alignment);
	}
	string read = alignment.Name;
	//now start to iterate over the bam file
	while ( bamFile.GetNextAlignment(alignment) ) {
		string currentRead = alignment.Name;
		if (currentRead.compare(read) != 0) { // new read under consideration
			if(currentReads.size() > 2) {
				cout << "Error, there should not be read with more than two entries\n";
				cout << currentReads.size() << " " << read << "\n";
				return;
			}
			if(currentReads.size() == 2) {
				BamAlignment read_1 = currentReads[0];
				BamAlignment read_2 = currentReads[1];
				if( (read_1.IsFirstMate() and read_2.IsSecondMate()) or (read_1.IsSecondMate() and read_2.IsFirstMate()) ) {
					// This should not be needed but I want to be sure to work always with read1 and read2
					readStatus read_1_status = computeReadType(read_1, max_insert, outtie); //there is no difference between working with the first or second reads
					if(read_1_status == pair_wrongChrs or read_1_status == pair_wrongDistance) { //read on different contigs or too far away
						// possible trans-location event found
						uint32_t startRead_1 			= read_1.Position;
						uint32_t chromosomeRead_1		= read_1.RefID;
						uint16_t qualityAligRead_1		= read_1.MapQuality;

						uint32_t startRead_2 			= read_2.Position;
						uint32_t chromosomeRead_2		= read_2.RefID;
						uint16_t qualityAligRead_2		= read_2.MapQuality;

						if(qualityAligRead_1 >= minimum_mapping_quality and qualityAligRead_2 >= minimum_mapping_quality) {
							if(chromosomeRead_1 < chromosomeRead_2) {
								Trans->insertConnection(chromosomeRead_1,startRead_1, chromosomeRead_2,startRead_2);
							} else {
								Trans->insertConnection(chromosomeRead_2,startRead_2, chromosomeRead_1,startRead_1);
							}
						}
					} else if (read_1_status == pair_wrongDistance) { // reads are too far away

					}
				} else {
					cout << "Error, both reads are either first or second read in pair, this cannot be possible probably some error in the way alignment has been performed\n";
					cout << currentReads.size() << " " << read << "\n";
					return;
				}

			}
			currentReads.clear();
			if(computeReadType(alignment, max_insert, outtie) != lowQualty) {
				currentReads.push_back(alignment);
			}
			read = currentRead;
		} else {
			if(computeReadType(alignment, max_insert, outtie) != lowQualty) {
				currentReads.push_back(alignment);
			}
		}
	}

	ofstream outputFileDescriptor;
	string fileName = outputFileHeader + "_transloacations.bed";
	outputFileDescriptor.open (fileName.c_str());
	float minCov = coverage/10;
	float maxCov = coverage*10;
	for (uint32_t i = 0; i<= Trans->chromosomesNum; i++) {
		for(uint32_t j = i +1; j<= Trans->chromosomesNum; j++) {
			Trans->findEvents(outputFileDescriptor, i,j, minimumSupportingPairs, minCov, maxCov, windowSize, windowStep);
		}
	}
	outputFileDescriptor.close();


	fileName = outputFileHeader + "_deletions.bed";
	outputFileDescriptor.open (fileName.c_str());
	for (uint32_t i = 0; i<= Trans->chromosomesNum; i++) {
		Trans->findEvents(outputFileDescriptor, i,i, minimumSupportingPairs, minCov, maxCov, windowSize, windowStep);
	}
	outputFileDescriptor.close();

	/*

	ofstream outputDeletionFile;
	outputDeletionFile.open ("output_deletions.bed");
	minimumNumberOfSupportingPairs = 2;
	minCov = 0;
	maxCov = 100;
	windowSize = 8000;
	windowStep = 1000;
	for (uint32_t i = 0; i<= Trans->chromosomesNum; i++) {
		Trans->findEvents(outputDeletionFile, i, i, minimumNumberOfSupportingPairs, minCov, maxCov, windowSize, windowStep);
	}
	*/
}



/*
 * Options.h
 *
 *  Created on: 08/mar/2011
 *      Author: Cristian Del Fabbro
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <string>
using namespace std;

#include <common.h>

#include "io/Fasta.h"
using namespace fasta;

#include "data_structures/Mask.h"


class Options {
public:
	Options() { set_defaults(); }
	Options(int argc, char *argv[]);
	virtual ~Options() { }

	bool process(int argc, char *argv[]);

	enum program_mode_t { program_unknown, program_create, program_search, program_filter_for_assembly};

	program_mode_t program_mode;

	int argc;
	char **argv;

	// input options
	string reference_file;
	bool auto_errors;
	t_errors errors_rate;
	t_errors common_errors_allowed;

	vector<string> input_files;

	Fasta::FASTQ_type fastqformat;
	bool force_fastqformat;

	bool contamination_check;
	string contamination_file;

	bool quality_check;
	bool trim;
	Mask::t_min_phred_value_CLC min_phred_value_CLC;
	Mask::t_min_phred_value_CLC min_mean_quality;
	t_pattern_length min_size;

	t_pattern_length seed_sizes;
	t_errors seed_errors;
	t_length max_gap;
	bool gap;

	string sample;
	bool bam_format;

	unsigned long int dot_every;
	bool verbose;
	bool debug;
	int threads_number;

	string query1;
	string query2;
	bool paired_ends;

	string vectors_file;

	int k;
	int blockLength;

	unsigned int delta;
	bool indels;

	// output options
	bool printAll;
	unsigned int toBePrinted;

	string output_file;



protected:
	void set_defaults();

};

#endif /* OPTIONS_H_ */

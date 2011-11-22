/*
 * read_functions.h
 *
 *  Created on: 09/mar/2011
 *      Author: cdf
 */

#ifndef READ_FUNCTIONS_H_
#define READ_FUNCTIONS_H_

#include <boost/thread/pthread/mutex.hpp>
using namespace boost;

#include "data_structures/Mask.h"
#include "data_structures/Reads.h"
#include "io/Fasta.h"
using namespace fasta;
using namespace reads;

int read_sequences(istream & input, int num_seq, Mask sequences[], Fasta::FASTQ_type format_type);
int read_sequences_no_threads(istream & input, int num_seq, Reads sequences[], Fasta::FASTQ_type format_type, long limit);

#endif /* READ_FUNCTIONS_H_ */

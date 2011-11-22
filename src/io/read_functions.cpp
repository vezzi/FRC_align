/*
 * read_functions.cpp
 *
 *  Created on: 09/mar/2011
 *      Author: cdf
 */

#include "io/read_functions.h"

mutex read_mutex;

int read_sequences(istream & input, int num_seq, Mask sequences[], Fasta::FASTQ_type format_type) {
	{
		mutex::scoped_lock lock(read_mutex);
		Fasta read;
		read.set_FASTQ_type(format_type);
		int n_seq = 0;
		while (not input.eof() and n_seq < num_seq) {
			input >> read;
			Mask & r = sequences[n_seq];
			r.set_id(read.get_id());
			r.set_sequence(read.get_sequence());
			r.set_quality(read.get_quality());
			n_seq++;
		}
		return n_seq;
	}
}

int read_sequences_no_threads(istream & input, int num_seq, Reads sequences[], Fasta::FASTQ_type format_type, long limit) {
	Fasta read;
	read.set_FASTQ_type(format_type);
	int n_seq = 0;
	while (input.tellg()< limit and (not input.eof()) and n_seq < num_seq) {
		input >> read;
		string seq = read.get_sequence();
		int numN=0;
		for(int i =0; i < seq.length(); i++) {
			if(seq.at(i) == 'N' or seq.at(i) == 'n')
				numN++;
		}
		if(numN == 0) {
			Reads & r = sequences[n_seq];
			r.initialize(read.get_sequence());
			n_seq++;
		}
	}
	return n_seq;

}

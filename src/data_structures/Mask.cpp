#include "Mask.h"

Mask::~Mask() {
}

Mask::Mask() {
	good_region_start = 0;
	good_region_stop = 0;
	globalPosition = 0;
	algn = 0;
	NM = 0;
	low_complexity = false;
	low_quality = false;
	discarded = false;
	masked = false;
	gapped = false;
	contaminated = false;
	type = unknown_alignment;
	strand = true; // needed
}


Mask::Mask(Fasta const & fasta) {
	good_region_start = 0;
	good_region_stop = 0;
	globalPosition = 0;
	algn = 0;
	NM = 0;
	masked = false;
	low_complexity = false;
	low_quality = false;
	discarded = false;
	sequence = fasta.get_sequence();
	quality = fasta.get_quality();
	gapped = false;
	contaminated = false;
	type = unknown_alignment;
	strand = true; // needed
}

ostream& operator<<(ostream& channel, const Mask& mask) {
	if (mask.masked) {
		t_pattern_length l = mask.sequence.length();
		for (t_pattern_length i=1; i <= l; i++)
			if ((i < mask.good_region_start) or (i > mask.good_region_stop))
				channel << (char)(mask.sequence[i-1] - 'A' + 'a');
			else
				channel << mask.sequence[i-1];
	} else
		channel << mask.sequence;
	return channel;
}

string Mask::to_string() const {
	stringstream ss;
	ss << *this;
	return ss.str();
}

string Mask::get_good_sequence() const {
	if (discarded)
		return string();
	if (masked)
		return sequence.substr(good_region_start-1,good_region_stop-good_region_start+1);
	else
		return sequence;
}

string Mask::get_good_quality() const {
	if (discarded)
		return string();
	if (masked)
		return quality.substr(good_region_start-1,good_region_stop-good_region_start+1);
	else
		return quality;
}


string Mask::get_masked() const {
	if (discarded) {
		string temp = sequence;
		for (t_pattern_length i = 0; i < sequence.length(); i++)
			temp[i] = tolower(temp[i]);
		return temp;
	} else if (masked) {
		string temp = sequence;
		for (t_pattern_length i = 1; i < good_region_start; i++)
			temp[i-1] = tolower(temp[i-1]);
		for (t_pattern_length i = good_region_stop+1; i <= temp.size(); i++)
			temp[i-1] = tolower(temp[i-1]);
		return temp;
	} else
		return sequence;
}

int Mask::get_five_prime() {
	if (strand) {
		if (masked)
			return position - good_region_start + 1;
		else
			return position;
	} else {
		if (masked)
			return position - good_region_stop - 1;
		else
			return position + sequence.size() - 1;
	}
}

string Mask::get_MD(const char * reference) {
	stringstream MD;
	int e = 0;

	string forORrev;
	if (strand)
		forORrev = get_good_sequence();
	else
		forORrev = reverse_complement_standalone_str(get_good_sequence());

	const char * fOr = forORrev.c_str();
	int l = get_good_length();
	reference += globalPosition;
	for (int i = 0; i < l ; i++) {
		if (*reference == *fOr)
			e++;
		else {
			if (e > 0)
				MD << e;
			MD << *fOr;
			e = 0;
		}
		reference++;
		fOr++;
	}

	if (e != 0)
		MD << e;
	return MD.str();
}

Mask::t_CIGAR Mask::get_CIGAR() {
	stringstream CIGAR;
	unsigned int operations = 0;
	unsigned int i;
	if (strand) {
		if (masked and good_region_start > 1) {
			//CIGAR << (good_region_start-1) << "S";
			i = (good_region_start-1)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
		//CIGAR << get_good_length() << "M";
		i = get_good_length()*16 + BAM_CMATCH;
		CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
		operations++;
		if (masked and good_region_stop < sequence.size()) {
			//CIGAR << (sequence.size() - good_region_stop) << "S";
			i = (sequence.size() - good_region_stop)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
	} else {
		if (masked and good_region_stop < sequence.size()) {
			//CIGAR << (sequence.size() - good_region_stop) << "S";
			i = (sequence.size() - good_region_stop)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
		//CIGAR << get_good_length() << "M";
		i = get_good_length() * 16 + BAM_CMATCH;
		CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
		operations++;
		if (masked and good_region_start > 1) {
			//CIGAR << (good_region_start-1) << "S";
			i = (good_region_start-1)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
	}
	return t_CIGAR(CIGAR.str(),operations);
}


Mask::t_CIGAR Mask::get_CIGAR_SW() {
	stringstream CIGAR;
	unsigned int operations = 0;
	unsigned int i;
	if (strand) {
		if (masked and good_region_start > 1) {
			//CIGAR << (good_region_start-1) << "S";
			i = (good_region_start-1)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
		//I have to perform SW in order to compute best score and backtrack to create CIGAR


		for(int k=0; k< cigarVector.size(); k++) {
			switch( cigarVector.at(k).first) {
			case 'M':
				i = cigarVector.at(k).second*16 + BAM_CMATCH ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			case 'D':
				i = cigarVector.at(k).second*16 + BAM_CDEL ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			case 'I':
				i = cigarVector.at(k).second*16 + BAM_CINS ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			}
			operations++;
		}
		if (masked and good_region_stop < sequence.size()) {
			//CIGAR << (sequence.size() - good_region_stop) << "S";
			i = (sequence.size() - good_region_stop)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
	} else {
		if (masked and good_region_stop < sequence.size()) {
			//CIGAR << (sequence.size() - good_region_stop) << "S";
			i = (sequence.size() - good_region_stop)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}

		for(int k=0; k< cigarVector.size(); k++) {
			switch( cigarVector.at(k).first) {
			case 'M':
				i = cigarVector.at(k).second*16 + BAM_CMATCH ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			case 'D':
				i = cigarVector.at(k).second*16 + BAM_CDEL ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			case 'I':
				i = cigarVector.at(k).second*16 + BAM_CINS ;
				CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
				break;
			}
			operations++;
		}
		if (masked and good_region_start > 1) {
			//CIGAR << (good_region_start-1) << "S";
			i = (good_region_start-1)*16 + BAM_CSOFT_CLIP;
			CIGAR << (char)(i%256) << (char)(i/256) << '\0' << '\0';
			operations++;
		}
	}
	return t_CIGAR(CIGAR.str(),operations);
}

void Mask::extract_fasta_masked(Fasta & output) const {
	output.set_sequence(get_good_sequence().c_str());
	output.set_quality(get_good_quality().c_str());
}


void Mask::set_fasta_masked(Fasta & output) const {
	output.set_id(get_id().c_str());
	output.set_sequence(get_good_sequence().c_str());
	output.set_quality(get_good_quality().c_str());
}

string Mask::print_information() const {
	stringstream ss;
	if (masked)
		ss << "=== sequence is masked === " << endl;
	else
		ss << "=== sequence is NOT masked ===" << endl;
	ss << "REAL_SEQUENCE:\t" << sequence << endl;
	ss << "REAL_QUALITY:\t" << quality << endl;
	ss << "GOOD SEQUENCE:\t" << get_good_sequence() << endl;
	ss << "GOOD QUALITY:\t" << get_good_quality() << endl;
	if (masked)
		ss << "Sequence is MASKED" << endl;
	if (low_complexity)
		ss << "Sequence is LOW COMPLEXITY" << endl;
	if (low_quality)
		ss << "Sequence is LOW QUALITY" << endl;
	if (trimmed)
		ss << "Sequence is TRIMMED" << endl;
	if (discarded)
		ss << "Sequence is DISCARDED" << endl;
	if (contaminated)
		ss << "Sequence is CONTAMINATED" << endl;
	ss << "5':\t " << good_region_start << endl;
	ss << "3':\t " << good_region_stop << endl;

	return ss.str();
}

char * Mask::quality_conversion() {
	const char * original = quality.c_str();
	size_t size = quality.size();
	char * quality_converted = new char[size+1];
	quality_converted[size]= '\0';
	if (strand) {
		for (size_t i = 0; i < size; i++)
			quality_converted[i] = original[i] - 33;
	} else {
		for (size_t i = 0; i < size; i++)
			quality_converted[size-i-1] = original[i] - 33;
	}
	return quality_converted;
}

void Mask::quality_trimming_CLC(t_min_phred_value_CLC min_phred_value_CLC, t_min_phred_value_CLC min_mean_quality, t_min_phred_value_CLC min_size) {
	const char *q = quality.c_str();
	size_t l = quality.size();

	// finding left start position
	size_t left = 0;
	bool not_found = true;
	while (not_found and (left < l)) {
		if ((q[left] - 33) >= min_phred_value_CLC)
			not_found = false;
		else
			left++;
	}

	if (not_found) {
		discarded = true;
		low_quality = true;
	} else {
		// finding right end position
		int a = q[left] - 33 - min_phred_value_CLC;
		int max = a;
		size_t right = left;
		unsigned int mean_sum = a;
		unsigned int mean_sum_temp = a;
		for (size_t i = left + 1; i < l; i++) {
			a += q[i] - 33 - min_phred_value_CLC;
			mean_sum_temp += q[i] - 33;
			if (a < 0)
				a = 0;
			if (max <= a) {
				max = a;
				right = i;
				mean_sum = mean_sum_temp;
			}
		}

		if ((right - left + 1 < min_size) or (mean_sum / (right - left + 1) < min_mean_quality)) {
			discarded = true;
			low_quality = true;
		} else {
			// good region is [left,right]
			masked = true;
			good_region_start = left+1;
			good_region_stop = right+1;
		}
	}

}

void Mask::compact_DNA(string & compact) {
	const char * dna = sequence.c_str();
	size_t l = sequence.size();
	compact.clear();
	compact.reserve((size_t)((l+1)/2)+1);
	char first;
	char second;
	if (strand)
		for (size_t i = 0; i < l; i+=2) {
			switch (dna[i]) {
			case 'A' :
			case 'a' : first = 1; break;
			case 'C' :
			case 'c' : first = 2; break;
			case 'G' :
			case 'g' : first = 4; break;
			case 'T' :
			case 't' : first = 8; break;
			default  : first = 15;
			}
			switch ((i+1)<l ? dna[i+1] : 'N') {
			case 'A' :
			case 'a' : second = 1; break;
			case 'C' :
			case 'c' : second = 2; break;
			case 'G' :
			case 'g' : second = 4; break;
			case 'T' :
			case 't' : second = 8; break;
			default  : second = 15;
			}
			compact.push_back((char)(first*16+second));
		}
	else
		for (size_t i = 0; i < l; i+=2) {
			switch (dna[l-i-1]) {
			case 'A' :
			case 'a' : first = 8; break;
			case 'C' :
			case 'c' : first = 4; break;
			case 'G' :
			case 'g' : first = 2; break;
			case 'T' :
			case 't' : first = 1; break;
			default  : first = 15;
			}
			switch ((i+1)<l ? dna[l-i-2] : 'N') {
			case 'A' :
			case 'a' : second = 8; break;
			case 'C' :
			case 'c' : second = 4; break;
			case 'G' :
			case 'g' : second = 2; break;
			case 'T' :
			case 't' : second = 1; break;
			default  : second = 15;
			}
			compact.push_back((char)(first*16+second));
		}
}

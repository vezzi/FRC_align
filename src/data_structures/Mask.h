#ifndef MASK_H_
#define MASK_H_

#include <sstream>

#include "samtools/bam.h"

#include "io/Fasta.h"
using namespace fasta;

#include "common.h"
using namespace std;

class Mask {
public:

	typedef pair<string,unsigned int> t_CIGAR;
	typedef unsigned short int t_min_phred_value_CLC;

	virtual ~Mask();
	Mask();
	Mask(Fasta const & fasta);

	friend ostream& operator<<(ostream& channel, const Mask& mask);
	string to_string() const;
	string get_good_sequence() const;
	INLINE const string & get_original_sequence() const { return sequence ; }
	string get_masked() const;
	INLINE const string & get_original_quality() const { return quality; }
	string get_good_quality() const;
//	void merge(Mask &mask);
	string print_information() const;
	int get_five_prime();
	string get_MD(const char * reference);
	t_CIGAR get_CIGAR();
	t_CIGAR get_CIGAR_SW();
	void quality_trimming_CLC(t_min_phred_value_CLC min_phred_value_CLC, t_min_phred_value_CLC min_mean_quality, t_min_phred_value_CLC min_size);
	char * quality_conversion();
	void compact_DNA(string & compact);
	
	void extract_fasta_masked(Fasta & output) const;
	void set_fasta_masked(Fasta & output) const;

	void low_quality_discard() { low_quality = true; discarded = true; }

	INLINE bool status_masked() { return masked; }
	INLINE bool status_low_complexity() { return low_complexity; }
	INLINE bool status_low_quality() { return low_quality; }
	INLINE bool status_trimmed() { return trimmed; }
	INLINE bool status_discarded() { return discarded; }

	INLINE void set_sequence(string seq) { sequence = seq;  }
	INLINE void set_quality(string qual) { quality = qual; }
	INLINE int get_sequenceLength() {return sequence.length();}
	INLINE int get_qualityLength() {return quality.length();}


	INLINE void set_id(string id1) {  id = id1; }
	INLINE string get_id() const { return id; }

	//reads object variables for print the solution

	INLINE void set_type(t_alignment T) { type = T;}
	INLINE t_alignment get_type() { return type;}

	INLINE void set_algn(int A) {algn = A;}
	INLINE int get_algn() { return algn;}

	INLINE void set_contig(size_t C) { contig = C;}
	INLINE size_t get_contig() { return contig;}

	INLINE void set_strand(bool S) { strand = S;}
	INLINE bool get_strand() { return strand;}

	INLINE void set_NM(int nm) { NM = nm;}
	INLINE int get_NM() { return NM;}

	INLINE void set_position(unsigned long int pos) { position = pos;}
	INLINE unsigned long int get_position() { return position;}

	INLINE int get_good_length() { return masked ? good_region_stop - good_region_start + 1 : sequence.length();}

	INLINE void set_globalPosition(unsigned long int p) {globalPosition = p;}
	INLINE unsigned long int get_globalPosition() { return globalPosition;}

	INLINE t_pattern_length get_good_region_start() { return good_region_start; }
	INLINE t_pattern_length get_good_region_stop()  { return good_region_stop; }

//protected:

	// da trasferire come input e output
	string id;
	string sequence;
	string quality;

	// da trasferire come output
	int algn;  //in quante posizioni ho allineato
	bool strand;  //true = +, false = -
	int NM;
	unsigned long int globalPosition;

	bool masked;
	bool low_complexity;
	bool low_quality;
	bool trimmed;
	bool discarded;
	bool contaminated;

	t_pattern_length good_region_start;
	t_pattern_length good_region_stop;

	bool gapped;
	int contig;
	int NM_gap;
	unsigned long int globalPosition_gap;
	int length1_gap;
	int length2_gap;

	//Delta
	vector<pair<int,int> > DELTA;
	unsigned int XD;

	//indels
	bool indels;
	vector<pair<char, unsigned int> > cigarVector;

	//SAM BAM print options
	bool primary; // primary alignment
	unsigned int IH; // number of printed solutions
	unsigned int HI; // solution printed


	// Non serve trasferire

	t_alignment type;  //NO
	unsigned long int position; // NO
	unsigned long int position_gap; // NO

	friend class Auto_Trim;
	friend class Low_Complexity;
	friend class PolyA;
	friend class Quality_Check;
	friend class Trim;
	friend class Use_Bases;
	friend class Vector;

};


#endif /*MASK_H_*/

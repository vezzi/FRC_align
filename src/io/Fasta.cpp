#include "Fasta.h"

#define BUFFER_LEN 10240

namespace fasta
{

Fasta::~Fasta() {
}

/**
 * Simple constructor with empty data.
 * @ param cs is colorspace? (optional, default = false)
 */
Fasta::Fasta(bool cs) {
	colorspace = cs;
	have_quality = false;
	columns = DEFAULT_FASTA_COLUMNS;
	mask_lowercase = false;
	type = standard;
}

/**
 * Constructor from an id and a string.
 * @param i the id
 * @param s the sequence
 * @param cs is colorspace? (optional, default = false)
 */
Fasta::Fasta(const char *i, const char *s, bool cs) {
	id = string(i);
	sequence = string(s);
	colorspace = cs;
	have_quality = false;
	columns = DEFAULT_FASTA_COLUMNS;
	mask_lowercase = false;
	type = standard;
}

/**
 * Constructor from an id and a string.
 * @param i the id
 * @param s the sequence
 * @param cs is colorspace? (optional, default = false)
 */
Fasta::Fasta(const string &i, const string &s, bool cs) {
	id = string(i.c_str());
	sequence = string(s.c_str());
	colorspace = cs;
	have_quality = false;
	columns = DEFAULT_FASTA_COLUMNS;
	mask_lowercase = false;
	type = standard;
}

/**
 * Constructor from an id and a string.
 * @param i the id
 * @param s the sequence
 * @param c the comment
 * @param cs is colorspace? (optional, default = false)
 */
Fasta::Fasta(const char *i, const char *s, const char *c, bool cs) {
	id = string(i);
	sequence = string(s);
	comment = string(c);
	colorspace = cs;
	have_quality = false;
	columns = DEFAULT_FASTA_COLUMNS;
	mask_lowercase = false;
	type = standard;
}

/**
 * Constructor from an id and a string.
 * @param i the id
 * @param s the sequence
 * @param c the comment
 * @param cs is colorspace? (optional, default = false)
 */
Fasta::Fasta(const string &i, const string &s, const string &c, bool cs) {
	id = string(i.c_str());
	sequence = string(s.c_str());
	comment = string(c.c_str());
	colorspace = cs;
	have_quality = false;
	columns = DEFAULT_FASTA_COLUMNS;
	mask_lowercase = false;
	type = standard;
}

/*
 * Return a short rappresentation of the fasta.
 * @param n how many characters shows, if 0 shows entire sequence
 */
string Fasta::description(size_t n) const {
	stringstream s;
	s << "ID : " << id << endl;
	s << "STR: ";
	if (n == 0 or n >= sequence.length())
		s << sequence;
	else {
		s << sequence.substr(0,n-1);
		s << "[...]";
	}
	s << endl << "LEN: " << sequence.length() << endl;;
	return s.str();
}

istream & Fasta::read_from_fasta(istream & buffer) throw (Incorrect_Format) {
	have_quality = false;
	char b[BUFFER_LEN+1];
	b[BUFFER_LEN] = '\0';
	stringstream buf;
	char c = buffer.peek();
	while (!buffer.eof() and (c == ' ' or c == '\n')) {
		buffer.ignore(1);
		if (!buffer.eof()) c = buffer.peek();
	}
	if (!buffer.eof() and c != '>') {
		stringstream ss;
		ss << "next character is " << c;
		throw Incorrect_Format(ss.str());
	}
	buffer.getline(b, BUFFER_LEN);
	id = string(b+1);
	long int pos = id.find(' ');
	if (pos > 0) {
		comment = id.substr(pos+1,id.length());
		id = id.substr(0,pos);
	}
	string temp;
	char char_temp;
	while (!buffer.eof() and (buffer.peek() != '>')) {
		buffer >> temp;
		for (string::iterator iter = temp.begin(); iter != temp.end(); iter++) {
			if (colorspace and (*iter == '.'))
				buf << text_delimitator;
			else if ((*iter >= 'a') and (*iter <= 'z')) {
				if (mask_lowercase)
					buf << masked_base;
				else {
					char_temp = *iter;
					char_temp = char_temp - 'a' + 'A';
					buf << char_temp;
				}
			} else if ((*iter >= 'A') and (*iter <= 'Z'))
				buf << *iter;
			else
				throw Incorrect_Format("Incorrect FASTA sequence");
		}
		c = buffer.peek();
		while (!buffer.eof() and (c == ' ' or c == '\n')) {
			buffer.ignore(1);
			if (!buffer.eof()) c = buffer.peek();
		}
	}
	sequence = string(buf.str());
	quality = string(sequence.size(),(char)(DEFAULT_QUALITY_VALUE+33));
	return buffer;
}

istream & Fasta::read_from_fastq(istream & buffer) throw (Incorrect_Format) {
	have_quality = true;
	char b[BUFFER_LEN+1];
	b[BUFFER_LEN] = '\0';
	char c = buffer.peek();
	while (!buffer.eof() and (c == ' ' or c == '\n')) {
		buffer.ignore(1);
		if (!buffer.eof()) c = buffer.peek();
	}
	if (!buffer.eof() and c != '@') {
		stringstream ss;
		ss << "next character is " << c << " instead of '@'";
		throw Incorrect_Format(ss.str());
	}

	buffer.ignore(1);
	buffer.getline(b, BUFFER_LEN);
	id = string(b);
	long int pos = id.find(' ');
	if (pos > 0) {
		comment = id.substr(pos+1,id.length());
		id = id.substr(0,pos);
	}
	if (buffer.eof())
		throw Incorrect_Format("file ended prematurely");

	c = buffer.peek();
	stringstream buf;
	char char_temp;
	while (!buffer.eof() and c != '+') {
		string temp;
		buffer >> temp;
		for (string::iterator iter = temp.begin(); iter != temp.end(); iter++) {
			if (colorspace and (*iter == '.'))
				buf << text_delimitator;
			else if ((*iter >= 'a') and (*iter <= 'z')) {
				if (mask_lowercase)
					buf << masked_base;
				else {
					char_temp = *iter;
					char_temp = char_temp - 'a' + 'A';
					buf << char_temp;
				}
			} else if ((*iter >= 'A') and (*iter <= 'Z'))
				buf << *iter;
			else
				throw Incorrect_Format("Incorrect FASTQ sequence");
		}
		if (buffer.eof())
			throw Incorrect_Format("file ended prematurely");

		c = buffer.peek();
		while (!buffer.eof() and (c == '\n')) {
			buffer.ignore(1);
			c = buffer.peek();
		}
	}
	sequence = buf.str();
	if (buffer.eof())
		throw Incorrect_Format("file ended prematurely");
	if (!buffer.eof() and c != '+') {
		stringstream ss;
		ss << "next character is " << c << " instead of '+'";
		throw Incorrect_Format(ss.str());
	}

	if (buffer.eof())
		throw Incorrect_Format("file ended prematurely");
	buffer.getline(b, BUFFER_LEN); //only for remove line starts with '+'

	if (buffer.eof())
		throw Incorrect_Format("file ended prematurely");
	buffer >> quality; // this is real quality!
	if (type == illumina)
		for (size_t i = 0; i < quality.size(); i++)
			quality[i] -= 31;

	while (not buffer.eof() and (buffer.peek() == '\n'))
		buffer.ignore(1);

	if (not buffer.eof()) {
		c = buffer.peek();
		while (not buffer.eof() and c != '@') {
			string temp;
			buffer >> temp;
			quality.append(temp);

			if (not buffer.eof()) {
				c = buffer.peek();
				while (not buffer.eof() and (c == '\n')) {
					buffer.ignore(1);
					c = buffer.peek();
				}
			}
		}
	}

	if (sequence.size() != quality.size()) {
		stringstream ss;
		ss << "sequence and quality have different number of characters";
		throw Incorrect_Format(ss.str());
	}

	for (string::iterator iter = sequence.begin(); iter != sequence.end(); iter++) {
		if ((*iter >= 'a') and (*iter <= 'z')) {
			if (mask_lowercase)
				*iter = masked_base;
			else
				*iter = *iter - 'a' + 'A';
		}
		if (colorspace and (*iter == '.'))
			*iter = text_delimitator;
	}
	while (!buffer.eof() and (buffer.peek() != '@')) {
		buffer.ignore(1);
	}

	return buffer;
}


/**
 * Reading a fasta sequence from a stream.
 * @param buffer the streamer
 * @param fasta the object where store the read sequence
 */
istream & operator>> (istream & buffer, Fasta & fasta) throw (Incorrect_Format) {
	if (not buffer.eof()) {
		wchar_t c;

		while (not buffer.eof()) {
			c = buffer.peek();
			switch (c) {
			case '>' :
				fasta.read_from_fasta(buffer);
				return buffer;
			case '@' :
				fasta.read_from_fastq(buffer);
				return buffer;
			default:
				buffer.ignore(1);
			}
		}
	}

	return buffer;
}

ostream & operator<< (ostream & buffer, const Fasta & fasta) {
	if (fasta.have_quality) {
		buffer << '@' << fasta.id;
		if (fasta.comment.length() > 0)
			buffer << " " << fasta.comment;
		buffer << endl;
		/*
		for (size_t i = 0; i < fasta.sequence.length(); i += fasta.columns) {
			buffer << fasta.sequence.substr(i,fasta.columns) << endl;
		}
		*/
		buffer << toupper(fasta.sequence) << endl;
		buffer << '+' << endl;
		/*
		for (size_t i = 0; i < fasta.quality.length(); i += fasta.columns) {
			buffer << fasta.quality.substr(i,fasta.columns) << endl;
		}
		*/
		buffer << fasta.quality << endl;
	} else {
		buffer << '>' << fasta.id;
		if (fasta.comment.length() > 0)
			buffer << " " << fasta.comment;
		buffer << endl;
		for (size_t i = 0; i < fasta.sequence.length(); i += fasta.columns) {
			buffer << toupper(fasta.sequence.substr(i,fasta.columns)) << endl;
		}
	}
	return buffer;
}


string Fasta::reverse_complement() const {
	return reverse_complement_standalone_str(sequence.c_str());
}

string Fasta::reverse_complement(size_t start, size_t stop) const throw (Data_Exception) {
	if (start >=sequence.length())
		throw Data_Exception(0,sequence.length()-1,start,"Wrong start position");
	if (stop >= sequence.length())
		throw Data_Exception(0,sequence.length()-1,stop,"Wrong stop position");
	if (start > stop)
		throw Data_Exception(0,stop,start,"start is greater than stop");
	string reverse;
	if (colorspace)
		for (size_t i = stop-1; i >= start ; i--)
			reverse.push_back(sequence[i]);
	else
		for (size_t i = stop-1; i >= start ; i--)
			reverse.push_back(reverse_complement_standalone_char(sequence[i]));
	return reverse;
}

/**
 * Return the sequence in the reverse direction (useful for colorspace)
 * @return the reversed sequence
 */
string Fasta::reverse() const {
	return reverse_standalone_str(sequence.c_str());
}

string Fasta::reverse(size_t start, size_t stop) const throw (Data_Exception) {
	if (start >=sequence.length())
		throw Data_Exception(0,sequence.length()-1,start,"Wrong start position");
	if (stop >= sequence.length())
		throw Data_Exception(0,sequence.length()-1,stop,"Wrong stop position");
	if (start > stop)
		throw Data_Exception(0,stop,start,"start is greater than stop");
	string reverse;
	for (size_t i = stop-1; i >= start ; i--)
		reverse.push_back(sequence[i]);
	return reverse;
}

/*
void Fasta::masking(Masked_Sequence& ms, char masking_character) {
	ms.mask_sequence(sequence,masking_character);
}
*/

t_quality Fasta::get_quality(size_t position) const throw (Data_Exception) {
	if (position >= quality.length())
		throw Data_Exception(0,quality.length()-1,position,"Wrong position for quality");
	return quality[position] - 33;
}

Fasta::t_quality_vector Fasta::get_quality_vector() const {
	t_quality_vector v;
	for (size_t position = 0; position < quality.length(); position++)
		v.push_back(quality[position] - 33);
	return v;
}

void Fasta::loweringNs() {
	for (size_t i = 0; i < sequence.size(); i++)
		if (sequence[i] == 'N')
			sequence[i] = 'n';
}

void Fasta::convertToColorspace() {
	const size_t length = sequence.size();
	const char * old = sequence.c_str();
	char new_sequence[length];
	new_sequence[length-1] = '\0';

	for (size_t i = 0; i < (length-1); i++) {
		unsigned short int first;
		unsigned short int second;
		switch (toupper(old[i])) {
		case 'A' : first = 0; break;
		case 'C' : first = 1; break;
		case 'G' : first = 2; break;
		case 'T' : first = 3; break;
		default:   first = 4; break;
		}
		switch (toupper(old[i+1])) {
		case 'A' : second = 0; break;
		case 'C' : second = 1; break;
		case 'G' : second = 2; break;
		case 'T' : second = 3; break;
		default:   second = 4; break;
		}
		if (first == 4 or second == 4)
			new_sequence[i] = '4';
		else
			new_sequence[i] = colorspace_conversion_table[first][second];
	}
}

Fasta::FASTQ_type Fasta::check_FASTQ_type_file(const string &file) {
	Auto_Unzip is(file);
	istream & filtered = is.filtered();
	string line;

	if (not filtered.eof() and filtered.peek() == '>')
		return standard; // it is not a FASTQ file
	while (not filtered.eof()) {
		getline(filtered,line); // skip id
		if (filtered.eof())
			return unknown;
		getline(filtered,line); // skip seq
		if (filtered.eof())
			return unknown;
		getline(filtered,line); // skip delimiter
		if (filtered.eof())
			return unknown;
		getline(filtered,line); // get quality!
		size_t l = line.size();
		const char *str = line.c_str();
		for (size_t i = 0; i < l; i++)
			if (str[i] < 59)
				return standard;
			else if (str[i] > 73)
				return illumina;
	}
	return unknown;
}

}

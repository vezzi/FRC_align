/*
 * Auto_Unzip.cpp
 *
 *  Created on: 3-apr-2009
 *      Author: cdf
 */

#include "Auto_Unzip.h"

namespace useful {

Auto_Unzip::Auto_Unzip() {
	input_stream = NULL;
	filtered_stream = NULL;
}

Auto_Unzip::Auto_Unzip(const string & filename) throw (File_Not_Found) {
	input_stream = NULL;
	filtered_stream = NULL;
	open(filename);
}

Auto_Unzip::Auto_Unzip(const char * filename) throw (File_Not_Found) {
	input_stream = NULL;
	filtered_stream = NULL;
	open(filename);
}

Auto_Unzip::~Auto_Unzip() {
	if (filtered_stream != NULL)
		delete filtered_stream;
	if (input_stream != NULL)
		delete input_stream;
}

void Auto_Unzip::open(const char * filename) throw (File_Not_Found) {
	string temp(filename);
	open(temp);
}

void Auto_Unzip::open(const string & filename) throw (File_Not_Found) {
	if (filtered_stream != NULL)
		delete filtered_stream;
	if (input_stream != NULL)
		delete input_stream;

	if (filename.compare("-") == 0) {
		streambuf * psbuf  = cin.rdbuf();
		filtered_stream = new iostreams::filtering_istream();
		filtered_stream->rdbuf(psbuf);
		input_stream = NULL;
	} else {
		input_stream = new ifstream(filename.c_str());
		if (not (*input_stream))
			throw File_Not_Found(filename);

		filtered_stream = new iostreams::filtering_istream();
		size_t pos = filename.rfind(".gz");
		if ((pos != string::npos) and (pos == (filename.size()-3)))
			filtered_stream->push(iostreams::gzip_decompressor());
		else {
			pos = filename.rfind(".bz2");
			if ((pos != string::npos) and (pos == (filename.size()-4)))
				filtered_stream->push(iostreams::bzip2_decompressor());
		}
		filtered_stream->push(*input_stream);
	}
}

streampos Auto_Unzip::tellg() {
	if (input_stream != NULL)
		return input_stream->tellg();
	else
		return -1;
}

void Auto_Unzip::seekg(streampos pos) {
	cout << "inside seekg " << pos << "\n";
	input_stream->seekg  (0, ios::end);;
	cout << input_stream->tellg() << "\n";
}

streampos Auto_Unzip::length() {
	if (input_stream != NULL) {
		streampos pos = input_stream->tellg();
		input_stream->seekg (0, ios::end);
		streampos length = input_stream->tellg();
		input_stream->seekg (pos);
		return length;
	} else
		return 0;
}

}

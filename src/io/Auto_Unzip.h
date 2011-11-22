/*
 * Auto_Unzip.h
 *
 *  Created on: 3-apr-2009
 *      Author: cdf
 */

#ifndef AUTO_UNZIP_H_
#define AUTO_UNZIP_H_

#include <fstream>
#include <iostream>
using namespace std;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
using namespace boost;

#include "errors/File_Not_Found.h"
using namespace errors;

namespace useful {

class Auto_Unzip {
public:
	Auto_Unzip();
	Auto_Unzip(const char * filename) throw (File_Not_Found);
	Auto_Unzip(const string & filename) throw (File_Not_Found);
	virtual ~Auto_Unzip();

	void open(const char * filename) throw (File_Not_Found);
	void open(const string & filename) throw (File_Not_Found) ;
	istream & filtered() {
		if (filtered_stream == NULL)
			return *input_stream;
		else
			return *filtered_stream;
	}
	ifstream & file() { return *input_stream; }
	streampos tellg();
	streampos length();
	void seekg(streampos pos);

protected:
	ifstream * input_stream;
	iostreams::filtering_istream * filtered_stream;

};

}

#endif /* AUTO_UNZIP_H_ */

#ifndef READS_H_
#define READS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
using namespace std;



#include "common.h"




namespace reads
{



class Reads
{
public:
	virtual ~Reads();
	Reads() {};
	Reads(string s);
	Reads(const Reads &r); // copy constructor
	Reads & operator=(const Reads & r);


	void initialize(string s);
	unsigned short length() const {	return l; } //* Return the length of the sequence.
	unsigned int * get_read() { return readP; } //* Return the read
	string  toString();


protected:
	unsigned short int l;
	unsigned int *readP; // read as given in input


};



}
#endif /*READS_H_*/

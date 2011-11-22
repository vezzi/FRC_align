#include "Reads.h"



namespace reads
{



Reads::~Reads() {

	delete [] readP;
	//delete readM;
}



Reads::Reads(const string s) {
	l = s.length();

	unsigned short int bufferREAD_length = ceil((s.length())/(float)16);
	readP = new unsigned int[bufferREAD_length];
	memset(readP , 0 , bufferREAD_length);
	int actualByte = 0;
	int actualBit = 0;

	for (unsigned int i = 0; i < s.length() ; i++) {
		switch (s.at(i)) {
		case 'A' : case 'a' : readP[actualByte] &= ~(1 << (31-actualBit)); readP[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		case 'C' : case 'c' : readP[actualByte] &= ~(1 << (31-actualBit)); readP[actualByte] |= (1 << (31-(actualBit+1)));  break;
		case 'G' : case 'g' : readP[actualByte] |= (1 << (31-actualBit)); readP[actualByte] &= ~(1 << (31-(actualBit+1)));   break;
		case 'T' : case 't' : readP[actualByte] |= (1 << (31-actualBit)); readP[actualByte] |= (1 << (31-(actualBit+1)));  break;
		}
		actualBit+=2;
		if(actualBit == 32) {
			actualByte++;
			actualBit=0;
		}
	}

}


void Reads::initialize(const string s) {
	l = s.length();

	unsigned short int bufferREAD_length = ceil((s.length())/(float)16);
	readP = new unsigned int[bufferREAD_length];
	memset(readP , 0 , bufferREAD_length);
	int actualByte = 0;
	int actualBit = 0;

	for (unsigned int i = 0; i < s.length() ; i++) {
		switch (s.at(i)) {
		case 'A' : case 'a' : readP[actualByte] &= ~(1 << (31-actualBit)); readP[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		case 'C' : case 'c' : readP[actualByte] &= ~(1 << (31-actualBit)); readP[actualByte] |= (1 << (31-(actualBit+1)));  break;
		case 'G' : case 'g' : readP[actualByte] |= (1 << (31-actualBit)); readP[actualByte] &= ~(1 << (31-(actualBit+1)));   break;
		case 'T' : case 't' : readP[actualByte] |= (1 << (31-actualBit)); readP[actualByte] |= (1 << (31-(actualBit+1)));  break;
		}
		actualBit+=2;
		if(actualBit == 32) {
			actualByte++;
			actualBit=0;
		}
	}

}



Reads::Reads(const Reads &r) {
	readP = r.readP;
	l= r.l;
}

Reads & Reads::operator=(const Reads & r) {
	if (this != &r) { // protect against invalid self-assignment
		this->l = r.l;
		this->readP = r.readP;
	}
	return *this;
}



string  Reads::toString() {
	unsigned int Mask = 0;
	unsigned int  block= 0;
	Mask |= 1 << 31;
	Mask |= 1 << 30;

	string out="";

	unsigned int R = readP[0];

	for(unsigned int i =0 ; i < l ; i++) {
		if(i%16 == 0 && i>0) {
			block++;
			R = readP[block];
		}
		unsigned int c = (R & Mask) >> 30;
		R = R << 2;
		switch(c) {
			case 0: out.append("A"); break;
			case 1: out.append("C"); break;
			case 2: out.append("G"); break;
			case 3: out.append("T"); break;
			default: cout << "strange "<< cout << "\n";
		}
	}

	return out;
}

}


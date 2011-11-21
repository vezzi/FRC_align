/* The MIT License                      
   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOS
*/

/*Contact <paul.igor.costea@scilifelab.se>*/

#include <string>
#include "util.h"

namespace pol_util
{

  ///> Defines
  #define MAX_LINE 1000

  /**
   * @description Fastq entry reading/writing + manipulation class
   */
  class FastqEntry{
  public:
    /**
     * @brief Entry reader and generator. The only way to get a new instance of a FastqEntry.
     * @brief Call this function to read another entry into memory and get a pointer to it.
     * @return Pointer to newly created FastqEntry
     */
    static FastqEntry* readEntry(FILE *file) {
      char line[4][MAX_LINE];
      int x;
      for (x=0; x<4; ++x) {//Fastq entry has 4 lines.
	if (NULL == fgets(line[x],MAX_LINE,file)) {
	  //EOF reached
	  return NULL;
	}
	//Remove end line character   
	line[x][strlen(line[x])-1] = '\0';

      }
      //DO some checking
      if (line[0][0] != '@') {
	fprintf(stderr,"[pol_util]:readEntry -> File seems to not contain fastq entries");
	return NULL;
      }
      //Entry read successfully. return FastqEntry
      return new FastqEntry(line);

    };

    /**
     * @brief Destructor
     */
    ~FastqEntry(){
      delete[] m_strName;
      delete[] m_strSeq;
      delete[] m_strPlus;
      delete[] m_strQual;
    };

    /**
     * @brief Return entry formated as a 4 line string.
     */
    std::string toString() {
      std::string entry = "";
      entry += m_strName;
      entry += '\n';
      entry += m_strSeq;
      entry += '\n';
      entry += m_strPlus;
      entry += '\n';
      entry += m_strQual;
      entry += '\n';
      return entry;
    };

    /**
     * @brief toString wrapper for direct file writing.
     */
    void write(FILE* file) {
      fprintf(file,"%s",this->toString().c_str());
    };

    /**
     * @brief Trim read the specified quality and minimum lenght
     * @return True if read is still usable
     */
    bool trim(int qual, int min_length)
    {
      int trim_qual = qual;
      int s = 0, l, max = 0, max_l = strlen(m_strQual) - 1;

      for (l = strlen(m_strQual) - 1; l >= 1; --l) {
	s += trim_qual - (m_strQual[l] - 64);
	if (s < 0) break;
	if (s > max) {
	  max = s; max_l = l;
	}
      }

      //Ignore this if len < cutoff
      if (max_l < min_length) {
	return false;
      }
      //Trim seq and quality strings
      m_strQual[max_l] = '\0';
      //qual[max_l+1] = '\0';
      m_strSeq[max_l] = '\0';
      //seq[max_l+1] = '\0';
      return true;
    };

    /**
     * @breif Cut N's from the end of read
     * @return True is read after cutting N is longer than or equal to min_length
     */
    bool removeNs(int min_length)
    {
      int l;
      for (l = strlen(m_strSeq) - 1; l >= min_length; --l) {
        if (m_strSeq[l] != 'N') {
	  break;
	}
      }

      //Ignore this if len < cutoff
      if (l < min_length) {
        return false;
      }

      m_strQual[l] = '\0';
      m_strSeq[l] = '\0';
      return true;
    };

    /**
     * @brief Check if entry contains given sequence.
     * @param seq -> Sequence to search.
     * @param mismatch -> maximum mismatches allowed.
     * @return True is seq found, false otherwise.
     */
    bool contains(std::string seq, int mismatch = 0)
    {
      int pos = 0;
      bool res = pol_util::find_substring(m_strSeq,seq,pos,mismatch);
      return res;
    };

  private:
    /**
     * @brief Constructor
     */
    FastqEntry(char entry[][MAX_LINE]){
      m_strName = new char[strlen(entry[0])+1];
      m_strSeq = new char[strlen(entry[1])+1];
      m_strPlus = new char[strlen(entry[2])+1];
      m_strQual = new char[strlen(entry[3])+1];

      strncpy(m_strName,entry[0],strlen(entry[0])+1);
      strncpy(m_strSeq,entry[1],strlen(entry[1])+1);
      strncpy(m_strPlus,entry[2],strlen(entry[2])+1);
      strncpy(m_strQual,entry[3],strlen(entry[3])+1);
    };

    ///> Sequence name
    char* m_strName;
    ///> Actual sequence of read
    char* m_strSeq;
    ///> Extra information line
    char* m_strPlus;
    ///> Quality string
    char* m_strQual;
  };

}

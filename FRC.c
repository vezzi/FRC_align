/*  
    qaTools - Just more qa tools.
    Copyright (C) 2011  P. Costea(paul.igor.costea@scilifelab.se)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <time.h>
#include <string>
#include "radix.h"
#include "sam.h"

typedef struct
{
  int doMedian,maxCoverage,minQual,maxInsert;
  bool spanCov,silent;
  FILE* detailed;
}Options;

#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)

#define EXIT_IF_NULL(P) \
  if (P == NULL) \
    return 1;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core)
{

  if (core->flag&BAM_FUNMAP) {
    return false;
  }

  return true;
}

/**
 * Print usage instructions
 */
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Version: 1.5\n");
  fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@scilifelab.se>\n\n");
  fprintf(stderr, "Usage:   qaCompute [options] <in.bam/sam> <output.out>\n");
  fprintf(stderr, "Options: \n");
  fprintf(stderr, "         -m            Also compute median coverage\n");
  fprintf(stderr, "         -q            Quality threshold. (min quality to consider) [1].\n");
  fprintf(stderr, "         -d            Print per-chromosome histogram [<output.out>.detail]\n");
  fprintf(stderr, "         -i            SIlent.Don't print too much stuff!\n");
  fprintf(stderr, "         -s [INT]      Compute 'span coverage' rather than base coverage, limiting insert size to INT. -1 -> consider all!\n");
  fprintf(stderr, "         -c [INT]      Maximum coverage to consider in histogram [30]\n");
  fprintf(stderr, "         -h [STR]      Use header from specified file. .sam OR .bam (Header must match the info in your input file. Otherwise, output is meaningless!)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: Input file should be sorted\n\n");
  return 1;
}

static void compute_print_cov(FILE* outputFile, Options userOpt, int* data, char* name,const uint32_t chrSize, uint64_t* coverageHist,const int currentTid)
{
  //clock_t start = clock();
  int32_t covVal = 0;
  uint64_t covSum = 0;
  uint32_t i;
  //Histogram vector.
  uint64_t* localCoverageHist = NULL;
  if (userOpt.detailed) {
    //Allocate
    localCoverageHist = new uint64_t[userOpt.maxCoverage+1];
    //Clear
    memset( localCoverageHist, 0, (userOpt.maxCoverage+1)*sizeof(uint64_t));
  }

  //Go through chromosome and count avarage covarage.
  for (i=0; i<chrSize; ++i){
    covVal += data[i];
    if (covVal < 0) {//int overrun!?
      fprintf(stderr,"Probably really good coverage, since variables overrun!\n");
    }
    //This will be sorted later.
    //If -m was not defined, this is useless, but cheaper than an 'if'
    data[i] = covVal;
    uint64_t prev = covSum;
    covSum += covVal;
    if (prev > covSum) {
      fprintf(stderr,"That's a big sum: %lu...\n",prev);
    }
    //Add value to histogram
    if (covVal > userOpt.maxCoverage) {
      ++coverageHist[userOpt.maxCoverage];
      if (localCoverageHist)
	++localCoverageHist[userOpt.maxCoverage];
    } else {
      ++coverageHist[covVal];
      if (localCoverageHist)
	++localCoverageHist[covVal];
    }

  }
  if (userOpt.doMedian)
    //Sort entireChr
    radix_sort(data, chrSize);

  if (localCoverageHist) {//Print details!
    fprintf(userOpt.detailed,"%s\t",name);
    for (int i=0; i<=userOpt.maxCoverage; ++i) {
        uint64_t coverage = 0;
        //All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition                                              
	for (int x = i; x<=userOpt.maxCoverage; ++x) coverage += localCoverageHist[x];
        fprintf(userOpt.detailed,"%3.5f\t",(double)(coverage)/chrSize*100);
    }
    fprintf(userOpt.detailed,"\n");
    fflush(userOpt.detailed);
    //Clean histogram!
    delete[] localCoverageHist;
    localCoverageHist = NULL;
  }

  //Printout avarage coverage over this chrom
  if (!userOpt.silent) {
    fprintf(stdout,"Coverage sum %lu ! \n", covSum);
    fprintf(stdout,"Average coverage over %s : %3.2f\n", name, (double)covSum / chrSize);
    if (userOpt.doMedian)
      fprintf(stdout,"Median coverage over %s : %d\n", name, data[chrSize/2]);
  }
  if (userOpt.doMedian)
    fprintf(outputFile, "%s\t%d\t%3.5f\t%d\n", name, chrSize, (double)covSum / chrSize, data[chrSize/2]);
  else
    fprintf(outputFile, "%s\t%d\t%3.5f\n", name, chrSize, (double)covSum / chrSize);

  //clock_t end = clock();
  //printf("time to compute this: %3.2f\n", (end-start)/CLOCKS_PER_SEC);
}

/**
 * Open a .sam/.bam file. 
 * @returns NULL is open failed.
 */
samfile_t * open_alignment_file(std::string path)
{
  samfile_t * fp = NULL;
  std::string flag = "r";
  if (path.substr(path.size()-3).compare("bam") == 0) {                                                                                                                                               
    //BAM file!                                                                                                                                    
    flag += "b";                                                                                                                                                                                                             
  }
  if ((fp = samopen(path.c_str(), flag.c_str() , 0)) == 0) {
    fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
  } 
  return fp;
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{
  samfile_t *fp;
  FILE *outputFile;
  Options userOpt;
  std::string headerFile = "";
  userOpt.doMedian = 0;
  userOpt.maxCoverage = 30;
  userOpt.spanCov = false;
  userOpt.silent = false;
  userOpt.detailed = NULL;
  userOpt.minQual = 1;
  userOpt.maxInsert = -1;
  bool doDetail = false;
  int arg;
  //Get args                                                                                                                                               
  while ((arg = getopt(argc, argv, "mdis:q:c:h:")) >= 0) {
    switch (arg) {
    case 'm': userOpt.doMedian = 1; break;
    case 'd': doDetail = true; break;
    case 'i': userOpt.silent = true; break;
    case 'q': userOpt.minQual = atoi(optarg); break;
    case 'c': userOpt.maxCoverage = atoi(optarg); break;
    case 'h': headerFile = optarg; break;
    case 's': userOpt.spanCov = true;
      userOpt.maxInsert = atoi(optarg);
      fprintf(stdout,"Max insert size %d\n",userOpt.maxInsert);
      break;
    }
  }

  if (argc-optind != 2) {
    print_usage();
    return 1;
  }

  bool outsideHeader = false;
  samfile_t * headerF = NULL;
  if (! headerFile.empty()) {
    outsideHeader = true;
    headerF = open_alignment_file(headerFile);
    EXIT_IF_NULL(headerF)
  }

  std::string alignFile(argv[optind]);
  fp = open_alignment_file(alignFile);
  EXIT_IF_NULL(fp);

  if (outsideHeader) {
    //Trick header!
    fp->header = headerF->header;
  }

  if ((outputFile = fopen(argv[optind+1], "wt")) == 0) {
    fprintf(stderr, "qaCompute: Filed to create output file %s\n", argv[optind+1]);
    return 1;
  }
  if (doDetail) {//Create detailed output file.
    std::string fName = argv[optind+1];
    fName += ".detail";
    userOpt.detailed = fopen(fName.c_str(),"wt");
    if (userOpt.detailed == NULL) {
      fprintf(stderr,"qaCompute: Unable to create detailed output file %s. No details will be printed!\n",fName.c_str());
    }
    fprintf(stdout,"Printing details in %s!\n",fName.c_str());
  }

    //Initialize bam entity
    bam1_t *b = bam_init1();

    //All var declarations
    uint64_t totalGenomeLength = 0;
    uint32_t unmappedReads = 0;
    uint32_t zeroQualityReads = 0;
    uint32_t totalNumberOfReads = 0;
    uint32_t totalProperPaires = 0;
    uint32_t chrSize = 0;
    uint32_t interChr = 0;
    uint32_t duplicates = 0;
    uint32_t usedReads = 0;
 
    int *entireChr = NULL;
    //Keep header for further reference
    bam_header_t* head = fp->header;
    
    int32_t currentTid = -1;

    //Create "map" vector for histogram
    uint64_t* coverageHist= (uint64_t*)malloc((userOpt.maxCoverage+1)*sizeof(uint64_t)); 
    memset( coverageHist, 0, (userOpt.maxCoverage+1)*sizeof(uint64_t));

    //Write file table header
    if (userOpt.doMedian == 1)
      fprintf(outputFile, "Chromosome\tSeq_len\tAvg_Cov\tMedian_Cov\n");
    else
      fprintf(outputFile, "Chromosome\tSeq_lem\tAvg_Cov\n");

    while (samread(fp, b) >= 0) {
      
      //uint32_t* cigar = bam1_cigar(b);
      //Get bam core.
      const bam1_core_t *core = &b->core;

      if (core == NULL) {
	//There is something wrong with the read/file
	printf("Input file is corrupt!");
	//Leak everything and exit!
	return -1;
      }

      //BAM block has been read
      if (!is_mapped(core))
	++unmappedReads;
      else {

	if (core->tid != currentTid) {
	  
	  //Count coverage!
	  if (currentTid != -1) {
	    if (!userOpt.silent)
	      fprintf(stdout,"Basing coverage on %u reads\n",usedReads);
	    usedReads = 0;
	    compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);
	  }

	  //Get length of next section                                                                                       
          chrSize = head->target_len[core->tid];
	  if (chrSize < 1) {//We can't have such sizes! this can't be right
	    fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid],chrSize);
	  }
          totalGenomeLength += chrSize;
	  if (!userOpt.silent)
	    fprintf(stdout,"Computing %s of size %u... \n",head->target_name[core->tid],chrSize);

	  //Done with current section.
	  //Allocate memory
	  entireChr = (int*)realloc(entireChr, (chrSize+1)*sizeof(int));
	  
	  if (entireChr == NULL) {
	    fprintf(stderr,"Allocation failed! \n");
	    return -1;
	  }
	  memset(entireChr, 0, (chrSize+1)*sizeof(int));
	  
	  currentTid = core->tid;
	
	}
	
	//If read has quality == 0, we won't count it as mapped
	if (core->qual >= userOpt.minQual) {
	 if (core->flag&BAM_FPROPER_PAIR) {
	    //Is part of a proper pair
	    ++totalProperPaires;
	  }

	 if (core->flag&BAM_FDUP) {
	   //This is a duplicate. Don't count it!.
	   ++duplicates;
	 } else {
	   if (!userOpt.spanCov) {
	     //All entries in SAM file are represented on the forward strand! (See specs of SAM format for details)    
	     ++entireChr[core->pos];
	     ++usedReads;
	     if ((uint32_t)(core->pos+core->l_qseq) >= chrSize)
	       --entireChr[chrSize-1];
	     else
	       --entireChr[core->pos+core->l_qseq];
	   } else {
	     //Computing span coverage. 
	     //Only consider first read in pair! and extend a bit to the end of the insert
	     if ((core->flag&BAM_FREAD1) //First in pair
		 && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
		 && (core->tid == core->mtid) /*Mate on the same chromosome*/
		 ) {
	       int32_t start = MIN(core->pos,core->mpos);
	       int32_t end = start+abs(core->isize);
	       int32_t iSize = end-start;
	       if ((userOpt.maxInsert == -1) || (iSize <= userOpt.maxInsert)) {
		 ++entireChr[start];
		 if ((uint32_t)end >= chrSize)
		   --entireChr[chrSize-1];
		 else
		   --entireChr[end];
		 ++usedReads;
	       }
	     } else if (core->tid != core->mtid) {
	       //Count inter-chrom mates
	       ++interChr;
	     }
	   }
	 }

	} else {
	  //Count is as unmapped?
	  ++zeroQualityReads;
	}
      }

      ++totalNumberOfReads;
      
    }

    //Compute coverage for the last "chromosome"
    compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);

    bam_destroy1(b);
    free(entireChr);

    fprintf(stdout,"\nDuplicates:%u \n", duplicates);

    //Print header for next table in output file
    fprintf(outputFile,"\nCov*X\tPercentage\tNr. of bases\n");

    fprintf(stdout,"Total genome lenght %lu \n", totalGenomeLength);
    //Compute procentages of genome cover!.
    int i;
    for (i=0; i<=userOpt.maxCoverage; ++i) {
      if (i == 0) {
	//Non-covered!
	fprintf(stdout,"%3.2f of genome has not been covered\n", (double)(coverageHist[i])/totalGenomeLength*100);
      } else {
	uint64_t coverage = 0;
	//All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition
	for (int x = i; x<=userOpt.maxCoverage; ++x) coverage += coverageHist[x];
	fprintf(stdout,"%3.2f of genome has been covered at least %dX \n", (double)(coverage)/totalGenomeLength*100, i);
	fprintf(outputFile,"%d\t%3.5f\t%lu\n",i, (double)(coverage)/totalGenomeLength*100, coverage);
      }

    }

    fprintf(outputFile,"\nOther\n");

    //Printout procentage of mapped/unmapped reads                                                                                                     
    double procentageOfUnmapped = 100*((double)unmappedReads/totalNumberOfReads);
    double procentageOfZeroQuality = 100*((double)zeroQualityReads/totalNumberOfReads);
    fprintf(outputFile,"Total number of reads: %u\n", totalNumberOfReads);
    fprintf(outputFile,"Total number of duplicates found and ignored: %u\n", duplicates);
    fprintf(outputFile,"Percentage of unmapped reads: %3.5f\n", procentageOfUnmapped);
    fprintf(outputFile,"Percentage of sub-par quality mappings: %3.5f\n", procentageOfZeroQuality);
    int32_t nrOfPaires = totalNumberOfReads/2;
    double procOfProperPaires = (double)(100*(double)totalProperPaires/2)/nrOfPaires;
    fprintf(outputFile,"Number of proper paired reads: %u\n", totalProperPaires);
    fprintf(outputFile,"Percentage of proper pairs: %3.5f\n", procOfProperPaires);
    if (userOpt.spanCov) {
      fprintf(outputFile, "Number of interchromosomal pairs: %u\n",interChr);
    }

    printf("Out of %u reads, you have %3.5f unmapped reads\n and %3.5f sub-par quality mappings\n", totalNumberOfReads ,procentageOfUnmapped, procentageOfZeroQuality);
    

    free(coverageHist);

  
    fclose(outputFile);
    if (outsideHeader) {
      //Must force this back to NULL, otherwise we'll delete it twice.
      fp->header = NULL;
      samclose(headerF);
    }
    samclose(fp);
    
    if (userOpt.detailed)
      fclose(userOpt.detailed);
  
  return 0;
}

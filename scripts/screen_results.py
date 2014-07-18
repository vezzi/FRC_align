import sys, os, glob
import argparse


def main(args):
    
    toBeMaskedElements = {};
    with open(args.bed_file) as fin:
        rows = ( line.split('\t') for line in fin )
        for row in rows:
            if not row[0].startswith("#"):
                try:
                    toBeMaskedElements[row[1]].append([row[2], row[3], row[7]])
                except :
                    toBeMaskedElements[row[1]] = [row[2], row[3], row[7]]
    #print toBeMaskedElements
    
    with open(args.variations) as fin:
        rows = ( line.rstrip().split('\t') for line in fin )
        for row in rows:
            chr_1       = row[0]
            chr_1_start = row[1]
            chr_1_end   = row[2]
            chr_1_supp  = row[3]
            chr_2       = row[5] ##TOBE CHANCGED
            chr_2_start = row[6] ##TOBE CHANCGED
            chr_2_end   = row[7] ##TOBE CHANCGED
            chr_2_supp  = row[8] ##TOBE CHANCGED

            expected_Links  = row[10] ##TOBE CHANCGED
            coverage        = row[11] ##TOBE CHANCGED
            observed_links  = row[12] ##TOBE CHANCGED
    
            ## now chategorize it
            categorized = 0
            if chr_1 in toBeMaskedElements:
                for interval in toBeMaskedElements[chr_1]:
                    if 
            
            print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            chr_1,chr_1_start,chr_1_end,chr_1_supp,chr_2,chr_2_start,chr_2_end, chr_2_supp, expected_Links, coverage, observed_links )

    return 0






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed-file', help="bed files containing the intervals to use as markers", type=str)
    parser.add_argument('--variations', help="tab file containing variations", type=str)
    args = parser.parse_args()

    main(args)
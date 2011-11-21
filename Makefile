# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -Wall -g
# Samtools path
SAMTOOLS=/Users/vezzi/Documents/workspace/samtools
INCLUDE=include/

all: FRC

FRC: FRC.o
	$(CC) FRC.o -o FRC -L$(SAMTOOLS) -lbam -lz #-lefence

FRC.o: FRC.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) FRC.c


clean:
	rm -rf *o FRC

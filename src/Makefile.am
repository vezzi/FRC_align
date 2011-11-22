ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=FRC FRC_debug


FRC_SOURCES=FRC.cpp \
	data_structures/Mask.cpp \
	errors/Data_Exception.cpp \
	errors/File_Not_Found.cpp \
	errors/Incorrect_Format.cpp \
	io/Auto_Unzip.cpp \
	io/Fasta.cpp \
	io/read_functions.cpp \
	data_structures/Reads.cpp 


FRC_debug_SOURCES=$(FRC_SOURCES)

SAMTOOLS=samtools/
INCLUDES = -I$(SAMTOOLS) 

FRC_CXXFLAGS= -O3 $(BOOST_CPPFLAGS)
FRC_debug_CXXFLAGS= -O0 -g3 $(BOOST_CPPFLAGS)
FRC_CFLAGS= -O3 
FRC_debug_CFLAGS= -O0 -g3 


FRC_LDFLAGS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS)  $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
FRC_LIBS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
FRC_debug_LDFLAGS = $(FRC_LDFLAGS)
FRC_debug_LIBS = $(FRC_LIBS)


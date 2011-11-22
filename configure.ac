dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(FRC, 0.3)
 
AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

AM_PROG_CC_C_O
AC_PROG_CXX
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_HEADERS([pthread.h stdio.h math.h limits.h iostream fstream string vector exception sstream bzlib.h],
	[],
	[AC_MSG_ERROR(["Header not found"])])
AC_PROG_INSTALL
LT_INIT

#ACX_MPI()

CXXFLAGS=" -Wall -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64"
CFLAGS=" -Wall -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64"
if test mpic++ = "$MPICXX"; then
	CC="$MPICC" 
	CXX="$MPICXX"
	LIBS="$MPILIBS $LIBS"
	CXXFLAGS+=" -DHAVE_MPI"
	CFLAGS+=" -DHAVE_MPI"
fi

BOOST_STATIC()
BOOST_REQUIRE([1.40.0])
BOOST_IOSTREAMS
BOOST_PROGRAM_OPTIONS
BOOST_THREADS


echo $ECHO_N "checking for 64 bit architecture... "
case "$build_cpu" in
x86_64*) CXXFLAGS+=" -m64 -DCPU_64" 
	echo "yes"
	;;
*) echo "no" ;;
esac

echo "======================================"

case "$build_cpu" in
x86_64*) echo "Compiling with 64 bit" ;;
*) 	echo "Compiling with 32 bit" ;;
esac
echo "======================================"

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT


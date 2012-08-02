dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(FRC_align, 1.0)

AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])


AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()
LT_INIT

# standard system requirements
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_HEADERS([pthread.h stdio.h math.h limits.h bzlib.h zlib.h \
	iostream fstream string vector exception sstream],
	[], [AC_MSG_ERROR(["Header not found"])] )
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# common flags
CXXFLAGS+=" -Wall"
CFLAGS+=" -Wall -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64"

# checking for 64 bit system (optional)
echo $ECHO_N "checking for 64 bit architecture... "
case "$build_cpu" in
x86_64*)
	CXXFLAGS+=" -m64 -DCPU_64"
	CFLAGS+=" -m64 -DCPU_64"
	COMPILING_64="yes"
	echo "yes"
	;;
*) echo "no" ;;
esac

# checking for BOOST library (required)
echo "===== BOOST support (required) ====="
BOOST_STATIC()
BOOST_REQUIRE([1.40.0])
BOOST_IOSTREAMS
BOOST_PROGRAM_OPTIONS
BOOST_THREADS

# checking for BOOST library (required)
echo "===== checking for samtools (required) ====="
AC_MSG_CHECKING(for samtools library)
AC_ARG_WITH([samtools], AC_HELP_STRING([--with-samtools=DIR], [Location of samtools header and library files (both in DIR) (optional)]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" = x; then
  SAMTOOLS_CFLAGS=""
  SAMTOOLS_LDFLAGS=""
  SAMTOOLS_LIBS=""
  AC_MSG_RESULT(disabled)
  samtools=false
else
  if test x$have_zlib = xno; then
     SAMTOOLS_CFLAGS=""
     SAMTOOLS_LDFLAGS=""
     SAMTOOLS_LIBS=""
     AC_MSG_RESULT(enabled, but zlib not enabled)
     samtools=false
  else
     AC_DEFINE(HAVE_SAMTOOLS,1,[Define to 1 if you have a working Samtools library.])
     SAMTOOLS_CFLAGS="-I$withval"
     SAMTOOLS_LDFLAGS="-L$withval"
     SAMTOOLS_LIBS="-lbam"
     AC_MSG_RESULT(enabled)
     samtools=true
  fi
fi
AC_SUBST(SAMTOOLS_CFLAGS)
AC_SUBST(SAMTOOLS_LDFLAGS)
AC_SUBST(SAMTOOLS_LIBS)
AM_CONDITIONAL([SAMTOOLS], [test "x$samtools" = xtrue])



# enable debug programs
AC_ARG_ENABLE([maintainer-mode],
     [AS_HELP_STRING([--enable-maintainer-mode],[create debug programs [no]], [28])],
     [enable_debug=$enableval],
     [enable_debug=no])
AM_CONDITIONAL(DEBUG_ENABLED,test "$enable_debug" != "no")

# enable static programs
AC_ARG_ENABLE([static-FRC],
     [AS_HELP_STRING([--enable-static-FRC],[create static FRC version [no]], [28])],
     [enable_static_FRC=$enableval],
     [enable_static_FRC=no])
AM_CONDITIONAL(STATIC_ENABLED,test "$enable_static_FRC" != "no")


# print summary
echo         "================ Summary ================="
if test "$COMPILING_64" = "yes" ; then
	echo "= Compiling with 64 bit                  =";
else
	echo "= Compiling with 32 bit                  =";
fi
if test "$enable_debug" != "no" ; then
	echo "= Compiling with debug programs          =";
else
	echo "= Compiling without debug programs       =";
fi
if test "$enable_static_FRC" != "no" ; then
	echo "= Compiling static FRC version  =";
else
	echo "= Not compiling static FRC version       =";
fi
echo         "=========================================="


AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT


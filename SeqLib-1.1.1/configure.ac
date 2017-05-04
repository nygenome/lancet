# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(seqkit, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/BamReader.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([disable])
##m4_include([m4/m4_ax_openmp.m4])

# Checks for programs.
AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

# check for c++11
#AX_CXX_COMPILE_STDCXX_11(,[optional])

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g $fail_on_warning -Wno-unknown-pragmas")
AC_SUBST(CXXFLAGS, "$CXXFLAGS")
AC_SUBST(CFLAGS, "$CFLAGS")
AC_SUBST(LDFLAGS, "$LDFLAGS")

# We always need to specify to link in external libraries
AC_SUBST(LIBS, "$LIBS -pthread")

AC_CONFIG_FILES([Makefile
	         src/Makefile
	         src/seqtools/Makefile])	

AC_OUTPUT

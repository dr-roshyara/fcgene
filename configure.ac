#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([fcgene], [1.0.7],[roshyra@yahoo.com])
AC_CONFIG_SRCDIR([src/fcgene.cpp])
AC_CONFIG_HEADER([src/config.h]) 
AC_INIT_AUTOMAKE 
AC_PROG_INSTALL
##########################################################################
# use the C++ compiler for the following checks
AC_LANG([C++])
##########################################################################
 AM_INIT_AUTOMAKE
##########################################################################
# Checks for programs.
AC_SUBST([AC_LDFLAGS])
##########################################################################
# check for XXX headers
##########################################################################
# store current CFLAGS and merge with AC_CXXFLAGS for compilation check
AC_CHECK_HEADERS([stdlib.h sys/timeb.h zlib.h]) 
# check for XXX library headers   
AC_MSG_CHECKING([for the MATH library headers])
# try to compile a file that includes a header of the library XXX
AC_COMPILE_IFELSE([[#include <math.h>]],
    [AC_MSG_RESULT([MATH Library was found])
        HAVE_MATH_H=1; AC_LDFLAGS="-lm";],
    [AC_MSG_RESULT([MATH LIBRARY NOT FOUND not found])
        HAVE_MATH_H=0;])

###
AC_MSG_CHECKING([for the ZLIB library headers])
AC_COMPILE_IFELSE([[#include "zlib.h"]],
    [AC_MSG_RESULT([ZLIB Library was found ])
        HAVE_ZLIB_H=1; AC_LDFLAGS="$AC_LDFLAGS -lz";],
     [AC_MSG_RESULT([ZLLIB was not found])
        HAVE_ZLIB_H=0;])


# handle results
if test $HAVE_ZLIB_H != 1; then
   echo "------------------------------------------------------------------------"
    AC_MSG_NOTICE([The ZLIB  library was not found. 
	That means fcGENE won't be able to read gzipped files. 
	If you want to read gezipped files, please first install zlib library 
	and then start ./configure again.!])
	echo "------------------------------------------------------------------------"
else 
  echo "------------------------------------------------------------------------"
    AC_MSG_NOTICE([ fcGENE can read and write gezipped files ending with *.gz.])
	echo "------------------------------------------------------------------------"

fi

##########################################################################


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_HEADER_MAJOR
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow])
AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

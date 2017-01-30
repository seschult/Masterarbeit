dnl Process this file with autoconf to produce a configure script.
dnl
dnl This file is free software; as a special exception the author gives
dnl unlimited permission to copy and/or distribute it, with or without
dnl modifications, as long as this notice is preserved.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
dnl implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

define([SVN_VERSION], esyscmd([sh -c "svnversion  |  tr -d '\n'"]))dnl

AC_INIT([BestCode],
	[0.SVN_VERSION],
	[andrii.verbytskyi@desy.de], [],
	[http://mpp.mpg.de])
AC_CONFIG_MACRO_DIR([m4])
dnl Must come before AM_INIT_AUTOMAKE.

dnl Minimum Autoconf version required.
AC_PREREQ([2.62])

AM_INIT_AUTOMAKE([1.9.6 -Wall    no-define])

dnl Checks for programs.
dnl We need a C compiler.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_LIBTOOL


AC_CONFIG_HEADERS([config.h])
#AX_CHECK_COMPILE_FLAG([-fdiagnostics-color=never], [CFLAGS="$CFLAGS -fdiagnostics-color=never"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wextra], [CFLAGS="$CFLAGS -Wextra"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wno-long-long], [CFLAGS="$CFLAGS -Wno-long-long"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wall], [CFLAGS="$CFLAGS -Wall"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Werror], [CFLAGS="$CFLAGS -Werror"], [], []) 
#AX_CHECK_COMPILE_FLAG([-pedantic], [CFLAGS="$CFLAGS -pedantic"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wno-error=long-long], [CFLAGS="$CFLAGS -Wno-error=long-long"], [], []) 


#AX_CHECK_COMPILE_FLAG([-fdiagnostics-color=never], [CXXFLAGS="$CXXFLAGS -fdiagnostics-color=never"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wextra], [CFLAGS="$CXXFLAGS -Wextra"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wno-long-long], [CXXFLAGS="$CFLAGS -Wno-long-long"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wall], [CXXFLAGS="$CXXFLAGS -Wall"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Werror], [CXXFLAGS="$CXXFLAGS -Werror"], [], []) 
#AX_CHECK_COMPILE_FLAG([-pedantic], [CXXFLAGS="$CXXFLAGS -pedantic"], [], []) 
#AX_CHECK_COMPILE_FLAG([-Wno-error=long-long], [CXXFLAGS="$CXXFLAGS -Wno-error=long-long"], [], []) 




AC_ARG_WITH([root],
  [AS_HELP_STRING([--with-root],
  [use root for builds, defaults to yes])],
  [ZBESTLIB_U="$withval"],
  [ZBESTLIB_U='yes'])
AM_CONDITIONAL(ZBESTLIB_USE, test "$ZBESTLIB_U" = yes)

if test x$ZBESTLIB_U = "xyes"; then
AC_DEFINE(USEROOT,1,"Use ROOT matrices and classes")
fi


AC_ARG_WITH([root],
  [AS_HELP_STRING([--with-root],
  [use root for builds, defaults to yes])],
  [ROOT_U="$withval"],
  [ROOT_U='yes'])
AM_CONDITIONAL(ROOT_USE, test "$ROOT_U" = yes)

if test x$ROOT_U = "xyes"; then
AC_DEFINE(USEROOT,1,"Use ROOT matrices and classes")
fi

AC_ARG_WITH([root-cflags],
  [AS_HELP_STRING([--with-root-cflags],
    [location of the ROOT headers, defaults to root-config --cflags])],
  [ROOT_CFLAGS="$withval"],
  [ROOT_CFLAGS=`root-config --cflags`])
AC_SUBST([ROOT_CFLAGS])




AC_ARG_WITH([root-libs],
  [AS_HELP_STRING([--with-root-libs],
    [location of the ROOT headers, defaults to root-config --libs -lMinuit])],
  [ROOT_LIBS="$withval"],
  [ROOT_LIBS=`root-config --libs`" -lMinuit "])
AC_SUBST([ROOT_LIBS])






dnl Copyright will apply as long as these sources are in use, e.g., are
dnl being compiled, which is reasonable year to claim the copyright.
AC_DEFINE([COPYRIGHT_YEAR], [m4_esyscmd([date +%Y])], [year in copyright message])
# Define these substitions here to keep all version information in one place.
# For information on how to properly maintain the library version information,
# refer to the libtool manual, section "Updating library version information":
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
AC_SUBST([EXAMPLE_SO_VERSION], [0:0:0])
AC_SUBST([EXAMPLE_API_VERSION], [${PACKAGE_VERSION}])
AC_SUBST([PACKAGE_URL], [http://mpp.mpg.de])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

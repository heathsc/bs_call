AC_INIT([GEMTools],[1.0])

AC_CHECK_LIB([z], [gzopen], [HAVE_ZLIB=1])
AC_SUBST([HAVE_ZLIB])
AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], [HAVE_BZLIB=1])
AC_SUBST([HAVE_BZLIB])

AC_CHECK_HEADERS([omp.h],[HAVE_OPENMP=1])
AC_SUBST([HAVE_OPENMP])

AC_ARG_WITH(gemtools,
[AS_HELP_STRING([--with-gemtools=PATH],
		[specify prefix directory for installed gemtool package.])])
if test "x$with_gemtools" != "x"; then
      gemtoolslibs="-L$with_gemtools/lib -lgemtools"
		gemtoolsinc="-I$with_gemtools/include -I$with_gemtools/resources/include"
fi

AC_SUBST([gemtoolsinc])
AC_SUBST([gemtoolslibs])

AC_CONFIG_FILES([GEMTools/Makefile.mk])
AC_CONFIG_FILES([src/Makefile.mk])
AC_OUTPUT

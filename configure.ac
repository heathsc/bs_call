AC_INIT([GEMTools],[1.0])

AC_CHECK_LIB([z], [gzopen], [HAVE_ZLIB=1])
AC_SUBST([HAVE_ZLIB])
AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], [HAVE_BZLIB=1])
AC_SUBST([HAVE_BZLIB])

AC_ARG_WITH(gsl,
[AS_HELP_STRING([--with-gsl=PATH],
		[specify prefix directory for installed gsl library.])])
if test "x$with_gsl" != "x"; then
      GSLLIBS="-L$with_gsl/lib -lgsl -lgslcblas"
		GSLINC="-I$with_gsl/include"
else
      GSLLIBS="-lgsl -lgslcblas"
		GSLINC=""
fi

AC_SUBST([GSLLIBS])
AC_SUBST([GSLINC])

AC_ARG_WITH(htslib,
[AS_HELP_STRING([--with-htslib=PATH],
		[specify prefix directory for installed htslib library.])])
if test "x$with_htslib" != "x"; then
      HTSLIBS="-L$with_htslib/lib -lhts"
		HTSINC="-I$with_htslib/include"
else
      HTSLIBS="-lhts"
		HTSINC=""
fi

AC_SUBST([HTSLIBS])
AC_SUBST([HTSINC])

AC_CONFIG_FILES([gt/Makefile.mk])
AC_CONFIG_FILES([src/Makefile.mk])
AC_OUTPUT

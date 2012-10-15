dnl Process this file with autoconf to produce a configure script
AC_INIT(DESCRIPTION)

R_SCMD="${R_HOME}/bin/Rscript -e"

dnl Rcpp
RCPP_CXXFLAGS=`${R_SCMD} "Rcpp:::CxxFlags()"`
RCPP_LDFLAGS=`${R_SCMD} "Rcpp:::LdFlags()"`

dnl pbdSLAP
EXT_LDFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('EXT_LDFLAGS','"${R_ARCH}"','pbdSLAP')"`
R_SLAP=`${R_SCMD} "source('./R/get_lib.r');get.lib('R_SLAP','"${R_ARCH}"')"`
if test -f "$R_SLAP" ; then
  SLAP_LDFLAGS="$R_SLAP \$(LAPACK_LIBS) \$(BLAS_LIBS)"
else
  SLAP_LDFLAGS="$EXT_LDFLAGS"
fi

dnl pbdMPI
SPMD_CPPFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('PKG_CPPFLAGS','"${R_ARCH}"')"`
SPMD_LDFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('PKG_LIBS','"${R_ARCH}"')"`

dnl Start to substitute templates
AC_SUBST(RCPP_CXXFLAGS)
AC_SUBST(RCPP_LDFLAGS)
AC_SUBST(SLAP_LDFLAGS)
AC_SUBST(SPMD_CPPFLAGS)
AC_SUBST(SPMD_LDFLAGS)
AC_OUTPUT(src/Makevars)

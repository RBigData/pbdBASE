## autoconf configure.ac && autoreconf -if
AC_PREREQ([2.69])
AC_INIT([pbdBASE], [0.3-0])
AC_CONFIG_SRCDIR([src/base_scalapack.c])

R_SCMD="${R_HOME}/bin/Rscript -e"


# SEXPtools
SEXPTOOLS_CPPFLAGS=`${R_SCMD} "SEXPtools:::cppflags()"`
SEXPTOOLS_LDFLAGS=`${R_SCMD} "SEXPtools:::ldflags()"`


# pbdSLAP
EXT_LDFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('EXT_LDFLAGS','"${R_ARCH}"','pbdSLAP')"`
R_SLAP=`${R_SCMD} "source('./R/get_lib.r');get.lib('R_SLAP','"${R_ARCH}"')"`
if test -f "$R_SLAP" ; then
  SLAP_LDFLAGS="$R_SLAP \$(LAPACK_LIBS) \$(BLAS_LIBS)"
else
  SLAP_LDFLAGS="$EXT_LDFLAGS"
fi


# pbdMPI
SPMD_CPPFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('PKG_CPPFLAGS','"${R_ARCH}"')"`
SPMD_LDFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('PKG_LIBS','"${R_ARCH}"')"`


# Report
echo " "
echo "****************** Results of pbdBASE package configure *****************"
echo " "
echo ">> SEXPTOOLS_CPPFLAGS = ${SEXPTOOLS_CPPFLAGS}"
echo ">> SEXPTOOLS_LDFLAGS = ${SEXPTOOLS_LDFLAGS}"
echo ">> EXT_LDFLAGS = ${EXT_LDFLAGS}"
echo ">> R_SLAP = ${R_SLAP}"
echo ">> SLAP_LDFLAGS = ${SLAP_LDFLAGS}"
echo ">> SPMD_CPPFLAGS = ${SPMD_CPPFLAGS}"
echo ">> SPMD_LDFLAGS = ${SPMD_LDFLAGS}"
echo " "
echo "*************************************************************************"
echo " "


# Substitute templates
AC_SUBST(SEXPTOOLS_CPPFLAGS)
AC_SUBST(SEXPTOOLS_LDFLAGS)
AC_SUBST(SLAP_LDFLAGS)
AC_SUBST(SPMD_CPPFLAGS)
AC_SUBST(SPMD_LDFLAGS)

AC_OUTPUT(src/Makevars)


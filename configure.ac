## autoconf configure.ac && autoreconf -if
AC_PREREQ([2.69])
AC_INIT([pbdBASE], [0.3-0])
AC_CONFIG_SRCDIR([src/base_scalapack.c])

AC_PROG_CC_C99


# Determine R_HOME.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

R_SCMD="${R_HOME}/bin/Rscript -e"


### Ability to disable src/export_blacs/ --- needed for Intel mpi
AC_ARG_ENABLE(blacsexport,
  AC_HELP_STRING(
    [--enable-blacsexport],
    [Enable BLACS export, default yes. Disable for Intel MPI.]
  ),
  [ENABLE_BLACSEXPORT="no"],
  [ENABLE_BLACSEXPORT="yes"]
)


if test "X$ENABLE_BLACSEXPORT" = "Xyes"; then
  EXPORT_BLACS="export_blacs/pkg_tools.o"
else
  EXPORT_BLACS=""
fi


# pbdSLAP
EXT_LDFLAGS=`${R_SCMD} "source('./R/get_conf.r');get.conf('EXT_LDFLAGS','"${R_ARCH}"','pbdSLAP')"`
R_SLAP=`${R_SCMD} "source('./R/get_lib.r');get.lib('R_SLAP','"${R_ARCH}"')"`
if test -f "$R_SLAP"; then
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
echo ">> EXT_LDFLAGS = ${EXT_LDFLAGS}"
echo ">> EXPORT_BLACS = ${EXPORT_BLACS}"
echo ">> R_SLAP = ${R_SLAP}"
echo ">> SLAP_LDFLAGS = ${SLAP_LDFLAGS}"
echo ">> SPMD_CPPFLAGS = ${SPMD_CPPFLAGS}"
echo ">> SPMD_LDFLAGS = ${SPMD_LDFLAGS}"
echo " "
echo "*************************************************************************"
echo " "


# Substitute templates
AC_SUBST(SLAP_LDFLAGS)
AC_SUBST(SPMD_CPPFLAGS)
AC_SUBST(SPMD_LDFLAGS)
AC_SUBST(EXPORT_BLACS)

AC_OUTPUT(src/Makevars)

#' Functions to set and get BLACS_APTS
#'
#' To set and get BLACS array/object/whatever pointers needed in and from R.
#' Because other packages has it's own memory stack vision that may not be
#' visiable by this package or vice versa.
#' 
#' The `set.blacs.apts()` is for advanced users. This one is needed to be
#' called within R from `pbdBASE` package to set the pointers to the memory
#' where BLACS had initialized so that the pointers are set to the right
#' address of the memory stack.
#'
#' The `get.blacs.apts()` is for debugging only. The advanced user mainly calls
#' the C version `get_BLACS_APTS_from_R()` in `src/export_blacs/pkg_ools.c`.
#'
#' I am lazy to use .C(), but should not hurt performance here.
#' Eventually, .pbdBASEEnv should pass to .C() and set/get pointers from it
#' instead of .GlobalEnv.
#' 
#' @name blacs_apts
#' @rdname blacs_apts
#' @export
set.blacs.apts <- function(){
  .C("set_BLACS_APTS_in_R", PACKAGE="pbdBASE")
  invisible()
} # End of set.blacs.apts()

#' @name blacs_apts
#' @rdname blacs_apts
#' @export
get.blacs.apts <- function(){
  .C("get_BLACS_APTS_from_R", PACKAGE="pbdBASE")
  invisible()
} # End of get.blacs.apts()


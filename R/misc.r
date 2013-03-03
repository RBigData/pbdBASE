# 'next best divisor'
base.nbd <- function(n, d)
{
  .Call("R_nbd", as.integer(n), as.integer(d), PACKAGE="pbdBASE")
}


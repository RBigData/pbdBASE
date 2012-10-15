# ##################################################
# --------------------------------------------------
# Classes
# --------------------------------------------------
# ##################################################

# Distributed Dense Matrix
ddmatrix <- 
  setClass(
           "ddmatrix", 
            representation(
                           Data="matrix",
                           dim="numeric",
                           ldim="numeric",
                           bldim="numeric",
                           CTXT="numeric"
            ),
            prototype(
                      Data=matrix(0),
                      dim=c(1,1),
                      ldim=c(1,1),
                      bldim=c(1,1),
                      CTXT=0
            )
)



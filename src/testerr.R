cat(.libPaths(), "\n")
#cat("No error")

#stop ("Error!\n")
pck.rcq <- require("Rcurl")
if (!pck.rcq)
    stop ("Package missing\n")

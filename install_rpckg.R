args=(commandArgs(TRUE))
eval(parse(text=args[[1]]))
cat("RLIB=", RLIB, "\n")

if (!("caTools" %in% installed.packages(lib.loc=RLIB)))
    install.packages("caTools", repos="http://cran.us.r-project.org", lib=RLIB)

devtools::document()
devtools::build()
# devtools::check()
system("R CMD check --as-cran ../eventPred_0.2.2.tar.gz")


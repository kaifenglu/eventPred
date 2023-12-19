devtools::document()
devtools::build()
# devtools::check()
system("R CMD check --as-cran ../eventPred_0.2.3.tar.gz")


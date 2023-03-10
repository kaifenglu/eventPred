devtools::document()
devtools::build()

system("R CMD check --as-cran ../eventPred_0.0.2.tar.gz")







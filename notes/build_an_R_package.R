devtools::document()
devtools::build()

system("R CMD check --as-cran ../eventPred_0.2.1.tar.gz")


devtools::document()
devtools::build()

system("R CMD check --as-cran ../eventPred_0.1.3.tar.gz")


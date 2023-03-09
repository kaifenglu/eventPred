#' Observed enrollment and event data
#'
#' Simulated enrollment and event data used in the examples.
#'
#' @format ## `observedData`
#' A data frame with 300 rows and 5 columns:
#' \describe{
#'   \item{\code{randdt}}{Randomization date}
#'   \item{\code{cutoffdt}}{Cutoff date of the data set}
#'   \item{\code{time}}{Day of event or censoring since randomization}
#'   \item{\code{event}}{Event indicator: 1 for event, 0 for non-event}
#'   \item{\code{dropout}}{Dropout indicator: 1 for dropout, 0 for non-dropout}
#' }
#' For ongoing subjects, both \code{event} and \code{dropout} are equal to 0.
"observedData"

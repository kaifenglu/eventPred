#' Observed enrollment and event data
#'
#' Simulated enrollment and event data used in the examples.
#'
#' @format ## `observedData`
#' A data frame with 300 rows and 5 columns:
#' \describe{
#'   \item{\code{randdt}}{The randomization date}
#'   \item{\code{cutoffdt}}{The cutoff date}
#'   \item{\code{time}}{The day of event or censoring since randomization}
#'   \item{\code{event}}{The event indicator: 1 for event, 0 for non-event}
#'   \item{\code{dropout}}{The dropout indicator: 1 for dropout, 0 for non-dropout}
#' }
#' For ongoing subjects, both \code{event} and \code{dropout} are equal to 0.
"observedData"

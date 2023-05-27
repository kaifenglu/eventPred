#' Interim enrollment and event data before enrollment completion
#'
#' A data frame with 225 rows and 7 columns:
#' \describe{
#'   \item{\code{trialsdt}}{The trial start date}
#'   \item{\code{randdt}}{The randomization date}
#'   \item{\code{cutoffdt}}{The cutoff date}
#'   \item{\code{treatment}}{The treatment group}
#'   \item{\code{time}}{The day of event or censoring since randomization}
#'   \item{\code{event}}{The event indicator: 1 for event, 0 for non-event}
#'   \item{\code{dropout}}{The dropout indicator: 1 for dropout,
#'   0 for non-dropout}
#' }
#' For ongoing subjects, both \code{event} and \code{dropout} are equal to 0.
"interimData1"


#' Interim enrollment and event data after enrollment completion
#'
#' A data frame with 300 rows and 7 columns:
#' \describe{
#'   \item{\code{trialsdt}}{The trial start date}
#'   \item{\code{randdt}}{The randomization date}
#'   \item{\code{cutoffdt}}{The cutoff date}
#'   \item{\code{treatment}}{The treatment group}
#'   \item{\code{time}}{The day of event or censoring since randomization}
#'   \item{\code{event}}{The event indicator: 1 for event, 0 for non-event}
#'   \item{\code{dropout}}{The dropout indicator: 1 for dropout,
#'   0 for non-dropout}
#' }
#' For ongoing subjects, both \code{event} and \code{dropout} are equal to 0.
"interimData2"


#' Final enrollment and event data after achieving the target number of events
#'
#' A data frame with 300 rows and 7 columns:
#' \describe{
#'   \item{\code{trialsdt}}{The trial start date}
#'   \item{\code{randdt}}{The randomization date}
#'   \item{\code{cutoffdt}}{The cutoff date}
#'   \item{\code{treatment}}{The treatment group}
#'   \item{\code{time}}{The day of event or censoring since randomization}
#'   \item{\code{event}}{The event indicator: 1 for event, 0 for non-event}
#'   \item{\code{dropout}}{The dropout indicator: 1 for dropout,
#'   0 for non-dropout}
#' }
#' For ongoing subjects, both \code{event} and \code{dropout} are equal to 0.
"finalData"

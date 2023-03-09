#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction by utilizing
#'   observed data and specified enrollment and event models.
#'
#' @param df The observed subject-level enrollment and event data,
#'   including \code{randdt}, \code{cutoffdt}, \code{time}, \code{event},
#'   and \code{dropout}. If not provided, defaults to \code{NULL} for
#'   enrollment and event prediction at the design stage.
#' @param target_n The target number of subjects to enroll in the study.
#' @param target_d The target number of events to reach in the study.
#' @param to_predict Specified what to predict: enrollment only, event
#'   only, or enrollment and event. If not provided, defaults to
#'   enrollment and event.
#' @param enroll_model The available enrollment models are Poisson,
#'   time-decay, and B-spline. If not provided, defaults to B-spine,
#'   in which case the number of knots and day lags need to be specified
#'   to compute the average enrollment rate to carry forward.
#' @param nknots The number of inner knots for the B-spline enrollment
#'   model. If not provided, defaults 1.
#' @param lags The day lags to compute the average enrollment rate to
#'   carry forward for the B-spline enrollment model. If not provided,
#'   defaults to 30.
#' @param parameter_enroll_model The enrollment model parameters for
#'   design stage enrollment prediction.
#' @param event_model The available event models are exponential,
#'   Weibull, log-normal, piecewise exponential, and model averaging (of
#'   Weibull and log-normal). If the piecewise exponential model is
#'   chosen, the number of pieces needs to be specified. If the event
#'   model is not provided, it is set to model averaging.
#' @param npieces The number of pieces for the piecewise exponential
#'   event model. If not provided, defaults to 3.
#' @param parameter_event_model The event model parameters for design
#'   stage event prediction.
#' @param dropout_model The available dropout models are none,
#'   exponential, Weibull, and log-normal. If not provided, defaults
#'   to Weibull.
#' @param parameter_dropout_model The dropout model parameters for
#'   design stage event prediction.
#' @param fixedFollowup A Boolean variable indicating whether a fixed
#'   follow-up design is used. If not provided, defaults to \code{FALSE}
#'   for a variable follow-up design.
#' @param followupTime The follow-up time for a fixed follow-up design,
#'   in days. If not provided, defaults to 365.
#' @param pilevel The prediction interval level. If not provided,
#'   defaults to 0.90.
#' @param nreps The number of replications for simulation. If not
#'   provided, defaults to 500.
#'
#' @details
#' For the time-decay model, the mean function is
#' \code{mu(t) = (mu/delta) (t - (1/delta)(1 - exp(-delta*t)))}
#' and the rate function is
#' \code{lambda(t) = (mu/delta) (1 - exp(-delta*t))}
#' For the B-spline model, the daily enrollment rate is approximated using
#' the exponential of a B-spline function:
#' \code{lambda(t) = exp(B(t)*theta)},
#' where B(t) represents the B-spline basis functions.
#'
#' For the \code{parameter_enroll_model}, it can be used to specify
#' a piecewise Poisson model that is parameterized through
#' \code{accrualTime} and \code{accrualIntensity} (which are treated
#' as fixed). For the homogeneous Poisson and time-decay models,
#' \code{parameter_enroll_model} is used to specify the prior
#' distribution of model parameters, with a very small variance
#' being used to fix the parameter values. It should be noted
#' that the B-spline model is not appropriate for use during
#' the design stage.
#'
#' For the \code{parameter_event_model}, it should be a list that
#' includes the \code{model} used to specify the event process,
#' \code{ngroups} to indicate the number of treatment groups,
#' \code{prob} to indicate the randomization probabilities
#' for each group, \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix, both of which
#' have \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the \code{j}-th
#' treatment group. For the piecewise exponential event model,
#' this should also include \code{knots} to indicate the location
#' of inner knots. It is important to note that the
#' \code{model averaging} event model cannot be used at the
#' design stage.
#'
#' For the \code{parameter_dropout_model}, it should be a list that
#' includes the \code{model} to specify the dropout process,
#' \code{ngroups} to indicate the number of treatment groups,
#' \code{prob} to indicate the randomization probabilities
#' for each group, \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix, both of which
#' have \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the \code{j}-th
#' treatment group.
#'
#' For analysis stage enrollment and event prediction, the
#' \code{parameter_enroll_model}, \code{parameter_event_model}, and
#' \code{parameter_dropout_model} are either set to \code{NULL} to
#' use the observed data only, or specify the prior distribution
#' of model parameters to be combined with observed data likelihood
#' for enhanced modeling flexibility.
#'
#' @return A list that includes the fits of observed data models,
#' as well as simulated enrollment data for new subjects and
#' simulated event data for ongoing and new subjects.
#'
#' @examples
#'
#' # Example 1: Enrollment prediction at analysis stage
#'
#' pred <- getPrediction(
#'   df = observedData, target_n = 480,
#'   to_predict = "enrollment only",
#'   enroll_model = "b-spline", nknots = 1, lags = 30,
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' # Example 2: Enrollment prediction at design stage
#'
#' parameter_enroll_model <- list(
#'   model = "piecewise poisson",
#'   accrualTime = seq(0, 8)*30.4375,
#'   accrualIntensity = 26/9*seq(1, 9)/30.4375)
#'
#' pred <- getPrediction(
#'   target_n = 480,
#'   to_predict = "enrollment only",
#'   parameter_enroll_model = parameter_enroll_model,
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' # Example 3: Event prediction at analysis stage after enrollment ends
#'
#' pred <- getPrediction(
#'   df = observedData, target_d = 200,
#'   to_predict = "event only",
#'   event_model = "piecewise exponential", npieces = 3,
#'   dropout_model = "exponential",
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' # Example 4: Event prediction at analysis stage before enrollment ends
#'
#' pred <- getPrediction(
#'   df = observedData, target_n = 480, target_d = 200,
#'   to_predict = "enrollment and event",
#'   enroll_model = "b-spline", nknots = 1, lags = 30,
#'   event_model = "piecewise exponential", npieces = 3,
#'   dropout_model = "exponential",
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' # Example 5: Event prediction at design stage
#'
#' parameter_enroll_model <- list(
#'   model = "piecewise poisson",
#'   accrualTime = seq(0, 8)*30.4375,
#'   accrualIntensity = 26/9*seq(1, 9)/30.4375)
#'
#' parameter_event_model <- list(
#'   model = "piecewise exponential",
#'   ngroups = 2,
#'   prob = c(0.5, 0.5),
#'   theta = log(c(0.0533, 0.0309, 0.0533, 0.0533)/30.4375),
#'   vtheta = diag(4)*1e-8,
#'   knots = 6*30.4375)
#'
#' parameter_dropout_model <- list(
#'   model = "exponential",
#'   ngroups = 2,
#'   prob = c(0.5, 0.5),
#'   theta = log(rep(-log(1-0.05)/12, 2)/30.4375),
#'   vtheta = diag(2)*1e-8)
#'
#' pred <- getPrediction(
#'   target_n = 480, target_d = 200,
#'   to_predict = "enrollment and event",
#'   parameter_enroll_model = parameter_enroll_model,
#'   parameter_event_model = parameter_event_model,
#'   parameter_dropout_model = parameter_dropout_model,
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' @export
#'
getPrediction <- function(
    df = NULL, target_n = NA, target_d = NA,
    to_predict = "enrollment and event",
    enroll_model = "B-spline", nknots = 1, lags = 30,
    parameter_enroll_model = NULL,
    event_model = "model averaging", npieces = 3,
    parameter_event_model = NULL,
    dropout_model = "weibull",
    parameter_dropout_model = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nreps = 500) {

  if (!is.null(df)) erify::check_class(df, "data.frame")
  if (!is.na(target_n)) erify::check_n(target_n)
  if (!is.na(target_d)) erify::check_n(target_d)
  if (is.na(target_n) & is.na(target_d))
    stop("At least one of target_n and target_d must be specified.")
  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))
  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline"))
  erify::check_n(nknots)
  erify::check_n(lags, zero=TRUE)
  if (!is.null(parameter_enroll_model))
    erify::check_class(parameter_enroll_model, "list")
  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))
  erify::check_n(npieces)
  if (!is.null(parameter_event_model))
    erify::check_class(parameter_event_model, "list")
  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-normal"))
  if (!is.null(parameter_dropout_model))
    erify::check_class(parameter_dropout_model, "list")
  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)

  if (!is.null(df)) {
    names(df) <- tolower(names(df))
    trialsdt = min(df$randdt)
    cutoffdt = df$cutoffdt[1]

    # summarize observed data
    observed <- summarizeObserved(df, to_predict, dropout_model)
  }

  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(df = observed$adsl, enroll_model, nknots)

      # combine prior and likelihood to yield posterior
      if (!is.null(parameter_enroll_model)) {
        if (tolower(enroll_model) == "poisson") {
          enroll_fit$theta <-
            1/(1/enroll_fit$vtheta +
                 1/parameter_enroll_model$vtheta)*
            (1/enroll_fit$vtheta*enroll_fit$theta +
               1/parameter_enroll_model$vtheta*parameter_enroll_model$theta)
          enroll_fit$vtheta <-
            1/(1/enroll_fit$vtheta +
                 1/parameter_enroll_model$vtheta)
        } else {
          enroll_fit$theta <-
            solve(solve(enroll_fit$vtheta) +
                    solve(parameter_enroll_model$vtheta),
                  solve(enroll_fit$vtheta, enroll_fit$theta) +
                    solve(parameter_enroll_model$vtheta,
                          parameter_enroll_model$theta))
          enroll_fit$vtheta <- solve(solve(enroll_fit$vtheta) +
                                       solve(parameter_enroll_model$vtheta))
        }
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        target_n, df = observed$adsl, enroll_fit = enroll_fit,
        lags, pilevel, nreps, showplot = FALSE)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        target_n, df = NULL, enroll_fit = parameter_enroll_model,
        lags, pilevel, nreps, showplot = FALSE)
    }
  }

  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      event_fit <- fitEvent(df = observed$adtte, event_model, npieces)

      # combine prior and likelihood to yield posterior
      if (!is.null(parameter_event_model)) {
        if (tolower(event_model) == "exponential") {
          event_fit$theta <- 1/(1/event_fit$vtheta +
                                  1/parameter_event_model$vtheta)*
            (1/event_fit$vtheta*event_fit$theta +
               1/parameter_event_model$vtheta*parameter_event_model$theta)
          event_fit$vtheta <- 1/(1/event_fit$vtheta +
                                   1/parameter_event_model$vtheta)
        } else {
          event_fit$theta <-
            solve(solve(event_fit$vtheta) +
                    solve(parameter_event_model$vtheta),
                  solve(event_fit$vtheta, event_fit$theta) +
                    solve(parameter_event_model$vtheta,
                          parameter_event_model$theta))
          event_fit$vtheta <- solve(solve(event_fit$vtheta) +
                                      solve(parameter_event_model$vtheta))
        }
      }


      # whether to include dropout model
      if (tolower(dropout_model) != "none") {
        dropout_fit <- fitDropout(df = observed$adtte, dropout_model)

        # combine prior and likelihood to yield posterior
        if (!is.null(parameter_dropout_model)) {
          if (tolower(dropout_model) == "exponential") {
            dropout_fit$theta <-
              1/(1/dropout_fit$vtheta +
                   1/parameter_dropout_model$vtheta)*
              (1/dropout_fit$vtheta*dropout_fit$theta +
                 1/parameter_dropout_model$vtheta*
                 parameter_dropout_model$theta)
            dropout_fit$vtheta <-
              1/(1/dropout_fit$vtheta +
                   1/parameter_dropout_model$vtheta)
          } else {
            dropout_fit$theta <-
              solve(solve(dropout_fit$vtheta) +
                      solve(parameter_dropout_model$vtheta),
                    solve(dropout_fit$vtheta, dropout_fit$theta) +
                      solve(parameter_dropout_model$vtheta,
                            parameter_dropout_model$theta))
            dropout_fit$vtheta <-
              solve(solve(dropout_fit$vtheta) +
                      solve(parameter_dropout_model$vtheta))
          }
        }


        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            target_d, df = observed$adtte,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit, dropout_fit = dropout_fit,
            fixedFollowup, followupTime, pilevel, nreps,
            showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            target_d, df = observed$adtte, newSubjects = NULL,
            event_fit = event_fit, dropout_fit = dropout_fit,
            fixedFollowup, followupTime, pilevel, nreps,
            showplot = FALSE)
        }
      } else {
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            target_d, df = observed$adtte,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit, dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nreps,
            showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            target_d, df = observed$adtte, newSubjects = NULL,
            event_fit = event_fit, dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nreps,
            showplot = FALSE)
        }
      }
    } else { # event prediction at design stage
      if (!is.null(parameter_dropout_model)) {
        event_pred <- predictEvent(
          target_d, df = NULL, newSubjects = enroll_pred$newSubjects,
          event_fit = parameter_event_model,
          dropout_fit = parameter_dropout_model,
          fixedFollowup, followupTime, pilevel, nreps,
          showplot = FALSE)
      } else {
        event_pred <- predictEvent(
          target_d, df = NULL, newSubjects = enroll_pred$newSubjects,
          event_fit = parameter_event_model, dropout_fit = NULL,
          fixedFollowup, followupTime, pilevel, nreps,
          showplot = FALSE)
      }
    }
  }


  # output results
  if (!is.null(df)) { # analysis stage prediction
    if (tolower(to_predict) == "enrollment only") {
      dfs <- enroll_pred$plotdata

      # separate data into observed and predicted
      dfa <- dfs %>% filter(is.na(.data$lower))
      dfb <- dfs %>% filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g2 <- flabel(dfs, trialsdt)

      # plot the enrollment data with month as x-axis label
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_step(data=dfa, aes(x=.data$date, y=.data$n), color="black") +
        geom_line(data=dfb, aes(x=.data$date, y=.data$n), color="blue") +
        geom_vline(xintercept = cutoffdt, linetype = 2) +
        scale_x_date(name = NULL,
                     labels = scales::date_format("%b"),
                     breaks = scales::breaks_width(bw),
                     minor_breaks = NULL,
                     expand = c(0.01, 0.01)) +
        labs(y = "Subjects", title = "Predicted subject enrollment") +
        theme_bw()

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)

      list(observed = observed, enroll_fit = enroll_fit,
           enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      dfs <- event_pred$plotdata

      # separate data into observed and predicted
      dfa <- dfs %>% filter(is.na(.data$lower))
      dfb <- dfs %>% filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g2 <- flabel(dfs, trialsdt)

      # plot the enrollment and time to event data with month as x-axis label
      # generate the plot
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_step(data=dfa, aes(x=.data$date, y=.data$n), color="black") +
        geom_line(data=dfb, aes(x=.data$date, y=.data$n), color="blue") +
        geom_vline(xintercept = cutoffdt, linetype = 2) +
        geom_hline(yintercept = target_d, linetype = 2) +
        scale_x_date(name = NULL,
                     labels = scales::date_format("%b"),
                     breaks = scales::breaks_width(bw),
                     minor_breaks = NULL,
                     expand = c(0.01, 0.01)) +
        labs(y = "Events", title = "Predicted events") +
        theme_bw()

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)

      if (tolower(dropout_model) != "none") {
        list(observed = observed, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(observed = observed, event_fit = event_fit,
             event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      df1 <- enroll_pred$plotdata %>%
        mutate(parameter = "subjects")
      df1last <- last(df1)

      df2 <- event_pred$plotdata %>%
        mutate(parameter = "events")
      df2last <- last(df2)

      # extend enrollment prediction time to event prediction time
      df3 <- df2last %>%
        mutate(n = df1last$n,
               lower = df1last$lower,
               upper = df1last$upper,
               parameter = df1last$parameter)

      dfs <- df1 %>%
        bind_rows(df3) %>%
        bind_rows(df2)

      # separate data into observed and predicted
      dfa <- dfs %>% filter(is.na(.data$lower))
      dfb <- dfs %>% filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g2 <- flabel(dfs, trialsdt)

      # plot the enrollment and time to event data with month as x-axis label
      # generate the plot
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$date, ymin=.data$lower, ymax=.data$upper,
                        group=.data$parameter),
                    alpha=0.5, fill="lightblue") +
        geom_step(data=dfa, aes(x=.data$date, y=.data$n,
                                group=.data$parameter), color="black") +
        geom_line(data=dfb, aes(x=.data$date, y=.data$n,
                                group=.data$parameter), color="blue") +
        geom_vline(xintercept = cutoffdt, linetype = 2) +
        geom_hline(yintercept = target_d, linetype = 2) +
        scale_x_date(name = NULL,
                     labels = scales::date_format("%b"),
                     breaks = scales::breaks_width(bw),
                     minor_breaks = NULL,
                     expand = c(0.01, 0.01)) +
        labs(y = "Subjects / Events",
             title = "Predicted cumulative subjects and events over time") +
        theme_bw()

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)

      if (tolower(dropout_model) != "none") {
        list(observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             event_pred = event_pred)
      }
    }
  } else { # design stage prediction
    if (tolower(to_predict) == "enrollment only") {
      dfb <- enroll_pred$plotdata

      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$time, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_line(data=dfb, aes(x=.data$time, y=.data$n), color="blue") +
        scale_x_continuous(name = "Days since randomization",
                           expand = c(0.01, 0.01)) +
        labs(y = "Subjects", title = "Predicted subject enrollment") +
        theme_bw()

      print(g1)

      list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      dfb <- event_pred$plotdata

      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$time, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_line(data=dfb, aes(x=.data$time, y=.data$n), color="blue") +
        geom_hline(yintercept = target_d, linetype = 2) +
        scale_x_continuous(name = "Days since randomization",
                           expand = c(0.01, 0.01)) +
        labs(y = "Events", title = "Predicted events") +
        theme_bw()

      print(g1)

      if (!is.null(parameter_dropout_model)) {
        list(event_fit = parameter_event_model,
             dropout_fit = parameter_dropout_model, event_pred = event_pred)
      } else {
        list(event_fit = parameter_event_model, event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      df1 <- enroll_pred$plotdata %>%
        mutate(parameter = "subjects")
      df1last <- last(df1)

      df2 <- event_pred$plotdata %>%
        mutate(parameter = "events")
      df2last <- last(df2)

      # extend enrollment prediction time to event prediction time
      df3 <- df2last %>%
        mutate(n = df1last$n,
               lower = df1last$lower,
               upper = df1last$upper,
               parameter = df1last$parameter)

      dfb <- df1 %>%
        bind_rows(df3) %>%
        bind_rows(df2)

      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$time, ymin=.data$lower, ymax=.data$upper,
                        group=.data$parameter),
                    alpha=0.5, fill="lightblue") +
        geom_line(data=dfb, aes(x=.data$time, y=.data$n,
                                group=.data$parameter), color="blue") +
        geom_hline(yintercept = target_d, linetype = 2) +
        scale_x_continuous(name = "Days since randomization",
                           expand = c(0.01, 0.01)) +
        labs(y = "Subjects / Events",
             title = "Predicted cumulative subject and events over time") +
        theme_bw()

      print(g1)


      if (!is.null(parameter_dropout_model)) {
        list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred,
             event_fit = parameter_event_model,
             dropout_fit = parameter_dropout_model, event_pred = event_pred)
      } else {
        list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred,
             event_fit = parameter_event_model, event_pred = event_pred)
      }
    }
  }

}



library(shiny)
library(shinyMatrix)
library(shinyFeedback)
library(shinyjs)
library(shinybusy)
library(readxl)
library(openxlsx)
library(dplyr)
library(prompter)
library(plotly)
library(lubridate)
library(eventPred)

# determine the placement of major ticks of x-axis
fbw <- function(n_months) {
  if (n_months <= 15) {
    bw = 'M1'
  } else if (n_months <= 30) {
    bw = 'M2'
  } else if (n_months <= 45) {
    bw = 'M3'
  } else if (n_months <= 90) {
    bw = 'M6'
  } else {
    bw = 'M12'
  }

  bw
}


observedPanel <- tabPanel(
  title = "Observed Data",
  value = "observed_data_panel",


  htmlOutput("statistics"),

  plotlyOutput("cum_accrual_plot"),

  conditionalPanel(
    condition = "input.stage != 'Real-time after enrollment completion'",
    plotlyOutput("daily_accrual_plot")
  ),

  conditionalPanel(
    condition = "input.to_predict == 'Enrollment and event' ||
    input.stage == 'Real-time after enrollment completion'",

    plotlyOutput("event_km_plot"),

    plotlyOutput("dropout_km_plot")
  ),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",
    dataTableOutput("inputdf"))
)


enrollmentPanel <- tabPanel(
  title = "Enrollment Model",
  value = "enroll_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "enroll_model_parameter",
        "Which enrollment model to use?",
        choices = c("Poisson",
                    "Time-decay",
                    "Piecewise Poisson"),
        selected = "Piecewise Poisson",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.enroll_model_parameter == 'Poisson'",

               numericInput(
                 "poisson_rate",
                 "Daily enrollment rate",
                 value = 1,
                 min = 0, max = 100, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_model_parameter == 'Time-decay'",

               fluidRow(
                 column(6,
                        numericInput(
                          "mu",
                          "Base rate, mu",
                          value = 1,
                          min = 0, max = 100, step = 1)
                 ),

                 column(6,
                        numericInput(
                          "delta",
                          "Decay rate, delta",
                          value = 1,
                          min = 0, max = 100, step = 1)
                 )
               )
             ),

             conditionalPanel(
               condition = "input.enroll_model_parameter ==
               'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "piecewise_poisson_rate",
                 "Daily enrollment rate by time interval",
                 value = matrix(c(0,1), ncol = 2,
                                dimnames = list("Interval 1",
                                                c("Starting time",
                                                  "Enrollment rate"))),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_piecewise_poisson_rate",
                            label=NULL, icon=icon("plus")),
               actionButton("del_piecewise_poisson_rate",
                            label=NULL, icon=icon("minus"))
             )
      )
    )
  ),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "enroll_model",
        "Which enrollment model to use?",
        choices = c("Poisson",
                    "Time-decay",
                    "B-spline",
                    "Piecewise Poisson"),
        selected = "B-spline",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.enroll_model == 'B-spline'",

               numericInput(
                 "nknots",
                 "How many inner knots to use?",
                 value = 0,
                 min = 0, max = 10, step = 1),

               numericInput(
                 "lags",
                 paste("How many days before cutoff to",
                       "average the enrollment rate over for prediction?"),
                 value = 30,
                 min = 0, max = 365, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_model == 'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "accrualTime",
                 "What is the starting time of each time interval?",
                 value = matrix(0, ncol = 1,
                                dimnames = list("Interval 1",
                                                "Starting time")),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_accrualTime",
                            label=NULL, icon=icon("plus")),
               actionButton("del_accrualTime",
                            label=NULL, icon=icon("minus"))
             )
      )

    ),

    plotlyOutput("enroll_fit")
  )
)


eventPanel <- tabPanel(
  title = "Event Model",
  value = "event_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "event_model_parameter",
        "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Piecewise exponential",
        inline = FALSE)
      ),


      column(6,
             conditionalPanel(
               condition = "input.event_model_parameter == 'Exponential'",

               shinyMatrix::matrixInput(
                 "exponential_survival",
                 "Hazard rate for each group",
                 value = matrix(c(0.0030, 0.0044),
                                nrow = 1,
                                dimnames = list(
                                  NULL, c("Treatment hazard rate",
                                          "Control hazard rate"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=FALSE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),

             conditionalPanel(
               condition = "input.event_model_parameter == 'Weibull'",

               shinyMatrix::matrixInput(
                 "weibull_survival",
                 "Weibull parameters",
                 value = matrix(c(1.42, 1.42, 392, 254),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list(
                                  c("Shape", "Scale"),
                                  c("Treatment group", "Control group"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),


             conditionalPanel(
               condition = "input.event_model_parameter == 'Log-normal'",

               shinyMatrix::matrixInput(
                 "lnorm_survival",
                 "Log-normal parameters",
                 value = matrix(c(5.4, 5, 1, 1),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list(
                                  c("Mean on log scale", "SD on log scale"),
                                  c("Treatment group", "Control group"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),


             conditionalPanel(
               condition = "input.event_model_parameter ==
               'Piecewise exponential'",

               shinyMatrix::matrixInput(
                 "piecewise_exponential_survival",
                 "Hazard rate by time interval for each group",
                 value = matrix(c(0, 0.0030, 0.0044),
                                nrow = 1,
                                dimnames = list(
                                  "Interval 1",
                                  c("Starting time",
                                    "Treatment hazard rate",
                                    "Control hazard rate"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               ),

               actionButton("add_piecewise_exponential_survival",
                            label=NULL, icon=icon("plus") ),
               actionButton("del_piecewise_exponential_survival",
                            label=NULL, icon=icon("minus"))
             )
      )
    )
  ),


  conditionalPanel(
    condition = "input.stage != 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "event_model",
        "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-normal",
                    "Piecewise exponential",
                    "Model averaging"),
        selected = "Model averaging",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.event_model == 'Piecewise exponential'",

               shinyMatrix::matrixInput(
                 "piecewiseSurvivalTime",
                 "What is the starting time of each time interval?",
                 value = matrix(0, ncol = 1,
                                dimnames = list("Interval 1",
                                                "Starting time")),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_piecewiseSurvivalTime",
                            label=NULL, icon=icon("plus")),
               actionButton("del_piecewiseSurvivalTime",
                            label=NULL, icon=icon("minus"))
             )
      )
    ),

    plotlyOutput("event_fit")
  )
)



dropoutPanel <- tabPanel(
  title = "Dropout Model",
  value = "dropout_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "dropout_model_parameter",
        "Which time-to-dropout model to use?",
        choices = c("None",
                    "Exponential",
                    "Weibull",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Exponential",
        inline = FALSE)
      ),

      column(6,

             conditionalPanel(
               condition = "input.dropout_model_parameter == 'Exponential'",

               shinyMatrix::matrixInput(
                 "exponential_dropout",
                 "Hazard rate for each group",
                 value = matrix(c(0.0003, 0.0003),
                                nrow = 1,
                                dimnames = list(
                                  NULL, c("Treatment hazard rate",
                                          "Control hazard rate"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=FALSE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),

             conditionalPanel(
               condition = "input.dropout_model_parameter == 'Weibull'",

               shinyMatrix::matrixInput(
                 "weibull_dropout",
                 "Weibull parameters",
                 value = matrix(c(1.25, 1.25, 2000, 1000),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list(
                                  c("Shape", "Scale"),
                                  c("Treatment group", "Control group"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),


             conditionalPanel(
               condition = "input.dropout_model_parameter == 'Log-normal'",

               shinyMatrix::matrixInput(
                 "lnorm_dropout",
                 "Log-normal parameters",
                 value = matrix(c(10, 8, 2.64, 2.64),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list(
                                  c("Mean on log scale", "SD on log scale"),
                                  c("Treatment group", "Control group"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               )
             ),


             conditionalPanel(
               condition = "input.dropout_model_parameter ==
               'Piecewise exponential'",

               shinyMatrix::matrixInput(
                 "piecewise_exponential_dropout",
                 "Hazard rate by time interval for each group",
                 value = matrix(c(0, 0.0003, 0.0003),
                                nrow = 1,
                                dimnames = list(
                                  "Interval 1",
                                  c("Starting time",
                                    "Treatment hazard rate",
                                    "Control hazard rate"))
                 ),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)
               ),

               actionButton("add_piecewise_exponential_dropout",
                            label=NULL, icon=icon("plus") ),
               actionButton("del_piecewise_exponential_dropout",
                            label=NULL, icon=icon("minus"))
             )
      )
    )
  ),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "dropout_model",
        "Which time-to-dropout model to use?",
        choices = c("None",
                    "Exponential",
                    "Weibull",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Weibull",
        inline = FALSE)
      ),


      column(6,
             conditionalPanel(
               condition = "input.dropout_model == 'Piecewise exponential'",

               shinyMatrix::matrixInput(
                 "piecewiseDropoutTime",
                 "What is the starting time of each time interval?",
                 value = matrix(0, ncol = 1,
                                dimnames = list("Interval 1",
                                                "Starting time")),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_piecewiseDropoutTime",
                            label=NULL, icon=icon("plus")),
               actionButton("del_piecewiseDropoutTime",
                            label=NULL, icon=icon("minus"))
             )
      )
    ),

    plotlyOutput("dropout_fit")
  )
)




predictPanel <- tabPanel(
  title = "Prediction Results",
  value = "prediction_results_panel",




  htmlOutput("pred_date"),
  plotlyOutput("pred_plot"),

  downloadButton("downloadSumdata", "Download summary data")
)



# user interface ----------------
ui <- fluidPage(

  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  prompter::use_prompt(),

  add_busy_spinner(),

  titlePanel("Enrollment and Event Prediction"),


  sidebarLayout(
    sidebarPanel(

      radioButtons(
        "stage",
        "Stage of the study",
        choices = c("Design stage",
                    "Real-time before enrollment completion",
                    "Real-time after enrollment completion"),
        selected = "Real-time after enrollment completion",
        inline = FALSE),

      conditionalPanel(
        condition = "input.stage == 'Design stage' ||
                 input.stage == 'Real-time before enrollment completion'",

        radioButtons(
          "to_predict",
          "What to predict?",
          choices = c("Enrollment only",
                      "Enrollment and event"),
          selected = "Enrollment only",
          inline = FALSE)
      ),



      conditionalPanel(
        condition =
          "input.stage == 'Real-time after enrollment completion'",

        radioButtons(
          "to_predict2",
          "What to predict?",
          choices = c("Event only"),
          selected = "Event only",
          inline = FALSE)
      ),




      fluidRow(
        column(6,
               conditionalPanel(
                 condition = "input.stage == 'Design stage' ||
                 input.stage == 'Real-time before enrollment completion'",

                 numericInput(
                   "target_n",
                   "Target enrollment #",
                   value = 500,
                   min = 1, max = 20000, step = 1)
               )
        ),

        column(6,
               conditionalPanel(
                 condition = "input.to_predict == 'Enrollment and event' ||
               input.stage == 'Real-time after enrollment completion'",

                 numericInput(
                   "target_d",
                   "Target event #",
                   value = 300,
                   min = 1, max = 10000, step = 1)
               )
        )
      ),


      conditionalPanel(
        condition = "input.stage != 'Design stage'",

        fileInput(
          "file1",
          "Upload subject level data",
          accept = ".xlsx"
        )
      ),


      fluidRow(
        column(6, radioButtons(
          "pilevel",
          "Prediction interval",

          choices = c("95%" = "0.95", "90%" = "0.90", "80%" = "0.80"),
          selected = "0.90",
          inline = TRUE)
        ),

        column(6, numericInput(
          "nyears",
          label = tags$span(
            "Years after cutoff",
            tags$span(icon(name = "question-circle")) %>%
              add_prompt(message = "Upper limit of x-axis on prediction plot",
                         position = "right")),
          value = 4,
          min = 1, max = 10, step = 1)
        )
      ),


      conditionalPanel(
        condition = "input.to_predict == 'Enrollment and event' ||
        input.stage == 'Real-time after enrollment completion'",

        checkboxGroupInput(
          "to_show",
          "What to show on prediction plot?",
          choices = c("Enrollment", "Event", "Dropout", "Ongoing"),
          selected = c("Enrollment", "Event"),
          inline = TRUE
        )
      ),


      fluidRow(
        column(6, numericInput(
          "nreps",
          label = "Simulation runs",
          value = 100,
          min = 100, max = 2000, step = 1)
        ),

        column(6, numericInput(
          "seed",
          label = "Seed",
          value = 2023,
          min = 0, max = 100000, step = 1
        ))
      ),


      fluidRow(
        column(6,
               conditionalPanel(
                 condition = "input.stage == 'Design stage'",

                 numericInput(
                   "allocation_ratio",
                   label = tags$span(
                     "Allocation ratio",
                     tags$span(icon(name = "question-circle")) %>%
                       add_prompt(message = "Treatment to control",
                                  position = "right")),
                   value=1, min=0.1, max=10, step=0.01)
               )
        ),

        column(6, actionButton(
          "predict", "Predict",
          style="color: #fff; background-color: #337ab7;
          border-color: #2e6da4")
        )
      )
    ),


    mainPanel(

      tabsetPanel(
        id = "results",
        observedPanel,
        enrollmentPanel,
        eventPanel,
        dropoutPanel,
        predictPanel
      )
    )
  )
)




# server function -------------
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })

  # whether to show or hide the observed data panel
  observeEvent(input$stage, {
    if (input$stage != 'Design stage') {
      showTab(inputId = "results", target = "observed_data_panel")
    } else {
      hideTab(inputId = "results", target = "observed_data_panel")
    }
  })


  # what to predict at different stages
  to_predict <- reactive({
    if (input$stage != 'Real-time after enrollment completion') {
      input$to_predict
    } else {
      input$to_predict2
    }
  })


  # whether to show or hide enrollment, event, and dropout panels
  observeEvent(to_predict(), {
    if (to_predict() == 'Enrollment only') {
      showTab(inputId = "results", target = "enroll_model_panel")
      hideTab(inputId = "results", target = "event_model_panel")
      hideTab(inputId = "results", target = "dropout_model_panel")
    } else if (to_predict() == 'Enrollment and event') {
      showTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
      showTab(inputId = "results", target = "dropout_model_panel")
    } else if (to_predict() == 'Event only') {
      hideTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
      showTab(inputId = "results", target = "dropout_model_panel")
    }
  })


  target_n <- reactive({
    req(input$target_n)
    valid = (input$target_n > 0 && input$target_n == round(input$target_n))
    shinyFeedback::feedbackWarning(
      "target_n", !valid,
      "Target enrollment must be a positive integer")
    req(valid)
    as.numeric(input$target_n)
  })


  target_d <- reactive({
    req(input$target_d)
    valid1 = (input$target_d > 0 && input$target_d == round(input$target_d))
    shinyFeedback::feedbackWarning(
      "target_d", !valid1,
      "Target events must be a positive integer")

    if (input$stage == "Enrollment and event") {
      valid2 = (input$target_d <= input$target_n)
      shinyFeedback::feedbackWarning(
        "target_d", !valid2,
        "Target events must be less than or equal to target enrollment")
    } else {
      valid2 = (input$target_d <= observed()$n0)
      shinyFeedback::feedbackWarning(
        "target_d", !valid2,
        "Target events must be less than or equal to sample size")
    }

    req(valid1 && valid2)
    as.numeric(input$target_d)
  })


  nyears <- reactive({
    req(input$nyears)
    valid = (input$nyears > 0)
    shinyFeedback::feedbackWarning(
      "nyears", !valid,
      "Years after cutoff must be a positive number")
    req(valid)
    as.numeric(input$nyears)
  })


  showEnrollment <- reactive({
    "Enrollment" %in% input$to_show
  })


  showEvent <- reactive({
    "Event" %in% input$to_show
  })


  showDropout <- reactive({
    "Dropout" %in% input$to_show
  })


  showOngoing <- reactive({
    "Ongoing" %in% input$to_show
  })


  nreps <- reactive({
    req(input$nreps)
    valid = (input$nreps > 0 && input$nreps == round(input$nreps))
    shinyFeedback::feedbackWarning(
      "nreps", !valid,
      "Number of simulations must be a positive integer")
    req(valid)
    as.numeric(input$nreps)
  })


  allocation_ratio <- reactive({
    req(input$allocation_ratio)
    valid = (input$allocation_ratio > 0)
    shinyFeedback::feedbackWarning(
      "allocation_ratio", !valid,
      "Allocation ratio must be a positive number")
    req(valid)
    as.numeric(input$allocation_ratio)
  })


  poisson_rate <- reactive({
    req(input$poisson_rate)
    valid = (input$poisson_rate > 0)
    shinyFeedback::feedbackWarning(
      "poisson_rate", !valid,
      "Daily enrollment rate must be a positive number")
    req(valid)
    as.numeric(input$poisson_rate)
  })


  mu <- reactive({
    req(input$mu)
    valid = (input$mu > 0)
    shinyFeedback::feedbackWarning(
      "mu", !valid,
      "Base rate must be a positive number")
    req(valid)
    as.numeric(input$mu)
  })


  delta <- reactive({
    req(input$delta)
    valid = (input$delta > 0)
    shinyFeedback::feedbackWarning(
      "delta", !valid,
      "Decay rate must be a positive number")
    req(valid)
    as.numeric(input$delta)
  })


  piecewise_poisson_rate <- reactive({
    req(input$piecewise_poisson_rate)
    t = as.numeric(input$piecewise_poisson_rate[,1])
    lambda = as.numeric(input$piecewise_poisson_rate[,2])

    valid1 = all(diff(t) > 0) && (t[1] == 0)
    if (!valid1) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }

    valid2 = all(lambda >= 0)
    if (!valid2) {
      showNotification(
        "Enrollment rate must be nonnegative"
      )
    }

    valid3 = any(lambda > 0)
    if (!valid3) {
      showNotification(
        "At least one enrollment rate must be positive"
      )
    }

    req(valid1 && valid2 && valid3)

    matrix(c(t, lambda), ncol = 2)
  })


  nknots <- reactive({
    req(input$nknots)
    valid = (input$nknots >= 0 && input$nknots == round(input$nknots))
    shinyFeedback::feedbackWarning(
      "nknots", !valid,
      "Number of inner knots must be a nonnegative integer")
    req(valid)
    as.numeric(input$nknots)
  })


  lags <- reactive({
    req(input$lags)
    valid = (input$lags >= 0 && input$lags == round(input$lags))
    shinyFeedback::feedbackWarning(
      "lags", !valid,
      "Number of day lags must be a nonnegative integer")
    req(valid)
    as.numeric(input$lags)
  })


  accrualTime <- reactive({
    t = as.numeric(input$accrualTime)
    valid = all(diff(t) > 0) && (t[1] == 0)
    if (!valid) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  exponential_survival <- reactive({
    lambda = as.numeric(input$exponential_survival)
    valid = all(lambda > 0)
    if (!valid) {
      showNotification(
        "Hazard rate must be positive"
      )
    }
    req(valid)
    lambda
  })


  weibull_survival <- reactive({
    shape = as.numeric(input$weibull_survival[1,])
    scale = as.numeric(input$weibull_survival[2,])

    valid1 = all(shape > 0)
    if (!valid1) {
      showNotification(
        "Weibull shape parameter must be positive"
      )
    }

    valid2 = all(scale > 0)
    if (!valid2) {
      showNotification(
        "Weibull scale parameter must be positive"
      )
    }

    req(valid1 && valid2)

    matrix(c(shape, scale), nrow = 2, byrow = TRUE)
  })


  lnorm_survival <- reactive({
    meanlog = as.numeric(input$lnorm_survival[1,])
    sdlog = as.numeric(input$lnorm_survival[2,])

    valid = all(sdlog > 0)
    if (!valid) {
      showNotification(
        "SD on the log scale must be positive"
      )
    }

    req(valid)

    matrix(c(meanlog, sdlog), nrow = 2, byrow = TRUE)
  })


  piecewise_exponential_survival <- reactive({
    t = as.numeric(input$piecewise_exponential_survival[,1])
    lambda = as.numeric(input$piecewise_exponential_survival[,-1])

    valid = all(lambda > 0)
    if (!valid) {
      showNotification(
        "Hazard rate must be positive"
      )
    }

    req(valid)

    matrix(c(t, lambda), nrow = length(t))
  })


  piecewiseSurvivalTime <- reactive({
    t = as.numeric(input$piecewiseSurvivalTime)
    valid = all(diff(t) > 0) && (t[1] == 0)
    if (!valid) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  exponential_dropout <- reactive({
    lambda = as.numeric(input$exponential_dropout)
    valid = all(lambda > 0)
    if (!valid) {
      showNotification(
        "Hazard rate must be positive"
      )
    }
    req(valid)
    lambda
  })


  weibull_dropout <- reactive({
    shape = as.numeric(input$weibull_dropout[1,])
    scale = as.numeric(input$weibull_dropout[2,])

    valid1 = all(shape > 0)
    if (!valid1) {
      showNotification(
        "Weibull shape parameter must be positive"
      )
    }

    valid2 = all(scale > 0)
    if (!valid2) {
      showNotification(
        "Weibull scale parameter must be positive"
      )
    }

    req(valid1 && valid2)

    matrix(c(shape, scale), nrow = 2, byrow = TRUE)
  })


  lnorm_dropout <- reactive({
    meanlog = as.numeric(input$lnorm_dropout[1,])
    sdlog = as.numeric(input$lnorm_dropout[2,])

    valid = all(sdlog > 0)
    if (!valid) {
      showNotification(
        "SD on the log scale must be positive"
      )
    }

    req(valid)

    matrix(c(meanlog, sdlog), nrow = 2, byrow = TRUE)
  })


  piecewise_exponential_dropout <- reactive({
    t = as.numeric(input$piecewise_exponential_dropout[,1])
    lambda = as.numeric(input$piecewise_exponential_dropout[,-1])

    valid = all(lambda > 0)
    if (!valid) {
      showNotification(
        "Hazard rate must be positive"
      )
    }

    req(valid)

    matrix(c(t, lambda), nrow = length(t))
  })


  piecewiseDropoutTime <- reactive({
    t = as.numeric(input$piecewiseDropoutTime)
    valid = all(diff(t) > 0) && (t[1] == 0)
    if (!valid) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  # input data set
  df <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    df <- readxl::read_excel(inFile$datapath)

    if (to_predict() == "Enrollment only") {
      required_columns <- c('randdt', 'cutoffdt')
    } else {
      required_columns <- c('randdt', 'cutoffdt', 'time', 'event', 'dropout')
    }

    column_names <- colnames(df)

    shiny::validate(
      need(all(required_columns %in% column_names),
           "You don't have the right data")
    )

    df
  })


  # summarize observed data
  observed <- reactive({
    if (!is.null(df()))
      summarizeObserved(df(), to_predict(), showplot = FALSE)
  })

  # enrollment fit
  enroll_fit <- reactive({
    if (!is.null(df()))
      fitEnrollment(df(), input$enroll_model, nknots(),
                    accrualTime(), showplot = FALSE)
  })

  # event fit
  event_fit <- reactive({
    if (!is.null(df()))
      fitEvent(df(), input$event_model,
               piecewiseSurvivalTime(), showplot = FALSE)
  })

  # dropout fit
  dropout_fit <- reactive({
    if (!is.null(df()) && input$dropout_model != "None") {
      shiny::validate(need(observed()$c0 > 0, paste(
        "The number of dropouts must be positive to fit a dropout model.")))
      fitDropout(df(), input$dropout_model,
                 piecewiseDropoutTime(), showplot = FALSE)
    }
  })


  # enrollment and event prediction
  pred <- eventReactive(input$predict, {
    set.seed(as.numeric(input$seed))

    if (to_predict() != "Enrollment only") {
      shiny::validate(need(
        showEnrollment() || showEvent() || showDropout() || showOngoing(),
        "Need at least one parameter to show on prediction plot"))
    }

    if (input$stage == "Design stage") {

      ngroups = 2
      allocation_ratio = allocation_ratio()
      prob = c(allocation_ratio, 1)/(allocation_ratio + 1)

      # enroll model specifications
      if (input$enroll_model_parameter == "Poisson") {
        theta = log(poisson_rate())
        vtheta = 1e-8
      } else if (input$enroll_model_parameter == "Time-decay") {
        theta = c(log(mu()), log(delta()))
        vtheta = diag(2)*1e-8
      } else if (input$enroll_model_parameter == "Piecewise Poisson") {
        theta = log(piecewise_poisson_rate()[,2])
        vtheta = diag(length(theta))*1e-8
        accrualTime = piecewise_poisson_rate()[,1]
      }

      enroll_model_parameter <- list(
        model = input$enroll_model_parameter,
        theta = theta, vtheta = vtheta)

      if (input$enroll_model_parameter == "Piecewise Poisson") {
        enroll_model_parameter$accrualTime = accrualTime
      }


      # event model specifications
      if (to_predict() == "Enrollment and event" ||
          to_predict() == "Event only") {
        if (input$event_model_parameter == "Exponential") {
          theta = log(exponential_survival())
        } else if (input$event_model_parameter == "Weibull") {
          theta = c(log(weibull_survival()[,1]),
                    log(weibull_survival()[,2]))
        } else if (input$event_model_parameter == "Log-normal") {
          theta = c(lnorm_survival()[1,1],
                    log(lnorm_survival()[2,1]),
                    lnorm_survival()[1,2],
                    log(lnorm_survival()[2,2]))
        } else if (input$event_model_parameter ==
                   "Piecewise exponential") {
          theta = c(
            log(piecewise_exponential_survival()[,2]),
            log(piecewise_exponential_survival()[,3]))
          piecewiseSurvivalTime =
            piecewise_exponential_survival()[,1]
        }

        event_model_parameter = list(
          ngroups = ngroups,
          prob = prob,
          model = input$event_model_parameter,
          theta = theta,
          vtheta = diag(length(theta))*1e-8)

        if (input$event_model_parameter == "Piecewise exponential") {
          event_model_parameter$piecewiseSurvivalTime = piecewiseSurvivalTime
        }


        # dropout model specifications
        if (input$dropout_model_parameter != "None") {
          if (input$dropout_model_parameter == "Exponential") {
            theta = log(exponential_dropout())
          } else if (input$dropout_model_parameter == "Weibull") {
            theta = c(log(weibull_dropout()[,1]),
                      log(weibull_dropout()[,2]))
          } else if (input$dropout_model_parameter == "Log-normal") {
            theta = c(lnorm_dropout()[1,1],
                      log(lnorm_dropout()[2,1]),
                      lnorm_dropout()[1,2],
                      log(lnorm_dropout()[2,2]))
          } else if (input$dropout_model_parameter ==
                     "Piecewise exponential") {
            theta = c(
              log(piecewise_exponential_dropout()[,2]),
              log(piecewise_exponential_dropout()[,3]))
            piecewiseDropoutTime =
              piecewise_exponential_dropout()[,1]
          }


          dropout_model_parameter = list(
            ngroups = ngroups,
            prob = prob,
            model = input$dropout_model_parameter,
            theta = theta,
            vtheta = diag(length(theta))*1e-8
          )

          if (input$dropout_model_parameter == "Piecewise exponential") {
            dropout_model_parameter$piecewiseDropoutTime =
              piecewiseDropoutTime
          }
        } else {
          dropout_model_parameter = NULL
        }
      }




      # get prediction results based on what to predict
      if (to_predict() == "Enrollment only") {
        getPrediction(
          to_predict = to_predict(),
          target_n = target_n(),
          enroll_model_parameter = enroll_model_parameter,
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE)
      } else if (to_predict() == "Enrollment and event") {
        getPrediction(
          to_predict = to_predict(),
          target_n = target_n(),
          target_d = target_d(),
          enroll_model_parameter = enroll_model_parameter,
          event_model_parameter = event_model_parameter,
          dropout_model_parameter = dropout_model_parameter,
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showDropout = showDropout(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE)
      } else if (to_predict() == "Event only") {
        getPrediction(
          to_predict = to_predict(),
          target_d = target_d(),
          event_model_parameter = event_model_parameter,
          dropout_model_parameter = dropout_model_parameter,
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showDropout = showDropout(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE)
      }
    } else { # real-time prediction
      shiny::validate(
        need(!is.null(df()),
             "Please upload data for real-time prediction."))

      if (to_predict() == "Enrollment only") {
        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_n = target_n(),
          enroll_model = input$enroll_model,
          nknots = nknots(),
          lags = lags(),
          accrualTime = accrualTime(),
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE)
      } else if (to_predict() == "Enrollment and event") {
        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        if (input$dropout_model != "None") {
          shiny::validate(need(observed()$c0 > 0, paste(
            "The number of dropouts must be positive",
            "to fit a dropout model.")))
        }

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_n = target_n(),
          target_d = target_d(),
          enroll_model = input$enroll_model,
          nknots = nknots(),
          accrualTime = accrualTime(),
          lags = lags(),
          event_model = input$event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          dropout_model = input$dropout_model,
          piecewiseDropoutTime = piecewiseDropoutTime(),
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showDropout = showDropout(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE)
      } else if (to_predict() == "Event only") {
        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        if (input$dropout_model != "None") {
          shiny::validate(need(observed()$c0 > 0, paste(
            "The number of dropouts must be positive",
            "to fit a dropout model.")))
        }

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_d = target_d(),
          event_model = input$event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          dropout_model = input$dropout_model,
          piecewiseDropoutTime = piecewiseDropoutTime(),
          pilevel = as.numeric(input$pilevel),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showDropout = showDropout(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE)
      }
    }
  })


  output$statistics <- renderText({
    if (!is.null(observed())) {
      str1 <- paste("Trial start date:", observed()$trialsdt)
      str2 <- paste("Data cutoff date:", observed()$cutoffdt)
      str3 <- paste("Days since trial start:", observed()$t0)
      str4 <- paste("Current number of subjects:", observed()$n0)

      if (to_predict() == "Enrollment and event" ||
          to_predict() == "Event only") {
        str5 <- paste("Current number of events:", observed()$d0)
        str6 <- paste("Current number of dropouts:", observed()$c0)
        str7 <- paste("Number of ongoing subjects:", observed()$r0)
        paste(str1, str2, str3, str4, str5, str6, str7, sep='<br/>')
      } else {
        paste(str1, str2, str3, str4, sep='<br/>')
      }
    }
  })



  output$inputdf <- renderDataTable(
    df(), options = list(pageLength = 10)
  )


  output$cum_accrual_plot <- renderPlotly({
    cum_accrual_plot <- observed()$cum_accrual_plot
    if (!is.null(cum_accrual_plot)) cum_accrual_plot
  })

  output$daily_accrual_plot <- renderPlotly({
    daily_accrual_plot <- observed()$daily_accrual_plot
    if (!is.null(daily_accrual_plot)) daily_accrual_plot
  })

  output$event_km_plot <- renderPlotly({
    event_km_plot <- observed()$event_km_plot
    if (!is.null(event_km_plot)) event_km_plot
  })

  output$dropout_km_plot <- renderPlotly({
    dropout_km_plot <- observed()$dropout_km_plot
    if (!is.null(dropout_km_plot)) dropout_km_plot
  })



  output$enroll_fit <- renderPlotly({
    if (!is.null(enroll_fit())) enroll_fit()$enroll_fit_plot
  })

  output$event_fit <- renderPlotly({
    if (!is.null(event_fit())) event_fit()$event_fit_plot
  })

  output$dropout_fit <- renderPlotly({
    if (!is.null(dropout_fit())) dropout_fit()$dropout_fit_plot
  })






  # enrollment and event predication date
  output$pred_date <- renderText({
    if (to_predict() == 'Enrollment only' ||
        to_predict() == 'Enrollment and event') {

      req(pred()$enroll_pred)


      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        if (!is.null(pred()$enroll_pred$enroll_pred_date)) {
          str1 <- paste0("Time from cutoff until ",
                         pred()$enroll_pred$target_n, " subjects: ",
                         pred()$enroll_pred$enroll_pred_date[1] -
                           observed()$cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ",
                         pred()$enroll_pred$enroll_pred_date[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$enroll_pred$enroll_pred_date[2], ", ",
                         pred()$enroll_pred$enroll_pred_date[3])
          text1 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text1 <- NULL
        }
      } else {
        if (!is.null(pred()$enroll_pred$enroll_pred_day)) {
          str1 <- paste0("Time from trial start until ",
                         pred()$enroll_pred$target_n, " subjects")
          str2 <- paste0("Median prediction day: ",
                         pred()$enroll_pred$enroll_pred_day[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$enroll_pred$enroll_pred_day[2], ", ",
                         pred()$enroll_pred$enroll_pred_day[3])
          text1 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text1 <- NULL
        }
      }
    }


    if (to_predict() == 'Enrollment and event' ||
        to_predict() == 'Event only') {

      req(pred()$event_pred)

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        if (input$dropout_model != "None") {
          shiny::validate(need(observed()$c0 > 0, paste(
            "The number of dropouts must be positive",
            "to fit a dropout model.")))
        }

        if (!is.null(pred()$event_pred$event_pred_date)) {
          str1 <- paste0("Time from cutoff until ",
                         pred()$event_pred$target_d, " events: ",
                         pred()$event_pred$event_pred_date[1] -
                           observed()$cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ",
                         pred()$event_pred$event_pred_date[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$event_pred$event_pred_date[2], ", ",
                         pred()$event_pred$event_pred_date[3])
          text2 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text2 <- NULL
        }
      } else {
        if (!is.null(pred()$event_pred$event_pred_day)) {
          str1 <- paste0("Time from trial start until ",
                         pred()$event_pred$target_d, " events")
          str2 <- paste0("Median prediction day: ",
                         pred()$event_pred$event_pred_day[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$event_pred$event_pred_day[2], ", ",
                         pred()$event_pred$event_pred_day[3])
          text2 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text2 <- NULL
        }
      }
    }


    if (to_predict() == 'Enrollment only') {
      text <- text1
    } else if (to_predict() == 'Event only') {
      text <- text2
    } else if (to_predict() == 'Enrollment and event') {
      text <- paste(text1, '', text2, sep='<br/>')
    }

    if (!is.null(text)) text

  })




  # enrollment prediction plot
  output$pred_plot <- renderPlotly({
    if (to_predict() == "Enrollment only") {
      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))
      }

      p1 <- pred()$enroll_pred$enroll_pred_plot
    } else {
      shiny::validate(need(
        showEnrollment() || showEvent() || showDropout() || showOngoing(),
        "Need at least one parameter to show on prediction plot"))

      req(pred()$event_pred)


      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        if (to_predict() == "Enrollment and event")
          shiny::validate(
            need(target_n() > observed()$n0,
                 "Target enrollment has been reached."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        if (input$dropout_model != "None") {
          shiny::validate(need(observed()$c0 > 0, paste(
            "The number of dropouts must be positive",
            "to fit a dropout model.")))
        }

      }


      dfs <- dplyr::tibble()

      if (showEnrollment()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$enroll_pred_df)
      }
      if (showEvent()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$event_pred_df)
      }
      if (showDropout()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$dropout_pred_df)
      }
      if (showOngoing()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)
      }

      dfs$parameter <- factor(dfs$parameter, levels = c(
        "Enrollment", "Event", "Dropout", "Ongoing"))


      if (input$stage != 'Design stage') {
        dfa <- dfs %>% dplyr::filter(is.na(lower))
        dfb <- dfs %>% dplyr::filter(!is.na(lower))

        n_months = lubridate::interval(min(dfs$date),
                                       max(dfs$date)) %/% months(1)
        bw = fbw(n_months)


        p1 <- plotly::plot_ly() %>%
          plotly::add_ribbons(
            data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", fillcolor = ~parameter,
            line = list(width=0)) %>%
          plotly::add_lines(
            data = dfb, x = ~date, y = ~n, color = ~parameter) %>%
          plotly::add_lines(
            data = dfa, x = ~date, y = ~n, color = ~parameter,
            line = list(shape="hv", width=2)) %>%
          plotly::add_lines(
            x = rep(observed()$cutoffdt, 2), y = range(dfs$n),
            name = "cutoff",
            line = list(dash="dash"),
            showlegend = FALSE) %>%
          plotly::layout(
            annotations = list(x = observed()$cutoffdt,
                               y = 0.1, yref = "paper",
                               text = 'cutoff', xanchor = "left",
                               font = list(size = 12),
                               showarrow = FALSE),
            xaxis = list(title = "",
                         zeroline = FALSE,
                         tickmode = "linear", dtick = bw),
            yaxis = list(zeroline = FALSE),
            legend = list(x = 0, y = 1.2, orientation = 'h'))

        if (showEvent()) {
          p1 <- p1 %>%
            plotly::add_lines(
              x = range(dfs$date), y = rep(target_d(), 2),
              name = 'target events', showlegend = FALSE,
              line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
            plotly::layout(
              annotations = list(x = 0.95, xref = "paper", y = target_d(),
                                 text = 'target events',
                                 xanchor = "right", yanchor = "bottom",
                                 font = list(size = 12),
                                 showarrow = FALSE))
        }


      } else {
        p1 <- plotly::plot_ly() %>%
          plotly::add_ribbons(
            data = dfs, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", fillcolor = ~parameter,
            line = list(width=0)) %>%
          plotly::add_lines(
            data = dfs, x = ~t, y = ~n, color = ~parameter) %>%
          plotly::layout(
            xaxis = list(title = "Days since trial start",
                         zeroline = FALSE),
            yaxis = list(zeroline = FALSE),
            legend = list(x = 0, y = 1.2, orientation = 'h'))

        if (showEvent()) {
          p1 <- p1 %>%
            plotly::add_lines(
              x = range(dfs$t), y = rep(target_d(), 2),
              name = 'target events', showlegend = FALSE,
              line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
            plotly::layout(
              annotations = list(x = 0.95, xref = "paper", y = target_d(),
                                 text = 'target events',
                                 xanchor = "right", yanchor = "bottom",
                                 font = list(size = 12),
                                 showarrow = FALSE))
        }
      }
    }

    p1

  })



  output$downloadSumdata <- downloadHandler(
    filename = function() {
      paste0("sumdata-", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      if (to_predict() == "Enrollment only") {
        sumdata <- pred()$enroll_pred$enroll_pred_df
      } else {
        sumdata <- pred()$event_pred$enroll_pred_df %>%
          dplyr::bind_rows(pred()$event_pred$event_pred_df) %>%
          dplyr::bind_rows(pred()$event_pred$dropout_pred_df) %>%
          dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)
      }
      openxlsx::write.xlsx(sumdata, file)
    }
  )



  observeEvent(input$add_accrualTime, {
    a = matrix(as.numeric(input$accrualTime),
               ncol=ncol(input$accrualTime))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$accrualTime)
    updateMatrixInput(session, "accrualTime", c)
  })


  observeEvent(input$del_accrualTime, {
    if (nrow(input$accrualTime) >= 2) {
      a = matrix(as.numeric(input$accrualTime),
                 ncol=ncol(input$accrualTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$accrualTime)
      updateMatrixInput(session, "accrualTime", b)
    }
  })


  observeEvent(input$add_piecewise_poisson_rate, {
    a = matrix(as.numeric(input$piecewise_poisson_rate),
               ncol=ncol(input$piecewise_poisson_rate))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$piecewise_poisson_rate)
    updateMatrixInput(session, "piecewise_poisson_rate", c)
  })


  observeEvent(input$del_piecewise_poisson_rate, {
    if (nrow(input$piecewise_poisson_rate) >= 2) {
      a = matrix(as.numeric(input$piecewise_poisson_rate),
                 ncol=ncol(input$piecewise_poisson_rate))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$piecewise_poisson_rate)
      updateMatrixInput(session, "piecewise_poisson_rate", b)
    }
  })


  observeEvent(input$add_piecewiseSurvivalTime, {
    a = matrix(as.numeric(input$piecewiseSurvivalTime),
               ncol=ncol(input$piecewiseSurvivalTime))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$piecewiseSurvivalTime)
    updateMatrixInput(session, "piecewiseSurvivalTime", c)
  })


  observeEvent(input$del_piecewiseSurvivalTime, {
    if (nrow(input$piecewiseSurvivalTime) >= 2) {
      a = matrix(as.numeric(input$piecewiseSurvivalTime),
                 ncol=ncol(input$piecewiseSurvivalTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$piecewiseSurvivalTime)
      updateMatrixInput(session, "piecewiseSurvivalTime", b)
    }
  })


  observeEvent(input$add_piecewise_exponential_survival, {
    a = matrix(as.numeric(input$piecewise_exponential_survival),
               ncol=ncol(input$piecewise_exponential_survival))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$piecewise_exponential_survival)
    updateMatrixInput(session, "piecewise_exponential_survival", c)
  })


  observeEvent(input$del_piecewise_exponential_survival, {
    if (nrow(input$piecewise_exponential_survival) >= 2) {
      a = matrix(as.numeric(input$piecewise_exponential_survival),
                 ncol=ncol(input$piecewise_exponential_survival))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$piecewise_exponential_survival)
      updateMatrixInput(session, "piecewise_exponential_survival", b)
    }
  })


  observeEvent(input$add_piecewiseDropoutTime, {
    a = matrix(as.numeric(input$piecewiseDropoutTime),
               ncol=ncol(input$piecewiseDropoutTime))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$piecewiseDropoutTime)
    updateMatrixInput(session, "piecewiseDropoutTime", c)
  })


  observeEvent(input$del_piecewiseDropoutTime, {
    if (nrow(input$piecewiseDropoutTime) >= 2) {
      a = matrix(as.numeric(input$piecewiseDropoutTime),
                 ncol=ncol(input$piecewiseDropoutTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$piecewiseDropoutTime)
      updateMatrixInput(session, "piecewiseDropoutTime", b)
    }
  })


  observeEvent(input$add_piecewise_exponential_dropout, {
    a = matrix(as.numeric(input$piecewise_exponential_dropout),
               ncol=ncol(input$piecewise_exponential_dropout))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1:nrow(c)))
    colnames(c) = colnames(input$piecewise_exponential_dropout)
    updateMatrixInput(session, "piecewise_exponential_dropout", c)
  })


  observeEvent(input$del_piecewise_exponential_dropout, {
    if (nrow(input$piecewise_exponential_dropout) >= 2) {
      a = matrix(as.numeric(input$piecewise_exponential_dropout),
                 ncol=ncol(input$piecewise_exponential_dropout))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1:nrow(b)))
      colnames(b) = colnames(input$piecewise_exponential_dropout)
      updateMatrixInput(session, "piecewise_exponential_dropout", b)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

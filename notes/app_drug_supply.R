library(shiny)
library(shinyMatrix)
library(shinyFeedback)
library(shinyjs)
library(shinybusy)
library(readxl)
library(writexl)
library(dplyr)
library(prompter)
library(plotly)
library(eventPred)


# conditional panels for number of treatments
f_treatment_allocation <- function(i) {
  conditionalPanel(
    condition = paste0("input.k == ", i),

    shinyMatrix::matrixInput(
      paste0("treatment_allocation_",i),
      label = tags$span(
        "Treatment allocation",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = "in a randomization block",
                     position = "right")
      ),

      value = matrix(rep(1,i),
                     ncol = 1,
                     dimnames = list(paste0("Treatment ", seq_len(i)),
                                     "Size")),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Exponential'",
                      "&&", "input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("exponential_survival_", i),
      "Hazard rate for each treatment",
      value = matrix(rep(0.0030, i),
                     nrow = 1,
                     dimnames = list(
                       NULL, paste("Treatment", 1:i))
      ),
      inputClass = "numeric",
      rows = list(names=FALSE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_weibull_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Weibull'",
                      "&&", "input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("weibull_survival_", i),
      "Weibull parameters",
      value = matrix(rep(c(1.42, 392), i),
                     nrow = 2,
                     byrow = FALSE,
                     dimnames = list(
                       c("Shape", "Scale"),
                       paste("Treatment", 1:i))
      ),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_lnorm_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Log-normal'",
                      "&&", "input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("lnorm_survival_", i),
      "Log-normal parameters",
      value = matrix(rep(c(5.4, 1), i),
                     nrow = 2,
                     byrow = FALSE,
                     dimnames = list(
                       c("Mean on log scale", "SD on log scale"),
                       paste("Treatment", 1:i))
      ),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}



f_piecewise_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste(
      "input.event_prior == 'Piecewise exponential'",
      "&&", "input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("piecewise_exponential_survival_", i),
      "Hazard rate by time interval for each treatment",
      value = matrix(c(0, rep(0.0030, i)),
                     nrow = 1,
                     dimnames = list(
                       "Interval 1",
                       c("Starting time", paste("Treatment", 1:i)))
      ),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    ),

    actionButton(paste0("add_piecewise_exponential_survival_", i),
                 label=NULL, icon=icon("plus")),
    actionButton(paste0("del_piecewise_exponential_survival_", i),
                 label=NULL, icon=icon("minus"))
  )
}



f_treatment_by_drug <- function(i, j) {
  conditionalPanel(
    condition = paste(
      "input.k ==", i, "&&", "input.l ==", j),

    shinyMatrix::matrixInput(
      paste0("treatment_by_drug_", i, "_", j),
      "Drugs contained in each treatment",
      value = matrix(1,
                     nrow = i,
                     ncol = j,
                     dimnames = list(
                       paste("Treatment", 1:i),
                       paste("Drug", 1:j))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_dosing_schedule <- function(j) {
  conditionalPanel(
    condition = paste("input.l ==", j),

    shinyMatrix::matrixInput(
      paste0("dosing_schedule_", j),
      "Dosing schedule for each drug",
      value = matrix(c(3, 10000, rep(100, j)),
                     nrow = 1,
                     ncol = j+2,
                     dimnames = list(
                       "Interval 1",
                       c("Frequency in weeks",
                         "Number of cycles",
                         paste("Drug", 1:j, "dose")))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    ),

    actionButton(paste0("add_dosing_schedule_", j),
                 label=NULL, icon=icon("plus")),
    actionButton(paste0("del_dosing_schedule_", j),
                 label=NULL, icon=icon("minus"))
  )
}


# cumulative dose for duration x and a drug with dosing schedule (w, N, d)
f_cum_dose <- function(x, w, N, d) {
  m = length(w)  # number of dosing intervals
  u = c(0, cumsum(7*w*N))
  i = findInterval(x, u)
  v = c(0, cumsum(N*d))
  n = floor((x - u[i])/(7*w[i])) + 1
  v[i] + n*d[i]
}


observedPanel <- tabPanel(
  title = "Observed Data",
  value = "observed_data_panel",


  htmlOutput("dates"),
  verbatimTextOutput("statistics"),

  plotlyOutput("cum_accrual_plot"),

  conditionalPanel(
    condition = "input.stage != 'Real-time after enrollment completion'",
    plotlyOutput("daily_accrual_plot")
  ),

  conditionalPanel(
    condition = "input.to_predict == 'Enrollment and event' ||
    input.stage == 'Real-time after enrollment completion'",

    plotlyOutput("event_km_plot")
  ),



  conditionalPanel(
    condition = "input.stage != 'Design stage'",
    dataTableOutput("input_df"))
)



enrollmentPanel <- tabPanel(
  title = "Enrollment Model",
  value = "enroll_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "enroll_prior",
        "Which enrollment model to use?",
        choices = c("Poisson",
                    "Time-decay",
                    "Piecewise Poisson"),
        selected = "Piecewise Poisson",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.enroll_prior == 'Poisson'",

               numericInput(
                 "poisson_rate",
                 "Daily enrollment rate",
                 value = 1,
                 min = 0, max = 100, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_prior == 'Time-decay'",

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
               condition = "input.enroll_prior ==
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
      column(4, radioButtons(
        "event_prior",
        "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Piecewise exponential",
        inline = FALSE)
      ),


      column(8,
             lapply(1:6, f_exponential_survival),
             lapply(1:6, f_weibull_survival),
             lapply(1:6, f_lnorm_survival),
             lapply(1:6, f_piecewise_exponential_survival)
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
                    "Model averaging",
                    "Spline"),
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
             ),

             conditionalPanel(
               condition = "input.event_model == 'Spline'",

               numericInput(
                 "spline_k",
                 "How many inner knots to use?",
                 value = 0,
                 min = 0, max = 10, step = 1),

               radioButtons(
                 "spline_scale",
                 "Which scale to model as a spline function?",
                 choices = c("hazard", "odds", "normal"),
                 selected = "hazard",
                 inline = TRUE)
             )
      )
    ),

    uiOutput("event_fit")
  )
)



dosingPanel <- tabPanel(
  title = "Dosing Model",
  value = "dosing_model_panel",

  selectInput(
    "l", "Number of drugs",
    choices = seq_len(12), selected = 2),

  lapply(c(outer(1:6, 1:12, function(i,j) 1000*i+j)),
         function(k) f_treatment_by_drug(k %/% 1000, k %% 1000)),

  lapply(1:12, f_dosing_schedule)

)



eventPredictPanel <- tabPanel(
  title = "Event Prediction",
  value = "event_prediction_panel",

  uiOutput("pred_date"),
  uiOutput("pred_plot"),

  downloadButton("downloadSumdata", "Download summary data"),
  downloadButton("downloadSimdata", "Download simulated data")
)


dosingPredictPanel <- tabPanel(
  title = "Dosing Prediction",
  value = "dosing_prediction_panel",

  uiOutput("dosing_plot"),
  downloadButton("downloadDosingdata", "Download dosing data")
)



# reduced style fileInput
fileInputNoExtra<-function(inputId, label, multiple = FALSE, accept = NULL,
                           width = NULL, buttonLabel = "Browse...",
                           placeholder = "No file selected"){

  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file",
                         style = "display: none;",
                         `data-restore` = restoredValue)
  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")

  tags$label(
    class = "input-group-btn",
    type="button",
    style=if (!is.null(width))
      paste0("width: ", validateCssUnit(width),";",
             "padding-right: 5px; padding-bottom: 0px; display:inline-block;"),

    span(class = "btn btn-default btn-file",type="button",
         buttonLabel, inputTag,
         style=if (!is.null(width))
           paste0("width: ", validateCssUnit(width),";",
                  "border-radius: 4px; padding-bottom:5px;"))
  )
}

# user interface ----------------
ui <- fluidPage(

  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  prompter::use_prompt(),

  add_busy_spinner(),

  titlePanel(tagList(
    span(HTML(paste(tags$span(style="font-size:14pt",
                              "Enrollment and Treatment Prediction"))),
         span(actionButton(
           "predict", "Predict",
           style="color: #fff; background-color: #337ab7;
          border-color: #2e6da4"),

           downloadButton("saveInputs", "Save inputs"),
           fileInputNoExtra("loadInputs", label=NULL, accept=".rds",
                            buttonLabel=list(icon("upload"), "Load inputs"),
                            width="116px"),
           tags$a(tags$span(icon(name = "question-circle")), target="_blank",
                  href="simplified_manual.pdf"),
           style="position:absolute;right:0.5em;",
           tags$style(type='text/css', "#saveInputs{margin-top: -5px;}")
         ))),
    windowTitle = "Enrollment and Treatment Prediction"),


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
        column(7,
               conditionalPanel(
                 condition = "input.stage == 'Design stage' ||
                 input.stage == 'Real-time before enrollment completion'",

                 numericInput(
                   "target_n",
                   "Target enrollment",
                   value = 300,
                   min = 1, max = 20000, step = 1)
               )
        ),

        column(5,
               conditionalPanel(
                 condition = "input.to_predict == 'Enrollment and event' ||
               input.stage == 'Real-time after enrollment completion'",

                 numericInput(
                   "target_d",
                   "Target events",
                   value = 200,
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
        column(7, radioButtons(
          "pilevel",
          "Prediction interval",

          choices = c("95%" = "0.95", "90%" = "0.90", "80%" = "0.80"),
          selected = "0.95",
          inline = TRUE)
        ),

        column(5, numericInput(
          "nyears",
          "Years after cutoff",
          value = 1,
          min = 1, max = 10, step = 1)
        )
      ),


      conditionalPanel(
        condition = "input.to_predict == 'Enrollment and event' ||
        input.stage == 'Real-time after enrollment completion'",

        checkboxGroupInput(
          "to_show",
          "What to show on prediction plot?",
          choices = c("Enrollment", "Event", "Ongoing"),
          selected = c("Enrollment", "Event"),
          inline = TRUE
        )
      ),



      fluidRow(
        column(7, checkboxInput(
          "by_treatment",
          "By treatment?",
          value = FALSE),

          conditionalPanel(
            condition = "input.by_treatment &&
            (input.to_predict == 'Enrollment and event' ||
            input.stage == 'Real-time after enrollment completion')",

            checkboxInput(
              "predict_dosing",
              "Predict dosing?",
              value = FALSE)
          )
        ),


        column(5, conditionalPanel(
          condition = "input.stage == 'Design stage' ||
                   (input.by_treatment &&
                   input.stage != 'Real-time after enrollment completion')",

          selectInput(
            "k", "Treatments",
            choices = seq_len(6), selected = 2))
        )

      ),


      conditionalPanel(
        condition = "input.stage == 'Design stage' ||
        (input.by_treatment &&
        input.stage != 'Real-time after enrollment completion')",

        lapply(2:6, f_treatment_allocation)
      ),


      fluidRow(
        column(7, numericInput(
          "nreps",
          label = "Simulation runs",
          value = 500,
          min = 100, max = 10000, step = 1)
        ),

        column(5, numericInput(
          "seed",
          label = "Seed",
          value = 2000,
          min = 0, max = 100000, step = 1
        ))
      )



    ),


    mainPanel(

      tabsetPanel(
        id = "results",
        observedPanel,
        enrollmentPanel,
        eventPanel,
        dosingPanel,
        eventPredictPanel,
        dosingPredictPanel
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


  # whether to show or hide enrollment and event panels
  observeEvent(to_predict(), {
    if (to_predict() == 'Enrollment only') {
      showTab(inputId = "results", target = "enroll_model_panel")
      hideTab(inputId = "results", target = "event_model_panel")
    } else if (to_predict() == 'Enrollment and event') {
      showTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
    } else if (to_predict() == 'Event only') {
      hideTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
    }
  })



  # whether to predict dosing
  predict_dosing <- reactive({
    if (input$by_treatment &&
        (input$to_predict == 'Enrollment and event' ||
         input$stage == 'Real-time after enrollment completion')) {
      input$predict_dosing
    } else {
      FALSE
    }
  })


  observeEvent(predict_dosing(), {
    if (predict_dosing()) {
      showTab(inputId = "results", target = "dosing_model_panel")
      showTab(inputId = "results", target = "dosing_prediction_panel")
    } else {
      hideTab(inputId = "results", target = "dosing_model_panel")
      hideTab(inputId = "results", target = "dosing_prediction_panel")
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

    if (input$to_predict == "Enrollment and event") {
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


  pilevel <- reactive(as.numeric(input$pilevel))


  showEnrollment <- reactive({
    "Enrollment" %in% input$to_show
  })


  showEvent <- reactive({
    "Event" %in% input$to_show
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


  k <- reactive(as.numeric(input$k))


  treatment_allocation <- reactive({
    req(k())
    if (k() > 1) {
      d = input[[paste0("treatment_allocation_", k())]]
      d <- as.numeric(d)

      valid = all(d > 0 & d == round(d))
      if (!valid) {
        showNotification("Treatment allocation must be positive integers")
      }
      req(valid)
      d
    } else {
      1
    }
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
    req(k())
    param = input[[paste0("exponential_survival_", k())]]
    lambda = as.numeric(param)
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
    req(k())
    param = input[[paste0("weibull_survival_", k())]]
    shape = as.numeric(param[1,])
    scale = as.numeric(param[2,])

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
    req(k())
    param = input[[paste0("lnorm_survival_", k())]]
    meanlog = as.numeric(param[1,])
    sdlog = as.numeric(param[2,])

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
    req(k())
    param = input[[paste0("piecewise_exponential_survival_", k())]]
    t = as.numeric(param[,1])
    lambda = as.numeric(param[,-1])

    valid1 = all(diff(t) > 0) && (t[1] == 0)
    if (!valid1) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }

    valid2 = all(lambda > 0)
    if (!valid2) {
      showNotification(
        "Hazard rate must be positive"
      )
    }

    req(valid1 && valid2)

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


  spline_k <- reactive({
    req(input$spline_k)
    valid = (input$spline_k >= 0 && input$spline_k == round(input$spline_k))
    shinyFeedback::feedbackWarning(
      "spline_k", !valid,
      "Number of inner knots must be a nonnegative integer")
    req(valid)
    as.numeric(input$spline_k)
  })




  l <- reactive(as.numeric(input$l))

  treatment_by_drug <- reactive({
    req(l())
    param = input[[paste0("treatment_by_drug_", k(), "_", l())]]
    t = as.numeric(param)

    valid = all(t==1 | t==0)
    if (!valid) {
      showNotification(
        "Entries of treatment by drug matrix must be 1 or 0"
      )
    }

    req(valid)

    dplyr::tibble(treatment = rep(1:k(), l()),
                  drug = rep(1:l(), each=k()),
                  included = as.logical(t)) %>%
      dplyr::filter(included) %>%
      dplyr::select(treatment, drug)
  })



  dosing_schedule <- reactive({
    req(l())
    param = input[[paste0("dosing_schedule_", l())]]
    w = as.numeric(param[,1])
    N = as.numeric(param[,2])
    d = as.numeric(param[,-c(1,2)])

    valid1 = all(w > 0)
    if (!valid1) {
      showNotification("Dosing frequency must be positive")
    }

    valid2 = all(N > 0 & N == round(N))
    if (!valid2) {
      showNotification("Number of cycles must be positive integers")
    }

    valid3 = all(d >= 0)
    if (!valid3) {
      showNotification("Doses must be nonnegative")
    }

    req(valid1 && valid2 && valid3)

    matrix(c(w, N, d), ncol = l() + 2)
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
      required_columns <- c('trialsdt', 'randdt', 'cutoffdt')
    } else {
      required_columns <- c('trialsdt', 'randdt', 'cutoffdt', 'time', 'event',
                            'dropout')
    }

    if (input$by_treatment) {
      required_columns <- c(required_columns, 'treatment')
    }

    column_names <- colnames(df)

    shiny::validate(
      need(all(required_columns %in% column_names),
           "You don't have the right data")
    )

    dplyr::tibble(df) %>%
      dplyr::mutate(trialsdt = as.Date(trialsdt),
                    randdt = as.Date(randdt),
                    cutoffdt = as.Date(cutoffdt))

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
    if (!is.null(df())) {
      if (!input$by_treatment) {
        event_fit <- fitEvent(
          df(), input$event_model, piecewiseSurvivalTime(),
          spline_k(), input$spline_scale, showplot = FALSE)
      } else {
        df_list <- split(df(), df()$treatment)

        event_fit <- lapply(df_list, function(df) fitEvent(
          df, input$event_model, piecewiseSurvivalTime(),
          spline_k(), input$spline_scale, showplot = FALSE))
      }

      event_fit
    }
  })




  # enrollment and event prediction
  pred <- eventReactive(input$predict, {
    set.seed(as.numeric(input$seed))

    if (to_predict() != "Enrollment only") {
      shiny::validate(
        need(showEnrollment() || showEvent() || showOngoing(),
             "Need at least one parameter to show on prediction plot"))
    }

    if (input$stage == "Design stage") {

      k = k()
      alloc = treatment_allocation()

      # enroll model specifications
      if (input$enroll_prior == "Poisson") {
        theta = log(poisson_rate())
      } else if (input$enroll_prior == "Time-decay") {
        theta = c(log(mu()), log(delta()))
      } else if (input$enroll_prior == "Piecewise Poisson") {
        theta = log(piecewise_poisson_rate()[,2])
        accrualTime = piecewise_poisson_rate()[,1]
      }

      enroll_prior <- list(
        model = input$enroll_prior,
        theta = theta,
        vtheta = diag(length(theta))*1e-8)

      if (input$enroll_prior == "Piecewise Poisson") {
        enroll_prior$accrualTime = accrualTime
      }


      # event model specifications
      if (to_predict() == "Enrollment and event" ||
          to_predict() == "Event only") {
        if (input$event_prior == "Exponential") {
          theta = log(exponential_survival())
        } else if (input$event_prior == "Weibull") {
          theta = c(log(weibull_survival()))
        } else if (input$event_prior == "Log-normal") {
          theta = c(matrix(c(lnorm_survival()[1,],
                             log(lnorm_survival()[2,])),
                           nrow = 2, byrow = TRUE))
        } else if (input$event_prior == "Piecewise exponential") {
          theta = c(log(piecewise_exponential_survival()[,-1]))
          piecewiseSurvivalTime = piecewise_exponential_survival()[,1]
        }

        event_prior = list(
          model = input$event_prior,
          ngroups = k,
          alloc = alloc,
          theta = theta,
          vtheta = diag(length(theta))*1e-8)

        if (input$event_prior == "Piecewise exponential") {
          event_prior$piecewiseSurvivalTime = piecewiseSurvivalTime
        }
      }


      # get prediction results based on what to predict
      if (!input$by_treatment || k == 1) {
        if (to_predict() == "Enrollment only") {
          getPrediction(
            to_predict = to_predict(),
            target_n = target_n(),
            enroll_prior = enroll_prior,
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showsummary = FALSE,
            showplot = FALSE)
        } else if (to_predict() == "Enrollment and event") {
          getPrediction(
            to_predict = to_predict(),
            target_n = target_n(),
            target_d = target_d(),
            enroll_prior = enroll_prior,
            event_prior = event_prior,
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showEnrollment = showEnrollment(),
            showEvent = showEvent(),
            showOngoing = showOngoing(),
            showsummary = FALSE,
            showplot = FALSE)
        }
      } else {  # by treatment

        if (to_predict() == "Enrollment only") {
          pred1 <- getPrediction(
            to_predict = to_predict(),
            target_n = target_n(),
            enroll_prior = enroll_prior,
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showsummary = FALSE,
            showplot = FALSE)

          n = target_n()
          nreps = nreps()
          blocksize = sum(alloc)
          nblocks = ceiling(n/blocksize)
          m = nblocks*blocksize
          treats = rep(1:k, alloc)
          index = rep(1:n, nreps) + rep((0:(nreps-1))*m, each=n)
          treatment = c(replicate(nreps*nblocks, sample(treats)))[index]

          newSubjects <- pred1$enroll_pred$newSubjects %>%
            dplyr::bind_cols(treatment=treatment)

          plower = (1 - pilevel())/2
          pupper = 1 - plower

          t = pred1$enroll_pred$enroll_pred_df$t

          # predicted number of subjects enrolled by treatment
          dfs <- dplyr::tibble(t = t) %>%
            dplyr::cross_join(newSubjects) %>%
            dplyr::group_by(treatment, t, draw) %>%
            dplyr::summarise(nenrolled = sum(arrivalTime <= t),
                             .groups = "drop_last") %>%
            dplyr::summarise(n = quantile(nenrolled, probs = 0.5),
                             lower = quantile(nenrolled, probs = plower),
                             upper = quantile(nenrolled, probs = pupper),
                             mean = mean(nenrolled),
                             var = var(nenrolled),
                             .groups = "drop_last")

          enroll_pred_df <- dfs %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, pred1$enroll_pred$enroll_pred_df)) %>%
            dplyr::arrange(treatment, t)

          enroll_pred_plot <- list()

          enroll_pred_plot[[1]] <- pred1$enroll_pred$enroll_pred_plot %>%
            plotly::layout(annotations = list(
              x = 0.5, y = 1, text = "<b>Overall</b>",
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))

          for (i in 1:k) {
            enroll_pred_plot[[i+1]] <- dfs %>%
              dplyr::filter(treatment == i) %>%
              plotly::plot_ly(x = ~t) %>%
              plotly::add_ribbons(ymin = ~lower, ymax = ~upper,
                                  name = "prediction interval",
                                  fill = "tonexty",
                                  line = list(width=0)) %>%
              plotly::add_lines(y = ~n, name = "median prediction",
                                line = list(width=1)) %>%
              plotly::layout(xaxis = list(title = "Days since trial start",
                                          zeroline = FALSE),
                             yaxis = list(title = "Subjects",
                                          zeroline = FALSE),
                             legend = list(x = 0, y = 1.05,
                                           yanchor = "bottom",
                                           orientation = "h"),
                             annotations = list(
                               x = 0.5, y = 1,
                               text = paste0("<b>treatment=", i, "</b>"),
                               xanchor = "center", yanchor = "bottom",
                               showarrow = FALSE, xref='paper', yref='paper'))
          }

          pred <- pred1
          pred$enroll_pred$newSubjects <- newSubjects
          pred$enroll_pred$enroll_pred_df <- enroll_pred_df
          pred$enroll_pred$enroll_pred_plot <- enroll_pred_plot
          pred

        } else if (to_predict() == "Enrollment and event") {

          pred1 <- getPrediction(
            to_predict = to_predict(),
            target_n = target_n(),
            target_d = target_d(),
            enroll_prior = enroll_prior,
            event_prior = event_prior,
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showEnrollment = showEnrollment(),
            showEvent = showEvent(),
            showOngoing = showOngoing(),
            showsummary = FALSE,
            showplot = FALSE)

          newEvents <- pred1$event_pred$newEvents

          newSubjects <- newEvents %>%
            dplyr::select(draw, arrivalTime, treatment)

          plower = (1 - pilevel())/2
          pupper = 1 - plower

          t = pred1$enroll_pred$enroll_pred_df$t

          # predicted number of subjects enrolled by treatment
          enroll_pred_df <- dplyr::tibble(t = t) %>%
            dplyr::cross_join(newSubjects) %>%
            dplyr::group_by(treatment, t, draw) %>%
            dplyr::summarise(nenrolled = sum(arrivalTime <= t),
                             .groups = "drop_last") %>%
            dplyr::summarise(n = quantile(nenrolled, probs = 0.5),
                             lower = quantile(nenrolled, probs = plower),
                             upper = quantile(nenrolled, probs = pupper),
                             mean = mean(nenrolled),
                             var = var(nenrolled),
                             .groups = "drop_last") %>%
            dplyr::mutate(parameter = "Enrollment") %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, pred1$event_pred$enroll_pred_df)) %>%
            dplyr::arrange(treatment, t)


          # number of events and subjects at-risk
          t = pred1$event_pred$event_pred_df$t

          df1 = dplyr::tibble(t = t) %>%
            dplyr::cross_join(newEvents) %>%
            dplyr::group_by(treatment, t, draw) %>%
            dplyr::summarise(nevents = sum(totalTime <= t & event == 1),
                             natrisk = sum(arrivalTime <= t & totalTime > t),
                             .groups = "drop_last")

          # predicted number of events after data cut
          event_pred_df <- df1 %>%
            dplyr::group_by(treatment, t) %>%
            dplyr::summarise(n = quantile(nevents, probs = 0.5),
                             lower = quantile(nevents, probs = plower),
                             upper = quantile(nevents, probs = pupper),
                             mean = mean(nevents),
                             var = var(nevents),
                             .groups = "drop_last") %>%
            dplyr::mutate(parameter = "Event") %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, pred1$event_pred$event_pred_df)) %>%
            dplyr::arrange(treatment, t)


          # predicted number of subjects at risk after data cut
          ongoing_pred_df <- df1 %>%
            dplyr::group_by(treatment, t) %>%
            dplyr::summarise(n = quantile(natrisk, probs = 0.5),
                             lower = quantile(natrisk, probs = plower),
                             upper = quantile(natrisk, probs = pupper),
                             mean = mean(natrisk),
                             var = var(natrisk),
                             .groups = "drop_last") %>%
            dplyr::mutate(parameter = "Ongoing") %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, pred1$event_pred$ongoing_pred_df)) %>%
            dplyr::arrange(treatment, t)


          pred <- pred1
          pred$event_pred$enroll_pred_df <- enroll_pred_df
          pred$event_pred$event_pred_df <- event_pred_df
          pred$event_pred$ongoing_pred_df <- ongoing_pred_df
          pred

        }


      }
    } else { # real-time prediction
      shiny::validate(
        need(!is.null(df()),
             "Please upload data for real-time prediction."))

      if (input$by_treatment) {
        k = length(unique(df()$treatment))
      } else {
        k = 1
      }

      if (!input$by_treatment || k == 1) {
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
            pilevel = pilevel(),
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


          getPrediction(
            df = df(),
            to_predict = to_predict(),
            target_n = target_n(),
            target_d = target_d(),
            enroll_model = input$enroll_model,
            nknots = nknots(),
            lags = lags(),
            accrualTime = accrualTime(),
            event_model = input$event_model,
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            k = spline_k(),
            scale = input$spline_scale,
            dropout_model = "none",
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showEnrollment = showEnrollment(),
            showEvent = showEvent(),
            showOngoing = showOngoing(),
            showsummary = FALSE,
            showplot = FALSE)
        } else if (to_predict() == "Event only") {
          shiny::validate(
            need(target_d() > observed()$d0,
                 "Target number of events has been reached."))

          getPrediction(
            df = df(),
            to_predict = to_predict(),
            target_d = target_d(),
            event_model = input$event_model,
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            k = spline_k(),
            scale = input$spline_scale,
            dropout_model = "none",
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showEnrollment = showEnrollment(),
            showEvent = showEvent(),
            showOngoing = showOngoing(),
            showsummary = FALSE,
            showplot = FALSE)
        }
      } else { # by treatment
        if (to_predict() == "Enrollment only") {
          shiny::validate(
            need(target_n() > observed()$n0,
                 "Target enrollment has been reached."))

          shiny::validate(
            need(k() == k,
                 "Number of treatments must match observed."))

          alloc = treatment_allocation()

          pred1 <- getPrediction(
            df = df(),
            to_predict = to_predict(),
            target_n = target_n(),
            enroll_model = input$enroll_model,
            nknots = nknots(),
            lags = lags(),
            accrualTime = accrualTime(),
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showsummary = FALSE,
            showplot = FALSE)

          trialsdt = pred1$observed$trialsdt
          cutoffdt = pred1$observed$cutoffdt
          t0 = pred1$observed$t0

          n = target_n() - pred1$observed$n0
          nreps = nreps()
          blocksize = sum(alloc)
          nblocks = ceiling(n/blocksize)
          m = nblocks*blocksize
          treats = rep(1:k, alloc)
          index = rep(1:n, nreps) + rep((0:(nreps-1))*m, each=n)
          treatment = c(replicate(nreps*nblocks, sample(treats)))[index]

          newSubjects <- pred1$enroll_pred$newSubjects %>%
            dplyr::bind_cols(treatment=treatment)

          # number of subjects already enrolled by treatment
          sum_by_trt <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(n0 = n())


          plower = (1 - pilevel())/2
          pupper = 1 - plower

          t = pred1$enroll_pred$enroll_pred_df$t
          t = unique(t[t >= pred1$observed$t0])

          # predicted number of subjects enrolled by treatment after cutoff
          dfb <- dplyr::tibble(t = t) %>%
            dplyr::cross_join(newSubjects) %>%
            dplyr::group_by(treatment, t, draw) %>%
            dplyr::summarise(nenrolled = sum(arrivalTime <= t),
                             .groups = "drop_last") %>%
            dplyr::left_join(sum_by_trt, by = "treatment") %>%
            dplyr::mutate(nenrolled = nenrolled + n0) %>%
            dplyr::summarise(n = quantile(nenrolled, probs = 0.5),
                             lower = quantile(nenrolled, probs = plower),
                             upper = quantile(nenrolled, probs = pupper),
                             mean = mean(nenrolled),
                             var = var(nenrolled),
                             .groups = "drop_last")

          # arrival time for subjects already enrolled before data cut
          dfa <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::arrange(randdt) %>%
            dplyr::mutate(t = as.numeric(randdt - trialsdt + 1),
                          n = dplyr::row_number()) %>%
            dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0) %>%
            dplyr::select(treatment, t, n, lower, upper, mean, var) %>%
            dplyr::group_by(treatment, t) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::group_by(treatment)

          # extend observed to cutoff date
          dfa1 <- dfa %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::mutate(t = t0)

          if (max(dfa$t) < t0) {
            dfa <- dfa %>%
              dplyr::bind_rows(dfa1)
          }


          # concatenate subjects enrolled before and after data cut
          # add day 1
          df0 <- dplyr::tibble(treatment = 1:k,
                               t = 1, n = 0, lower = NA, upper = NA,
                               mean = 0, var = 0)

          # concatenate subjects enrolled before and after data cut
          dfs <- df0 %>%
            dplyr::bind_rows(dfa) %>%
            dplyr::bind_rows(dfb) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt))

          # separate data into observed and predicted
          dfa <- dfs %>% dplyr::filter(is.na(lower))
          dfb <- dfs %>% dplyr::filter(!is.na(lower))

          enroll_pred_df <- dfs %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, pred1$enroll_pred$enroll_pred_df)) %>%
            dplyr::arrange(treatment, t)

          enroll_pred_plot <- list()

          enroll_pred_plot[[1]] <- pred1$enroll_pred$enroll_pred_plot %>%
            plotly::layout(annotations = list(
              x = 0.5, y = 1, text = "<b>Overall</b>",
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))

          for (i in 1:k) {
            dfsi <- dfs %>%
              dplyr::filter(treatment == i)

            dfbi <- dfb %>%
              dplyr::filter(treatment == i)

            dfai <- dfa %>%
              dplyr::filter(treatment == i)

            enroll_pred_plot[[i+1]] <- plotly::plot_ly() %>%
              plotly::add_ribbons(data = dfbi, x = ~date,
                                  ymin = ~lower, ymax = ~upper,
                                  name = "prediction interval",
                                  fill = "tonexty",
                                  line = list(width=0)) %>%
              plotly::add_lines(data = dfbi, x = ~date, y = ~n,
                                name = "median prediction",
                                line = list(width=1)) %>%
              plotly::add_lines(data = dfai, x = ~date, y = ~n,
                                name = "observed",
                                line = list(shape="hv", width=1)) %>%
              plotly::add_lines(x = rep(cutoffdt, 2), y = range(dfsi$n),
                                name = "cutoff",
                                line = list(dash="dash"),
                                showlegend = FALSE) %>%
              plotly::layout(
                xaxis = list(title = "", zeroline = FALSE),
                yaxis = list(title = "Subjects", zeroline = FALSE),
                legend = list(x = 0, y = 1.05, yanchor = "bottom",
                              orientation = "h"),
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>treatment=", i, "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))
          }

          pred <- pred1
          pred$enroll_pred$newSubjects <- newSubjects
          pred$enroll_pred$enroll_pred_df <- enroll_pred_df
          pred$enroll_pred$enroll_pred_plot <- enroll_pred_plot
          pred

        } else if (to_predict() == "Enrollment and event") {
          alloc = treatment_allocation()

          observed <- observed()
          trialsdt = observed$trialsdt
          cutoffdt = observed$cutoffdt
          t0 = observed$t0
          n0 = observed$n0
          d0 = observed$d0
          r0 = observed$r0

          df <- df() %>%
            dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                          totalTime = arrivalTime + time)

          shiny::validate(
            need(target_n() > n0,
                 "Target enrollment has been reached."))

          shiny::validate(
            need(target_d() > d0,
                 "Target number of events has been reached."))

          sum_by_trt <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(n0 = n(),
                             d0 = sum(event),
                             r0 = sum(!event))

          shiny::validate(
            need(k() == k,
                 "Number of treatments must match observed."))

          # predict enrollment
          enroll_fit <- fitEnrollment(
            df = observed$adsl,
            enroll_model = input$enroll_model,
            nknots = nknots(),
            accrualTime = accrualTime(),
            showplot = FALSE)

          enroll_pred <- predictEnrollment(
            df = observed$adsl,
            target_n = target_n(),
            enroll_fit = enroll_fit$enroll_fit,
            lags = lags(),
            pilevel = pilevel(),
            nyears = nyears(),
            nreps = nreps(),
            showsummary = FALSE,
            showplot = FALSE)



          # assign treatment for new subjects
          n = target_n() - n0
          nreps = nreps()
          blocksize = sum(alloc)
          nblocks = ceiling(n/blocksize)
          m = nblocks*blocksize
          treats = rep(1:k, alloc)
          index = rep(1:n, nreps) + rep((0:(nreps-1))*m, each=n)
          treatment = c(replicate(nreps*nblocks, sample(treats)))[index]

          newSubjects <- enroll_pred$newSubjects %>%
            dplyr::bind_cols(treatment=treatment)


          # enrollment by treatment
          plower = (1 - pilevel())/2
          pupper = 1 - plower

          t = enroll_pred$enroll_pred_df$t
          t = unique(t[t >= t0])

          # predicted number of subjects enrolled by treatment after cutoff
          dfb <- dplyr::tibble(t = t) %>%
            dplyr::cross_join(newSubjects) %>%
            dplyr::group_by(treatment, t, draw) %>%
            dplyr::summarise(nenrolled = sum(arrivalTime <= t),
                             .groups = "drop_last") %>%
            dplyr::left_join(sum_by_trt, by = "treatment") %>%
            dplyr::mutate(nenrolled = nenrolled + n0) %>%
            dplyr::summarise(n = quantile(nenrolled, probs = 0.5),
                             lower = quantile(nenrolled, probs = plower),
                             upper = quantile(nenrolled, probs = pupper),
                             mean = mean(nenrolled),
                             var = var(nenrolled),
                             .groups = "drop_last")

          # arrival time for subjects already enrolled before data cut
          dfa <- df %>%
            dplyr::group_by(treatment) %>%
            dplyr::arrange(randdt) %>%
            dplyr::mutate(t = as.numeric(randdt - trialsdt + 1),
                          n = dplyr::row_number()) %>%
            dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0) %>%
            dplyr::select(treatment, t, n, lower, upper, mean, var) %>%
            dplyr::group_by(treatment, t) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::group_by(treatment)

          # extend observed to cutoff date
          dfa1 <- dfa %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::mutate(t = t0)

          if (max(dfa$t) < t0) {
            dfa <- dfa %>%
              dplyr::bind_rows(dfa1)
          }

          # add day 1
          df0 <- dplyr::tibble(treatment = 1:k,
                               t = 1, n = 0, lower = NA, upper = NA,
                               mean = 0, var = 0)

          # concatenate subjects enrolled before and after data cut
          enroll_pred_df <- df0 %>%
            dplyr::bind_rows(dfa) %>%
            dplyr::bind_rows(dfb) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(parameter = "Enrollment") %>%
            dplyr::bind_rows(dplyr::bind_cols(
              treatment = 9999, parameter = "Enrollment",
              enroll_pred$enroll_pred_df)) %>%
            dplyr::arrange(treatment, t)


          # fit event by treatment
          event_fit <- list()
          for (i in 1:k) {
            event_fit[[i]] <- fitEvent(
              df = dplyr::filter(observed$adtte, treatment == i),
              event_model = input$event_model,
              piecewiseSurvivalTime = piecewiseSurvivalTime(),
              k = spline_k(), scale = input$spline_scale,
              showplot = FALSE)
          }


          # predict event by treatment
          event_pred <- list()
          for (i in 1:k) {
            event_pred[[i]] <- predictEvent(
              df = dplyr::filter(observed$adtte, treatment == i),
              target_d = target_d(),
              newSubjects = dplyr::filter(newSubjects, treatment == i),
              event_fit = event_fit[[i]]$event_fit,
              pilevel = pilevel(),
              nyears = nyears(),
              nreps = nreps(),
              showEnrollment = showEnrollment(),
              showEvent = showEvent(),
              showOngoing = showOngoing(),
              showsummary = FALSE,
              showplot = FALSE)
          }


          newEvents <- dplyr::tibble()
          for (i in 1:k) {
            newEvents <- newEvents %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$newEvents))
          }


          # survival distribution function for cumulative number of events
          sdf <- function(t, target_d, d0, newEvents) {
            sumdata <- newEvents %>%
              dplyr::group_by(draw) %>%
              dplyr::summarize(n = sum(totalTime <= t & event == 1) + d0)
            mean(sumdata$n < target_d)
          }

          # obtain the quantiles
          q = 1 - c(0.5, plower, pupper)
          pred_day = rep(NA, length(q))
          tmax = max(newEvents$totalTime[newEvents$event==1])
          for (j in 1:length(q)) {
            # check if the quantile can be estimated from observed data
            if (sdf(tmax, target_d(), d0, newEvents) <= q[j]) {
              pred_day[j] = uniroot(function(x)
                sdf(x, target_d(), d0, newEvents) - q[j],
                c(t0, tmax), tol = 1)$root
              pred_day[j] = ceiling(pred_day[j])
            }
          }

          pred_date <- as.Date(pred_day - 1, origin = trialsdt)

          str1 <- paste0("Time from cutoff until ", target_d(), " events: ",
                         pred_date[1] - cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ", pred_date[1])
          str3 <- paste0("Prediction interval: ", pred_date[2], ", ",
                         pred_date[3])
          s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")


          # set up time points for plotting
          t1 = t0 + 365*nyears()
          t = unique(c(seq(t0, t1, 30), t1))

          # number of events and subjects at-risk after data cut
          df1 = dplyr::tibble(t = t) %>%
            dplyr::cross_join(newEvents) %>%
            dplyr::group_by(t, draw) %>%
            dplyr::summarise(nevents = sum(totalTime <= t &
                                             event == 1) + d0,
                             natrisk = sum(arrivalTime <= t &
                                             totalTime > t),
                             .groups = "drop_last")

          # predicted number of events after data cut
          dfb = df1 %>%
            dplyr::summarise(n = quantile(nevents, probs = 0.5),
                             lower = quantile(nevents, probs = plower),
                             upper = quantile(nevents, probs = pupper),
                             mean = mean(nevents),
                             var = var(nevents))

          # predicted number of subjects at risk after data cut
          dff = df1 %>%
            dplyr::summarise(n = quantile(natrisk, probs = 0.5),
                             lower = quantile(natrisk, probs = plower),
                             upper = quantile(natrisk, probs = pupper),
                             mean = mean(natrisk),
                             var = var(natrisk))

          # observed number of events before data cut
          dfa <- df %>%
            dplyr::arrange(totalTime) %>%
            dplyr::mutate(t = totalTime,
                          n = cumsum(event),
                          lower = NA,
                          upper = NA) %>%
            dplyr::mutate(mean = n, var = 0) %>%
            dplyr::select(t, n, lower, upper, mean, var) %>%
            dplyr::group_by(t) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::ungroup()



          # observed number of subjects at risk before cutoff
          t2 = setdiff(sort(unique(c(df$arrivalTime, df$totalTime))), t0)

          dfe <- dplyr::tibble(t = t2) %>%
            dplyr::cross_join(df) %>%
            dplyr::group_by(t) %>%
            dplyr::summarise(n = sum(arrivalTime <= t &
                                       totalTime > t),
                             .groups = "drop_last") %>%
            dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0) %>%
            dplyr::select(t, n, lower, upper, mean, var) %>%
            dplyr::bind_rows(dplyr::tibble(t = t0, n = r0,
                                           lower = NA, upper = NA,
                                           mean = r0, var = 0))

          # day 1
          df0 <- dplyr::tibble(t = 1, n = 0, lower = NA, upper = NA,
                               mean = 0, var = 0)

          # add day 1 and concatenate events before and after data cut
          event_pred_df1 <- df0 %>%
            dplyr::bind_rows(dfa) %>%
            dplyr::bind_rows(dfb) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(treatment = 9999, parameter = "Event")

          # add day 1 and concatenate ongoing before and after data cut
          ongoing_pred_df1 <- df0 %>%
            dplyr::bind_rows(dfe) %>%
            dplyr::bind_rows(dff) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(treatment = 9999, parameter = "Ongoing")

          # add by treatment summary
          event_pred_df <- tibble()
          ongoing_pred_df <- tibble()
          for (i in 1:k) {
            event_pred_df <- event_pred_df %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$event_pred_df))

            ongoing_pred_df <- ongoing_pred_df %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$ongoing_pred_df))
          }

          event_pred_df <- event_pred_df %>%
            dplyr::bind_rows(event_pred_df1) %>%
            dplyr::arrange(treatment, t)

          ongoing_pred_df <- ongoing_pred_df %>%
            dplyr::bind_rows(ongoing_pred_df1) %>%
            dplyr::arrange(treatment, t)


          pred <- list(stage = "Real-time before enrollment completion",
                       to_predict = "Enrollment and event",
                       observed = observed,
                       enroll_fit = enroll_fit, enroll_pred = enroll_pred,
                       event_fit = event_fit,
                       event_pred = event_pred[[1]])

          pred$enroll_pred$newSubjects <- newSubjects
          pred$event_pred$newEvents <- newEvents
          pred$event_pred$enroll_pred_df <- enroll_pred_df
          pred$event_pred$event_pred_df <- event_pred_df
          pred$event_pred$ongoing_pred_df <- ongoing_pred_df
          pred$event_pred$event_pred_day <- pred_day
          pred$event_pred$event_pred_date <- pred_date
          pred$event_pred$event_pred_summary <- s1

          pred
        } else if (to_predict() == "Event only") {
          shiny::validate(
            need(target_d() > observed()$d0,
                 "Target number of events has been reached."))

          sum_by_trt <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(n0 = n(),
                             d0 = sum(event),
                             r0 = sum(!event))


          observed <- observed()
          trialsdt = observed$trialsdt
          cutoffdt = observed$cutoffdt
          t0 = observed$t0
          n0 = observed$n0
          d0 = observed$d0
          r0 = observed$r0

          df <- df() %>%
            dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                          totalTime = arrivalTime + time)

          # enrollment by treatment and overall
          dfpooled <- df %>% dplyr::mutate(treatment = 9999)

          dfa <- df %>%
            dplyr::bind_rows(dfpooled) %>%
            dplyr::group_by(treatment) %>%
            dplyr::arrange(randdt) %>%
            dplyr::mutate(t = as.numeric(randdt - trialsdt + 1),
                          n = dplyr::row_number()) %>%
            dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0) %>%
            dplyr::select(treatment, t, n, lower, upper, mean, var) %>%
            dplyr::group_by(treatment, t) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::group_by(treatment)

          # extend observed to cutoff date
          dfa1 <- dfa %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::mutate(t = t0)

          if (max(dfa$t) < t0) {
            dfa <- dfa %>%
              dplyr::bind_rows(dfa1)
          }

          # add predicted from data cut to specified years after data cut
          dfb1 <- dfa1 %>%
            dplyr::mutate(lower = n, upper = n, mean = n, var = 0)

          dfb2 <- dfb1 %>%
            dplyr::mutate(t = t0 + 365*nyears())

          # add day 1
          df0 <- dplyr::tibble(treatment = c(1:k, 9999),
                               t = 1, n = 0, lower = NA, upper = NA,
                               mean = 0, var = 0)

          enroll_pred_df <- df0 %>%
            dplyr::bind_rows(dfa) %>%
            dplyr::bind_rows(dfb1) %>%
            dplyr::bind_rows(dfb2) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(parameter = "Enrollment") %>%
            dplyr::arrange(treatment, t)


          # fit event by treatment
          event_fit <- list()
          for (i in 1:k) {
            event_fit[[i]] <- fitEvent(
              df = dplyr::filter(observed$adtte, treatment == i),
              event_model = input$event_model,
              piecewiseSurvivalTime = piecewiseSurvivalTime(),
              k = spline_k(), scale = input$spline_scale,
              showplot = FALSE)
          }


          # predict event by treatment
          event_pred <- list()
          for (i in 1:k) {
            event_pred[[i]] <- predictEvent(
              df = dplyr::filter(observed$adtte, treatment == i),
              target_d = target_d(),
              newSubjects = NULL,
              event_fit = event_fit[[i]]$event_fit,
              pilevel = pilevel(),
              nyears = nyears(),
              nreps = nreps(),
              showEnrollment = showEnrollment(),
              showEvent = showEvent(),
              showOngoing = showOngoing(),
              showsummary = FALSE,
              showplot = FALSE)
          }


          newEvents <- dplyr::tibble()
          for (i in 1:k) {
            newEvents <- newEvents %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$newEvents))
          }


          # survival distribution function for cumulative number of events
          sdf <- function(t, target_d, d0, newEvents) {
            sumdata <- newEvents %>%
              dplyr::group_by(draw) %>%
              dplyr::summarize(n = sum(totalTime <= t & event == 1) + d0)
            mean(sumdata$n < target_d)
          }


          plower = (1 - pilevel())/2
          pupper = 1 - plower

          # obtain the quantiles
          q = 1 - c(0.5, plower, pupper)
          pred_day = rep(NA, length(q))
          tmax = max(newEvents$totalTime[newEvents$event==1])
          for (j in 1:length(q)) {
            # check if the quantile can be estimated from observed data
            if (sdf(tmax, target_d(), d0, newEvents) <= q[j]) {
              pred_day[j] = uniroot(function(x)
                sdf(x, target_d(), d0, newEvents) - q[j],
                c(t0, tmax), tol = 1)$root
              pred_day[j] = ceiling(pred_day[j])
            }
          }

          pred_date <- as.Date(pred_day - 1, origin = trialsdt)

          str1 <- paste0("Time from cutoff until ", target_d(), " events: ",
                         pred_date[1] - cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ", pred_date[1])
          str3 <- paste0("Prediction interval: ", pred_date[2], ", ",
                         pred_date[3])
          s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")


          # set up time points for plotting
          t1 = t0 + 365*nyears()
          t = unique(c(seq(t0, t1, 30), t1))

          # number of events and subjects at-risk after data cut
          df1 = dplyr::tibble(t = t) %>%
            dplyr::cross_join(newEvents) %>%
            dplyr::group_by(t, draw) %>%
            dplyr::summarise(nevents = sum(totalTime <= t &
                                             event == 1) + d0,
                             natrisk = sum(arrivalTime <= t &
                                             totalTime > t),
                             .groups = "drop_last")

          # predicted number of events after data cut
          dfb = df1 %>%
            dplyr::summarise(n = quantile(nevents, probs = 0.5),
                             lower = quantile(nevents, probs = plower),
                             upper = quantile(nevents, probs = pupper),
                             mean = mean(nevents),
                             var = var(nevents))

          # predicted number of subjects at risk after data cut
          dff = df1 %>%
            dplyr::summarise(n = quantile(natrisk, probs = 0.5),
                             lower = quantile(natrisk, probs = plower),
                             upper = quantile(natrisk, probs = pupper),
                             mean = mean(natrisk),
                             var = var(natrisk))

          # observed number of events before data cut
          dfa <- df %>%
            dplyr::arrange(totalTime) %>%
            dplyr::mutate(t = totalTime,
                          n = cumsum(event),
                          lower = NA,
                          upper = NA) %>%
            dplyr::mutate(mean = n, var = 0) %>%
            dplyr::select(t, n, lower, upper, mean, var) %>%
            dplyr::group_by(t) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::ungroup()

          # observed number of subjects at risk before cutoff
          t2 = setdiff(sort(unique(c(df$arrivalTime, df$totalTime))), t0)

          dfe <- dplyr::tibble(t = t2) %>%
            dplyr::cross_join(df) %>%
            dplyr::group_by(t) %>%
            dplyr::summarise(n = sum(arrivalTime <= t &
                                       totalTime > t),
                             .groups = "drop_last") %>%
            dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0) %>%
            dplyr::select(t, n, lower, upper, mean, var) %>%
            dplyr::bind_rows(dplyr::tibble(t = t0, n = r0,
                                           lower = NA, upper = NA,
                                           mean = r0, var = 0))

          # day 1
          df0 <- dplyr::tibble(t = 1, n = 0, lower = NA, upper = NA,
                               mean = 0, var = 0)


          # add day 1 and concatenate events before and after data cut
          event_pred_df1 <- df0 %>%
            dplyr::bind_rows(dfa) %>%
            dplyr::bind_rows(dfb) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(treatment = 9999, parameter = "Event")

          # add day 1 and concatenate ongoing before and after data cut
          ongoing_pred_df1 <- df0 %>%
            dplyr::bind_rows(dfe) %>%
            dplyr::bind_rows(dff) %>%
            dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
            dplyr::mutate(treatment = 9999, parameter = "Ongoing")

          # add by treatment summary
          event_pred_df <- tibble()
          ongoing_pred_df <- tibble()
          for (i in 1:k) {
            event_pred_df <- event_pred_df %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$event_pred_df))

            ongoing_pred_df <- ongoing_pred_df %>%
              dplyr::bind_rows(dplyr::bind_cols(
                treatment = i, event_pred[[i]]$ongoing_pred_df))
          }

          event_pred_df <- event_pred_df %>%
            dplyr::bind_rows(event_pred_df1) %>%
            dplyr::arrange(treatment, t)

          ongoing_pred_df <- ongoing_pred_df %>%
            dplyr::bind_rows(ongoing_pred_df1) %>%
            dplyr::arrange(treatment, t)


          pred <- list(stage = "Real-time after enrollment completion",
                       to_predict = "Event only",
                       observed = observed,
                       event_fit = event_fit,
                       event_pred = event_pred[[1]])

          pred$event_pred$newEvents <- newEvents
          pred$event_pred$enroll_pred_df <- enroll_pred_df
          pred$event_pred$event_pred_df <- event_pred_df
          pred$event_pred$ongoing_pred_df <- ongoing_pred_df
          pred$event_pred$event_pred_day <- pred_day
          pred$event_pred$event_pred_date <- pred_date
          pred$event_pred$event_pred_summary <- s1

          pred

        }
      }

    }
  })




  dosing <- reactive({
    if (predict_dosing()) {
      req(pred()$event_pred)
      req(pred()$stage == input$stage &&
            pred()$to_predict == to_predict())

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
      }


      plower = (1 - pilevel())/2
      pupper = 1 - plower

      w = dosing_schedule()[,1]
      N = dosing_schedule()[,2]


      if (input$stage == 'Design stage') {

        newEvents <- pred()$event_pred$newEvents

        t = pred()$event_pred$event_pred_df$t

        df1 = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents)


        dfb <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df1 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime <= t) %>%
              dplyr::group_by(drug, t, draw) %>%
              dplyr::summarise(cum_dose = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d)),
                .groups = "drop_last")
          }))

        # predicted number of doses dispensed
        dosing_pred_df <- dfb %>%
          dplyr::group_by(drug, t) %>%
          dplyr::summarise(n = quantile(cum_dose, probs = 0.5),
                           lower = quantile(cum_dose, probs = plower),
                           upper = quantile(cum_dose, probs = pupper),
                           mean = mean(cum_dose),
                           var = var(cum_dose),
                           .groups = "drop_last")

        pred <- pred()
        pred$event_pred$dosing_pred_df <- dosing_pred_df
        pred
      } else if (input$stage == 'Real-time before enrollment completion') {
        df0 <- dplyr::tibble(drug = 1:l(),
                             t = 1, n = 0, lower = NA, upper = NA,
                             mean = 0, var = 0)

        df <- df() %>%
          dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                        totalTime = arrivalTime + time)

        trialsdt = pred()$observed$trialsdt
        t0 = pred()$observed$t0

        # observed dosing before data cut
        t = unique(c(seq(0, t0, 30), t0))

        df1 <- dplyr::tibble(t = t) %>%
          dplyr::cross_join(df)

        dosing_pred_obs_df <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df1 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime <= t) %>%
              dplyr::group_by(drug, t) %>%
              dplyr::summarise(n = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d)),
                .groups = "drop_last") %>%
              dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0)
          })) %>%
          dplyr::bind_rows(df0) %>%
          dplyr::arrange(drug, t)

        dosing_pred_t0 <- dosing_pred_obs_df %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::rename(cum_dose_t0 = n) %>%
          dplyr::select(drug, cum_dose_t0)


        # predicted dosing after data cut
        newEvents <- pred()$event_pred$newEvents
        t = pred()$event_pred$event_pred_df$t
        t = unique(t[t >= t0])

        df2 = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents)

        # from ongoing subjects
        dfa <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df2 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime <= t0 & totalTime > t0) %>%
              dplyr::group_by(drug, t, draw) %>%
              dplyr::summarise(cum_dose_a = sum(
                (f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d) -
                   f_cum_dose(t0 - arrivalTime, w, N, d))),
                .groups = "drop_last")
          }))

        # from new subjects
        dfb <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df2 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime > t0 & arrivalTime <= t) %>%
              dplyr::group_by(drug, t, draw) %>%
              dplyr::summarise(cum_dose_b = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d)),
                .groups = "drop_last")
          }))

        # dose dispense for time points after data cut
        dosing_pred_new_df <- dfa %>%
          dplyr::full_join(dfb, by = c("drug", "t", "draw")) %>%
          dplyr::left_join(dosing_pred_t0, by = 'drug') %>%
          dplyr::mutate(cum_dose = cum_dose_t0 +
                          ifelse(is.na(cum_dose_a), 0, cum_dose_a) +
                          ifelse(is.na(cum_dose_b), 0, cum_dose_b)) %>%
          dplyr::group_by(drug, t) %>%
          dplyr::summarise(n = quantile(cum_dose, probs = 0.5),
                           lower = quantile(cum_dose, probs = plower),
                           upper = quantile(cum_dose, probs = pupper),
                           mean = mean(cum_dose),
                           var = var(cum_dose),
                           .groups = "drop_last")

        # predicted dosing
        dosing_pred_df <- dosing_pred_obs_df %>%
          dplyr::bind_rows(dosing_pred_new_df) %>%
          dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
          dplyr::arrange(drug, t)

        pred <- pred()
        pred$event_pred$dosing_pred_df <- dosing_pred_df
        pred
      } else if (input$stage == "Real-time after enrollment completion") {
        df0 <- dplyr::tibble(drug = 1:l(),
                             t = 1, n = 0, lower = NA, upper = NA,
                             mean = 0, var = 0)

        df <- df() %>%
          dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                        totalTime = arrivalTime + time)

        trialsdt = pred()$observed$trialsdt
        t0 = pred()$observed$t0

        # observed dosing before data cut
        t = unique(c(seq(0, t0, 30), t0))

        df1 <- dplyr::tibble(t = t) %>%
          dplyr::cross_join(df)

        dosing_pred_obs_df <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df1 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime <= t) %>%
              dplyr::group_by(drug, t) %>%
              dplyr::summarise(n = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d)),
                .groups = "drop_last") %>%
              dplyr::mutate(lower = NA, upper = NA)
          })) %>%
          dplyr::bind_rows(df0) %>%
          dplyr::arrange(drug, t)

        dosing_pred_t0 <- dosing_pred_obs_df %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::rename(cum_dose_t0 = n) %>%
          dplyr::select(drug, cum_dose_t0)


        # predicted dosing after data cut
        newEvents <- pred()$event_pred$newEvents

        t0 = pred()$observed$t0
        t = pred()$event_pred$event_pred_df$t
        t = unique(t[t >= t0])

        df2 = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents)

        dfa <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            d = dosing_schedule()[,2+j]
            df2 %>%
              dplyr::right_join(treatment_by_drug() %>%
                                  dplyr::filter(drug == j),
                                by = "treatment") %>%
              dplyr::filter(arrivalTime <= t0 & totalTime > t0) %>%
              dplyr::group_by(drug, t, draw) %>%
              dplyr::summarise(cum_dose_a = sum(
                (f_cum_dose(pmin(totalTime, t) - arrivalTime, w, N, d) -
                   f_cum_dose(t0 - arrivalTime, w, N, d))),
                .groups = "drop_last")
          }))

        dosing_pred_new_df <- dfa %>%
          dplyr::left_join(dosing_pred_t0, by = 'drug') %>%
          dplyr::mutate(cum_dose = cum_dose_t0 +
                          ifelse(is.na(cum_dose_a), 0, cum_dose_a)) %>%
          dplyr::group_by(drug, t) %>%
          dplyr::summarise(n = quantile(cum_dose, probs = 0.5),
                           lower = quantile(cum_dose, probs = plower),
                           upper = quantile(cum_dose, probs = pupper),
                           mean = mean(cum_dose),
                           var = var(cum_dose),
                           .groups = "drop_last")

        # predicted dosing
        dosing_pred_df <- dosing_pred_obs_df %>%
          dplyr::bind_rows(dosing_pred_new_df) %>%
          dplyr::mutate(date = as.Date(t - 1, origin = trialsdt)) %>%
          dplyr::arrange(drug, t)

        pred <- pred()
        pred$event_pred$dosing_pred_df <- dosing_pred_df
        pred
      }
    }
  })



  output$dates <- renderText({
    if (!is.null(observed())) {
      str1 <- paste("Trial start date:", observed()$trialsdt)
      str2 <- paste("Data cutoff date:", observed()$cutoffdt)
      str3 <- paste("Days since trial start:", observed()$t0)
      paste(str1, str2, str3, sep='<br/>')
    }
  })


  output$statistics <- renderPrint({
    if (!is.null(df())) {

      if (input$by_treatment) {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          sum_by_trt <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(n0 = n(),
                             d0 = sum(event),
                             r0 = sum(!event))

          sum_overall <- sum_by_trt %>%
            dplyr::summarise(n0 = sum(n0),
                             d0 = sum(d0),
                             r0 = sum(r0))

          sumall <- sum_by_trt %>%
            dplyr::bind_rows(dplyr::bind_cols(treatment = 9999, sum_overall))

          table <- t(sumall %>%
                       dplyr::select(n0, d0, r0))
          colnames(table) <- paste("Treatment", sumall$treatment)
          colnames(table)[ncol(table)] <- "Overall"
          rownames(table) <- c("Current number of subjects",
                               "Current number of events",
                               "Number of ongoing subjects")
        } else {
          sum_by_trt <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(n0 = n())

          sum_overall <- sum_by_trt %>%
            dplyr::summarise(n0 = sum(n0))

          sumall <- sum_by_trt %>%
            dplyr::bind_rows(dplyr::bind_cols(treatment = 9999, sum_overall))

          table <- t(sumall %>% dplyr::select(n0))
          colnames(table) <- paste("Treatment", sumall$treatment)
          colnames(table)[ncol(table)] <- "Overall"
          rownames(table) <- c("Current number of subjects")
        }
      } else {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          table <- t(df() %>%
                       dplyr::summarise(n0 = n(),
                                        d0 = sum(event),
                                        r0 = sum(!event)))
          colnames(table) <- "Overall"
          rownames(table) <- c("Current number of subjects",
                               "Current number of events",
                               "Number of ongoing subjects")
        } else {
          table <- t(df() %>%
                       dplyr::summarise(n0 = n()))
          colnames(table) <- "Overall"
          rownames(table) <- c("Current number of subjects")
        }
      }

      print(table, quote=FALSE)
    }
  })


  output$cum_accrual_plot <- renderPlotly({
    cum_accrual_plot <- observed()$cum_accrual_plot
    if (!is.null(cum_accrual_plot)) {
      if (!input$by_treatment) {
        cum_accrual_plot
      } else {
        trialsdt = observed()$trialsdt
        cutoffdt = observed()$cutoffdt

        # enrollment data
        adsl <- df() %>%
          dplyr::group_by(treatment) %>%
          dplyr::arrange(randdt) %>%
          dplyr::mutate(adt = as.Date(time - 1, origin = randdt),
                        n = dplyr::row_number(),
                        parameter = "Enrollment",
                        date = randdt)

        # remove duplicate
        adslu <- adsl %>%
          dplyr::group_by(treatment, randdt) %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::select(treatment, n, parameter, date)

        # dummy subject to initialize time axis at trial start
        adsl0 <- dplyr::tibble(treatment = 1:k(),
                               n = 0,
                               parameter = "Enrollment",
                               date = trialsdt)

        # extend enrollment information to cutoff date
        adsl1 <- adsl %>%
          dplyr::group_by(treatment) %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::mutate(date = cutoffdt) %>%
          dplyr::select(treatment, n, parameter, date)

        if (grepl("event", to_predict(), ignore.case = TRUE)) {
          # time to event data
          adtte <- df() %>%
            dplyr::group_by(treatment) %>%
            dplyr::mutate(adt = as.Date(time - 1,
                                        origin = randdt)) %>%
            dplyr::arrange(adt) %>%
            dplyr::mutate(n = cumsum(event),
                          parameter = "Event",
                          date = adt)

          # remove duplicate
          adtteu <- adtte %>%
            dplyr::group_by(treatment, adt) %>%
            dplyr::slice(dplyr::n()) %>%
            dplyr::ungroup() %>%
            dplyr::select(treatment, n, parameter, date)

          # dummy subject to initialize time axis at trial start
          adtte0 <- dplyr::tibble(treatment = 1:k(),
                                  n = 0,
                                  parameter = "Event",
                                  date = trialsdt)

          # combine enrollment and time to event data
          ad <- adsl0 %>%
            dplyr::bind_rows(adslu) %>%
            dplyr::bind_rows(adsl1) %>%
            dplyr::bind_rows(adtte0) %>%
            dplyr::bind_rows(adtteu)
        } else {
          ad <- adsl0 %>%
            dplyr::bind_rows(adslu) %>%
            dplyr::bind_rows(adsl1)
        }


        ad <- ad %>%
          dplyr::mutate(treatmentc = paste0("treatment=", treatment))

        # plot cumulative enrollment and event data
        if (length(unique(ad$parameter)) > 1) {
          cum_accrual_plot <- plotly::plot_ly(
            ad, x=~date, y=~n, color=~parameter, colors=c("blue", "red"),
            linetype=~treatmentc) %>%
            plotly::add_lines(line = list(shape = "hv")) %>%
            plotly::layout(
              xaxis = list(title = ""),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'))
        } else {
          cum_accrual_plot <- plotly::plot_ly(
            ad, x=~date, y=~n, linetype=~treatmentc) %>%
            plotly::add_lines(line = list(shape = "hv")) %>%
            plotly::layout(
              xaxis = list(title = ""),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1, yanchor = "middle",
                            orientation = 'h'),
              title = list(text = "Cumulative enrollment"))
        }

        cum_accrual_plot
      }
    }
  })


  output$daily_accrual_plot <- renderPlotly({
    daily_accrual_plot <- observed()$daily_accrual_plot
    if (!is.null(daily_accrual_plot)) daily_accrual_plot
  })


  output$event_km_plot <- renderPlotly({
    event_km_plot <- observed()$event_km_plot
    if (!is.null(event_km_plot)) {
      if (!input$by_treatment) {
        event_km_plot
      } else {
        adtte <- observed()$adtte
        event_km_fit <- survival::survfit(survival::Surv(time, event) ~
                                            treatment, data = adtte)
        treatment <- attr(event_km_fit$strata, "names")

        event_km_df <- dplyr:: tibble(
          treatment = treatment,
          time = 0, surv = 1) %>%
          dplyr::bind_rows(dplyr::tibble(
            treatment = rep(treatment, event_km_fit$strata),
            time = event_km_fit$time,
            surv = event_km_fit$surv))

        event_km_plot <- plotly::plot_ly(event_km_df, x=~time, y=~surv,
                                         linetype=~treatment) %>%
          plotly::add_lines(line = list(shape = "hv")) %>%
          plotly::layout(xaxis = list(title = "Days since randomization",
                                      zeroline = FALSE),
                         yaxis = list(title = "Survival probability",
                                      zeroline = FALSE),
                         legend = list(x = 0, y = 1,  yanchor = "middle",
                                       orientation = 'h'),
                         title = list(
                           text = "Kaplan-Meier plot for time to event"))

        event_km_plot
      }
    }
  })



  output$input_df <- renderDataTable(
    df(), options = list(pageLength = 10)
  )



  output$enroll_fit <- renderPlotly({
    if (!is.null(enroll_fit())) enroll_fit()$enroll_fit_plot
  })


  output$event_fit1 <- renderPlotly({
    if (!is.null(event_fit())) {
      if (!input$by_treatment) {
        event_fit()$event_fit_plot
      } else {
        k = length(event_fit())
        event_fit_plot <- list()
        for (i in 1:k) {
          event_fit_plot[[i]] <- event_fit()[[i]]$event_fit_plot %>%
            plotly::layout(annotations = list(
              x = 0.5, y = 1, text = paste0("<b>treatment=", i, "</b>"),
              xanchor = "center", yanchor = "middle",
              showarrow = FALSE, xref='paper', yref='paper'))
        }

        plotly::subplot(event_fit_plot, nrows = k,
                        titleX = TRUE, titleY = TRUE, margin = 0.1)
      }
    }
  })


  output$event_fit <- renderUI({
    if (input$by_treatment && k() > 1) {
      plotlyOutput("event_fit1", height=240*k())
    } else {
      plotlyOutput("event_fit1")
    }
  })




  # enrollment predication date
  output$enroll_pred_date <- renderText({
    if (to_predict() == 'Enrollment only' ||
        to_predict() == 'Enrollment and event') {

      req(pred()$enroll_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())


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
    } else {
      text1 <- NULL
    }

    if (!is.null(text1)) text1
  })


  # event predication date
  output$event_pred_date <- renderText({
    if (to_predict() == 'Enrollment and event' ||
        to_predict() == 'Event only') {

      req(pred()$event_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))


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
    } else {
      text2 <- NULL
    }

    if (!is.null(text2)) text2

  })


  output$pred_date <- renderUI({
    if (to_predict() == 'Enrollment only') {
      htmlOutput("enroll_pred_date")
    } else if (to_predict() == 'Event only') {
      htmlOutput("event_pred_date")
    } else {
      fluidRow(column(6, htmlOutput("enroll_pred_date")),
               column(6, htmlOutput("event_pred_date")))
    }
  })


  # enrollment and event prediction plot
  output$pred_plot1 <- renderPlotly({
    if (to_predict() == "Enrollment only") {
      req(pred()$enroll_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))
      }

      enroll_pred_plot <- pred()$enroll_pred$enroll_pred_plot
      if ((!input$by_treatment || k() == 1) &&
          !("list" %in% class(enroll_pred_plot))) {
        p1 <- enroll_pred_plot
      } else if ((input$by_treatment && k() > 1) &&
                 ("list" %in% class(enroll_pred_plot)) &&
                 length(enroll_pred_plot) == k() + 1) {
        p1 <- plotly::subplot(enroll_pred_plot, nrows = k()+1, margin = 0.05)
      } else {
        p1 <- NULL
      }

    } else { # predict event only or predict enrollment and event
      shiny::validate(
        need(showEnrollment() || showEvent() || showOngoing(),
             "Need at least one parameter to show on prediction plot"))

      req(pred()$event_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())


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
      }


      dfs <- dplyr::tibble()

      if (showEnrollment()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$enroll_pred_df)
      }
      if (showEvent()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$event_pred_df)
      }
      if (showOngoing()) {
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)
      }

      dfs$parameter <- factor(dfs$parameter, levels = c(
        "Enrollment", "Event", "Ongoing"))


      if (input$stage != 'Design stage') {
        dfa <- dfs %>% dplyr::filter(is.na(lower))
        dfb <- dfs %>% dplyr::filter(!is.na(lower))

        if ((!input$by_treatment || k() == 1) &&
            !("treatment" %in% names(dfs))) {
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
                                 y = 0,
                                 text = 'cutoff', xanchor = "left",
                                 yanchor = "bottom",
                                 font = list(size = 12),
                                 showarrow = FALSE),
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'))

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
        } else if ((input$by_treatment && k() > 1) &&
                   ("treatment" %in% names(dfs)) &&
                   (length(unique(dfs$treatment)) == k() + 1)) {

          event_pred_plot <- list()

          # overall
          dfs0 <- dplyr::filter(dfs, treatment == 9999)
          dfb0 <- dplyr::filter(dfb, treatment == 9999)
          dfa0 <- dplyr::filter(dfa, treatment == 9999)

          event_pred_plot[[1]] <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfb0, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", fillcolor = ~parameter,
              line = list(width=0)) %>%
            plotly::add_lines(
              data = dfb0, x = ~date, y = ~n, color = ~parameter) %>%
            plotly::add_lines(
              data = dfa0, x = ~date, y = ~n, color = ~parameter,
              line = list(shape="hv", width=2)) %>%
            plotly::add_lines(
              x = rep(observed()$cutoffdt, 2), y = range(dfs0$n),
              name = "cutoff",
              line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              annotations = list(x = observed()$cutoffdt,
                                 y = 0,
                                 text = 'cutoff', xanchor = "left",
                                 yanchor = "bottom",
                                 font = list(size = 12),
                                 showarrow = FALSE),
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h')) %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1, text = "<b>Overall</b>",
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))

          if (showEvent()) {
            event_pred_plot[[1]] <- event_pred_plot[[1]] %>%
              plotly::add_lines(
                x = range(dfs0$date), y = rep(target_d(), 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
              plotly::layout(
                annotations = list(x = 0.95, xref = "paper", y = target_d(),
                                   text = 'target events',
                                   xanchor = "right", yanchor = "bottom",
                                   font = list(size = 12),
                                   showarrow = FALSE))
          }

          # by treatment
          for (i in 1:k()) {
            dfsi <- dplyr::filter(dfs, treatment == i)
            dfbi <- dplyr::filter(dfb, treatment == i)
            dfai <- dplyr::filter(dfa, treatment == i)

            event_pred_plot[[i+1]] <- plotly::plot_ly() %>%
              plotly::add_ribbons(
                data = dfbi, x = ~date, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", fillcolor = ~parameter,
                line = list(width=0)) %>%
              plotly::add_lines(
                data = dfbi, x = ~date, y = ~n, color = ~parameter) %>%
              plotly::add_lines(
                data = dfai, x = ~date, y = ~n, color = ~parameter,
                line = list(shape="hv", width=2)) %>%
              plotly::add_lines(
                x = rep(observed()$cutoffdt, 2), y = range(dfsi$n),
                name = "cutoff",
                line = list(dash="dash"),
                showlegend = FALSE) %>%
              plotly::layout(
                xaxis = list(title = "", zeroline = FALSE),
                yaxis = list(zeroline = FALSE),
                legend = list(x = 0, y = 1.05, yanchor = "bottom",
                              orientation = 'h'),
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>treatment=", i, "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))
          }

          p1 <- plotly::subplot(event_pred_plot, nrows = k()+1, margin = 0.05)

        } else {
          p1 = NULL
        }

      } else {  # Design stage

        if ((!input$by_treatment || k() == 1) &&
            !("treatment" %in% names(dfs))) {
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
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'))

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
        } else if ((input$by_treatment && k() > 1) &&
                   ("treatment" %in% names(dfs)) &&
                   (length(unique(dfs$treatment)) == k() + 1)) {

          event_pred_plot <- list()

          # overall
          event_pred_plot[[1]] <- dfs %>%
            dplyr::filter(treatment == 9999) %>%
            plotly::plot_ly() %>%
            plotly::add_ribbons(
              x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", fillcolor = ~parameter,
              line = list(width=0)) %>%
            plotly::add_lines(
              x = ~t, y = ~n, color = ~parameter) %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'),
              annotations = list(
                x = 0.5, y = 1, text = "<b>Overall</b>",
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))

          if (showEvent()) {
            event_pred_plot[[1]] <- event_pred_plot[[1]] %>%
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

          # by treatment
          for (i in 1:k()) {
            event_pred_plot[[i+1]] <- dfs %>%
              dplyr::filter(treatment == i) %>%
              plotly::plot_ly() %>%
              plotly::add_ribbons(
                x = ~t, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", fillcolor = ~parameter,
                line = list(width=0)) %>%
              plotly::add_lines(
                x = ~t, y = ~n, color = ~parameter) %>%
              plotly::layout(
                xaxis = list(title = "Days since trial start",
                             zeroline = FALSE),
                yaxis = list(zeroline = FALSE),
                legend = list(x = 0, y = 1.05, yanchor = "bottom",
                              orientation = 'h'),
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>treatment=", i, "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))
          }

          p1 <- plotly::subplot(event_pred_plot, nrows = k()+1, margin = 0.05)
        } else {
          p1 <- NULL
        }

      }
    }

    p1

  })


  output$pred_plot <- renderUI({
    if (input$by_treatment && k() > 1) {
      plotlyOutput("pred_plot1", height=250*(k()+1))
    } else {
      plotlyOutput("pred_plot1")
    }
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
          dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)
      }
      writexl::write_xlsx(sumdata, file)
    }
  )

  output$downloadSimdata <- downloadHandler(
    filename = function() {
      paste0("simdata-", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      if (to_predict() == "Enrollment only") {
        simdata <- pred()$enroll_pred$newSubjects
      } else {
        simdata <- pred()$event_pred$newEvents
      }
      writexl::write_xlsx(simdata, file)
    }
  )


  output$dosing_plot1 <- renderPlotly({
    if (predict_dosing()) {
      req(dosing()$event_pred)
      req(dosing()$stage == input$stage &&
            dosing()$to_predict == to_predict())

      df <- dosing()$event_pred$dosing_pred_df

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

        dosing_pred_plot <- list()

        for (j in 1:l()) {
          dfs <- dplyr::filter(df, drug == j)
          dfa <- dplyr::filter(df, drug == j & is.na(lower))
          dfb <- dplyr::filter(df, drug == j & !is.na(lower))

          dosing_pred_plot[[j]] <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval") %>%
            plotly::add_lines(
              data = dfb, x = ~date, y = ~n, name = "median prediction") %>%
            plotly::add_lines(
              data = dfa, x = ~date, y = ~n, name = "observed") %>%
            plotly::add_lines(
              x = rep(observed()$cutoffdt, 2), y = range(dfs$n),
              name = "cutoff",
              line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(title = "Doses to dispense",
                           zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'),
              annotations = list(
                x = 0.5, y = 1, text = paste0("<b>drug=", j, "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))

          if (j==1) {
            dosing_pred_plot[[j]] <- dosing_pred_plot[[j]] %>%
              plotly::layout(
                annotations = list(x = observed()$cutoffdt,
                                   y = 0,
                                   text = 'cutoff', xanchor = "left",
                                   yanchor = "bottom",
                                   font = list(size = 12),
                                   showarrow = FALSE))
          }
        }

        p1 <- plotly::subplot(dosing_pred_plot, nrows = l(),
                              titleX = TRUE, titleY = TRUE, margin = 0.1)
      } else {
        dosing_pred_plot <- list()

        for (j in 1:l()) {
          dosing_pred_plot[[j]] <- df %>%
            dplyr::filter(drug == j) %>%
            plotly::plot_ly() %>%
            plotly::add_ribbons(
              name = "prediction interval",
              x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0)) %>%
            plotly::add_lines(
              x = ~t, y = ~n, name = "median prediction") %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(title = "Doses to dispense",
                           zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'),
              annotations = list(
                x = 0.5, y = 1, text = paste0("<b>drug=", j, "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))
        }

        p1 <- plotly::subplot(dosing_pred_plot, nrows = l(),
                              titleX = TRUE, titleY = TRUE, margin = 0.1)
      }

    } else {
      p1 <- NULL
    }

    p1
  })



  output$dosing_plot <- renderUI({
    if (l() > 1) {
      plotlyOutput("dosing_plot1", height=250*l())
    } else {
      plotlyOutput("dosing_plot1")
    }
  })


  output$downloadDosingdata <- downloadHandler(
    filename = function() {
      paste0("dosingdata-", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      dosingdata <- dosing()$event_pred$dosing_pred_df
      writexl::write_xlsx(dosingdata, file)
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
    b = matrix(a[nrow(a),], nrow=1)
    b[1,1] = b[1,1] + 1
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


  lapply(1:6, function(i) {
    pwexp <- paste0("piecewise_exponential_survival_", i)
    observeEvent(input[[paste0("add_piecewise_exponential_survival_", i)]], {
      a = matrix(as.numeric(input[[pwexp]]), ncol=ncol(input[[pwexp]]))
      b = matrix(a[nrow(a),], nrow=1)
      b[1,1] = b[1,1] + 1
      c = rbind(a, b)
      rownames(c) = paste("Interval", seq(1:nrow(c)))
      colnames(c) = colnames(input[[pwexp]])
      updateMatrixInput(session, pwexp, c)
    })
  })


  lapply(1:6, function(i) {
    pwexp <- paste0("piecewise_exponential_survival_", i)
    observeEvent(input[[paste0("del_piecewise_exponential_survival_", i)]], {
      if (nrow(input[[pwexp]]) >= 2) {
        a = matrix(as.numeric(input[[pwexp]]), ncol=ncol(input[[pwexp]]))
        b = matrix(a[-nrow(a),], ncol=ncol(a))
        rownames(b) = paste("Interval", seq(1:nrow(b)))
        colnames(b) = colnames(input[[pwexp]])
        updateMatrixInput(session, pwexp, b)
      }
    })
  })



  lapply(1:12, function(j) {
    dosing <- paste0("dosing_schedule_", j)
    observeEvent(input[[paste0("add_dosing_schedule_", j)]], {
      a = matrix(as.numeric(input[[dosing]]), ncol=ncol(input[[dosing]]))
      b = matrix(a[nrow(a),], nrow=1)
      c = rbind(a, b)
      rownames(c) = paste("Interval", seq(1:nrow(c)))
      colnames(c) = colnames(input[[dosing]])
      updateMatrixInput(session, dosing, c)
    })
  })


  lapply(1:12, function(j) {
    dosing <- paste0("dosing_schedule_", j)
    observeEvent(input[[paste0("del_dosing_schedule_", j)]], {
      if (nrow(input[[dosing]]) >= 2) {
        a = matrix(as.numeric(input[[dosing]]), ncol=ncol(input[[dosing]]))
        b = matrix(a[-nrow(a),], ncol=ncol(a))
        rownames(b) = paste("Interval", seq(1:nrow(b)))
        colnames(b) = colnames(input[[dosing]])
        updateMatrixInput(session, dosing, b)
      }
    })
  })



  # save inputs
  output$saveInputs <- downloadHandler(
    filename = function() {
      paste0("inputs-", Sys.Date(), ".rds")
    },

    content = function(file) {
      x <- list(
        treatment_allocation = matrix(
          treatment_allocation(), ncol=1,
          dimnames = list(paste("Treatment", 1:k()), "Size")),
        exponential_survival = matrix(
          exponential_survival(), nrow = 1,
          dimnames = list(NULL, paste("Treatment", 1:k()))),
        weibull_survival = matrix(
          weibull_survival(), nrow = 2,
          dimnames = list(c("Shape", "Scale"), paste("Treatment", 1:k()))),
        lnorm_survival = matrix(
          lnorm_survival(), nrow = 2,
          dimnames = list(c("Mean on log scale", "SD on log scale"),
                          paste("Treatment", 1:k()))),
        piecewise_exponential_survival = matrix(
          piecewise_exponential_survival(), ncol = k()+1,
          dimnames = list(paste("Interval",
                                1:nrow(piecewise_exponential_survival())),
                          c("Starting time", paste("Treatment", 1:k())))),
        treatment_by_drug = matrix(
          as.numeric(input[[paste0("treatment_by_drug_", k(), "_", l())]]),
          ncol = l(),
          dimnames = list(paste("Treatment", 1:k()), paste("Drug", 1:l()))),
        dosing_schedule = matrix(
          dosing_schedule(), ncol = l() + 2,
          dimnames = list(paste("Interval", 1:nrow(dosing_schedule())),
                          c("Frequency in weeks",
                            "Number of cycles",
                            paste("Drug", 1:l(), "dose")))),
        enroll_prior = input$enroll_prior,
        poisson_rate = poisson_rate(),
        mu = mu(),
        delta = delta(),
        piecewise_poisson_rate = matrix(
          piecewise_poisson_rate(), ncol = 2,
          dimnames = list(paste("Interval", 1:nrow(piecewise_poisson_rate())),
                          c("Starting time", "Enrollment rate"))),
        enroll_model = input$enroll_model,
        nknots = nknots(),
        lags = lags(),
        accrualTime = matrix(
          accrualTime(), ncol = 1,
          dimnames = list(paste("Interval", 1:length(accrualTime())),
                          "Starting time")),
        event_prior = input$event_prior,
        event_model = input$event_model,
        piecewiseSurvivalTime = matrix(
          piecewiseSurvivalTime(), ncol = 1,
          dimnames = list(paste("Interval", 1:length(piecewiseSurvivalTime())),
                          "Starting time")),
        spline_k = spline_k(),
        spline_scale = input$spline_scale,
        l = l(),
        stage = input$stage,
        to_predict = input$to_predict,
        to_predict2 = input$to_predict2,
        target_n = target_n(),
        target_d = input$target_d,
        pilevel = pilevel(),
        nyears = nyears(),
        to_show = input$to_show,
        by_treatment = input$by_treatment,
        predict_dosing = predict_dosing(),
        k = k(),
        nreps = nreps(),
        seed = input$seed
      )

      save(x, file = file)
    }
  )


  # load inputs
  observeEvent(input$loadInputs, {
    file <- input$loadInputs
    ext <- tools::file_ext(file$datapath)

    req(file)

    valid <- (ext == "rds")
    if (!valid) showNotification("Please upload an rds file")
    req(valid)

    load(file=file$datapath)

    if (x$stage == 'Design stage' ||
        (x$by_treatment &&
         x$stage != 'Real-time after enrollment completion')) {
      updateMatrixInput(
        session, paste0("treatment_allocation_", x$k),
        value=matrix(x$treatment_allocation, ncol = 1,
                     dimnames = list(paste("Treatment", 1:x$k),
                                     "Size")))
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Exponential') {
      updateMatrixInput(
        session, paste0("exponential_survival_", x$k),
        value=matrix(x$exponential_survival, ncol = x$k,
                     dimnames = list(NULL, paste("Treatment", 1:x$k))))
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Weibull') {
      updateMatrixInput(
        session, paste0("weibull_survival_", x$k),
        value=matrix(x$weibull_survival, ncol = x$k,
                     dimnames = list(c("Shape", "Scale"),
                                     paste("Treatment", 1:x$k))))
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Log-normal') {
      updateMatrixInput(
        session, paste0("lnorm_survival_", x$k),
        value=matrix(x$lnorm_survival, ncol = x$k,
                     dimnames = list(c("Mean on log scale",
                                       "SD on log scale"),
                                     paste("Treatment", 1:x$k))))
    }

    if (x$stage == 'Design stage' &&
        x$event_prior == 'Piecewise exponential') {
      updateMatrixInput(
        session, paste0("piecewise_exponential_survival_", x$k),
        value=matrix(x$piecewise_exponential_survival, ncol = x$k + 1,
                     dimnames = list(
                       paste("Interval",
                             1:nrow(x$piecewise_exponential_survival)),
                       c("Starting time", paste("Treatment", 1:x$k)))))
    }


    if (x$predict_dosing) {
      updateMatrixInput(
        session, paste0("treatment_by_drug_", x$k, "_", x$l),
        value=matrix(x$treatment_by_drug, ncol = x$l,
                     dimnames = list(paste("Treatment", 1:x$k),
                                     paste("Drug", 1:x$l))))
      updateMatrixInput(
        session, paste0("dosing_schedule_", x$l),
        value=matrix(x$dosing_schedule, ncol = x$l + 2,
                     dimnames = list(paste("Interval",
                                           1:nrow(x$dosing_schedule)),
                                     c("Frequency in weeks",
                                       "Number of cycles",
                                       paste("Drug", 1:x$l, "dose")))))
      updateSelectInput(session, "l", selected=x$l)
    }

    if (x$stage == 'Design stage') {
      updateRadioButtons(session, "enroll_prior", selected=x$enroll_prior)

      if (x$enroll_prior == "Poisson") {
        updateNumericInput(session, "poisson_rate", value=x$poisson_rate)
      } else if (x$enroll_prior == "Time-decay") {
        updateNumericInput(session, "mu", value=x$mu)
        updateNumericInput(session, "delta", value=x$delta)
      } else if (x$enroll_prior == "Piecewise Poisson") {
        updateMatrixInput(
          session, "piecewise_poisson_rate",
          value=matrix(x$piecewise_poisson_rate, ncol = 2,
                       dimnames = list(
                         paste("Interval",
                               1:nrow(x$piecewise_poisson_rate)),
                         c("Starting time", "Enrollment rate"))))
      }

      if (x$to_predict == 'Enrollment and event') {
        updateRadioButtons(session, "event_prior", selected=x$event_prior)
      }

    } else {

      if (x$stage == 'Real-time before enrollment completion') {
        updateRadioButtons(session, "enroll_model", selected=x$enroll_model)

        if (x$enroll_model == "B-spline") {
          updateNumericInput(session, "nknots", value=x$nknots)
          updateNumericInput(session, "lags", value=x$lags)
        } else if (x$enroll_model == "Piecewise Poisson") {
          updateMatrixInput(
            session, "accrualTime",
            value=matrix(x$accrualTime, ncol = 1,
                         dimnames = list(
                           paste("Interval", 1:nrow(x$accrualTime)),
                           "Starting time")))
        }
      }


      if ((x$stage == 'Real-time before enrollment completion' &&
           x$to_predict == 'Enrollment and event') ||
          x$stage == 'Real-time after enrollment completion') {

        updateRadioButtons(session, "event_model", selected=x$event_model)

        if (x$event_model == "Piecewise exponential") {
          updateMatrixInput(
            session, "piecewiseSurvivalTime",
            value=matrix(x$piecewiseSurvivalTime, ncol = 1,
                         dimnames = list(
                           paste("Interval",
                                 1:nrow(x$piecewiseSurvivalTime)),
                           "Starting time")))
        } else if (x$event_model == "Spline") {
          updateNumericInput(session, "spline_k", value=x$spline_k)
          updateRadioButtons(session, "spline_scale", value=x$spline_scale)
        }
      }


    }


    updateRadioButtons(session, "stage", selected=x$stage)

    if (x$stage == 'Design stage' ||
        x$stage == 'Real-time before enrollment completion') {
      updateRadioButtons(session, "to_predict", selected=x$to_predict)
      updateNumericInput(session, "target_n", value=x$target_n)
    } else {
      updateRadioButtons(session, "to_predict2", selected=x$to_predict2)
    }

    if (x$to_predict == 'Enrollment and event' ||
        x$stage == 'Real-time after enrollment completion') {
      updateNumericInput(session, "target_d", value=x$target_d)
      updateCheckboxGroupInput(session, "to_show", selected=x$to_show)
    }

    updateNumericInput(session, "pilevel", value=x$pilevel)
    updateNumericInput(session, "nyears", value=x$nyears)

    updateCheckboxInput(session, "by_treatment", value=x$by_treatment)

    if (x$by_treatment &&
        (x$to_predict == 'Enrollment and event' ||
         x$stage == 'Real-time after enrollment completion')) {
      updateCheckboxInput(session, "predict_dosing", value=x$predict_dosing)
    }


    if (x$stage == 'Design stage' ||
        (x$by_treatment &&
         x$stage != 'Real-time after enrollment completion')) {
      updateSelectInput(session, "k", selected=x$k)
    }

    updateNumericInput(session, "nreps", value=x$nreps)
    updateNumericInput(session, "seed", value=x$seed)

  })
}

# Run the application
shinyApp(ui = ui, server = server)

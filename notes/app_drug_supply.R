library(shiny)
library(shinyMatrix)
library(shinyFeedback)
library(shinyjs, warn.conflicts = FALSE)
library(shinybusy)
library(readxl)
library(writexl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(prompter)
library(ggplot2)
library(plotly, warn.conflicts = FALSE)
library(eventPred)


rdirichlet <- function (n = 1, alpha) {
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}


# conditional panels for treatment allocation
f_treatment_allocation <- function(i) {
  conditionalPanel(
    condition = paste0("input.k == ", i),

    shinyMatrix::matrixInput(
      paste0("treatment_allocation_", i),
      label = tags$span(
        "Treatment allocation",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = "in a randomization block",
                     position = "right")),

      value = matrix(rep(1,i), ncol = 1,
                     dimnames = list(paste("Treatment", 1:i), "Size")),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE, editableNames=TRUE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Exponential' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("exponential_survival_", i),
      label = "Hazard rate for each treatment",
      value = matrix(rep(0.0030, i), nrow = 1,
                     dimnames = list(NULL, paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=FALSE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_weibull_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Weibull' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("weibull_survival_", i),
      label = "Weibull parameters",
      value = matrix(rep(c(1.42, 392), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Shape", "Scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_llogis_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Log-logistic' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("llogis_survival_", i),
      label = "Log-logistic parameters",
      value = matrix(rep(c(5.4, 1), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Location on log scale",
                                       "Scale on log scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}

f_lnorm_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Log-normal' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("lnorm_survival_", i),
      label = "Log-normal parameters",
      value = matrix(rep(c(5.4, 1), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Mean on log scale",
                                       "SD on log scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_piecewise_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste(
      "input.event_prior == 'Piecewise exponential' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("piecewise_exponential_survival_", i),
      label = "Hazard rate by time interval for each treatment",
      value = matrix(c(0, rep(0.0030, i)), nrow = 1,
                     dimnames = list(
                       "Interval 1",
                       c("Starting time", paste("Treatment", 1:i)))),
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


f_drug_description <- function(j) {
  conditionalPanel(
    condition = paste0("input.l == ", j),

    shinyMatrix::matrixInput(
      paste0("drug_description_", j),
      label = tags$span(
        "Drug names and dose units",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = paste(
            "enter drugs determined by treatment arms",
            "followed by optional chemotherapy drugs if any"),
                     position = "right")),
      value = matrix(c(paste("Drug", seq_len(j)), rep("mg", j)),
                     nrow = j, ncol = 2,
                     dimnames = list(NULL, c("Drug Name", "Dose Unit"))),
      inputClass = "character",
      rows = list(names=FALSE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_treatment_by_drug <- function(i, j) {
  conditionalPanel(
    condition = paste("input.k ==", i, "&&", "input.l - input.l2 ==", j),

    shinyMatrix::matrixInput(
      paste0("treatment_by_drug_", i, "_", j),
      label = "Drugs contained in each treatment",
      value = matrix(1, nrow = i, ncol = j,
                     dimnames = list(paste("Treatment", 1:i),
                                     paste("Drug", 1:j))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


# conditional panels for chemotherapy occupancy
f_chemotherapy_occupancy <- function(i) {
  conditionalPanel(
    condition = paste0("input.stage == 'Design stage' && input.m == ", i),

    shinyMatrix::matrixInput(
      paste0("chemotherapy_occupancy_", i),
      label = tags$span(
        "Chemotherapy occupancy",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = paste(
            "no need to enter the probability for the last chemotherapy",
            "as it will be calculated automatically"),
            position = "left")),
      value = matrix(rep(1/i,i), ncol = 1,
                     dimnames = list(paste("Chemotherapy", 1:i),
                                     "Probability")),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE, editableNames=TRUE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_chemotherapy_by_drug <- function(i, j) {
  conditionalPanel(
    condition = paste("input.m ==", i, "&&", "input.l2 ==", j),

    shinyMatrix::matrixInput(
      paste0("chemotherapy_by_drug_", i, "_", j),
      label = "Drugs contained in each chemotherapy",
      value = matrix(1, nrow = i, ncol = j,
                     dimnames = list(paste("Chemotherapy", 1:i),
                                     paste("Drug", 1:j))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_dosing_schedule <- function(j) {
  lapply(1:j, function(i) {
    conditionalPanel(
      condition = paste("input.l ==", j),

      htmlOutput(paste0("label_dosing_schedule_", j, "_", i)),

      shinyMatrix::matrixInput(
        paste0("dosing_schedule_", j, "_", i),
        label = NULL,
        value = matrix(c(21, 100, 10000),
                       nrow = 1, ncol = 3,
                       dimnames = list(
                         "Interval 1",
                         c("Days per Cycle", "Dose per Cycle",
                           "Number of Cycles"))),
        inputClass = "numeric",
        rows = list(names=TRUE, extend=FALSE, editableNames=TRUE),
        cols = list(names=TRUE, extend=FALSE)
      ),

      actionButton(paste0("add_dosing_schedule_", j, "_", i),
                   label=NULL, icon=icon("plus")),
      actionButton(paste0("del_dosing_schedule_", j, "_", i),
                   label=NULL, icon=icon("minus"))
    )
  })
}


# cumulative dose for duration x and a drug with dosing schedule (w, d, N)
# here x can be a vector
f_cum_dose <- function(x, w, d, N) {
  m = length(w)
  u = c(0, cumsum(w*N))
  i = pmin(findInterval(x, u), m)
  z = pmin(x, u[m+1]) - u[i]
  n = pmin(floor(z/w[i]) + 1, N[i])
  v = c(0, cumsum(d*N))
  v[i] + n*d[i]
}


# detailed dosing for duration x and a drug with dosing schedule (w, d, N)
# here we assume that x is a scalar
g_cum_dose <- function(x, w, d, N) {
  m = length(w)
  u = c(0, cumsum(w*N))
  i = min(findInterval(x, u), m)
  z = min(x, u[m+1]) - u[i]
  n = min(floor(z/w[i]) + 1, N[i])

  N1 = N[1:i]
  N1[i] = n # number of cycles in the last interval

  w1 <- w0 <- rep(w[1:i], times = N1)
  w1[length(w1)] = z %% w[i] # days in last cycle
  w2 = c(0, cumsum(w1)) + 1 # relative day at start of cycle

  dplyr::tibble(
    day_rel_randdt = w2[-length(w2)],
    days_per_cycle = w0,
    dose_to_dispense = rep(d[1:i], times = N1),
    cum_dose = cumsum(dose_to_dispense))
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
        label = "Which enrollment model to use?",
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
                 label = "Daily enrollment rate",
                 value = 1,
                 min = 0, max = 100, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_prior == 'Time-decay'",

               fluidRow(
                 column(6, numericInput(
                   "mu",
                   label = "Base rate, mu",
                   value = 1.5,
                   min = 0, max = 100, step = 1)
                 ),

                 column(6, numericInput(
                   "delta",
                   label = "Decay rate, delta",
                   value = 2,
                   min = 0, max = 100, step = 1)
                 )
               )
             ),

             conditionalPanel(
               condition = "input.enroll_prior == 'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "piecewise_poisson_rate",
                 label = "Daily enrollment rate by time interval",
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
        label = "Which enrollment model to use?",
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
                 label = "How many inner knots to use?",
                 value = 0,
                 min = 0, max = 10, step = 1),

               numericInput(
                 "lags",
                 label = paste("How many days before the last enrollment",
                               "date to average",
                               "the enrollment rate over for prediction?"),
                 value = 30,
                 min = 0, max = 365, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_model == 'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "accrualTime",
                 label = "What is the starting time of each time interval?",
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
        label = "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-logistic",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Piecewise exponential",
        inline = FALSE)
      ),

      column(8,
             lapply(1:6, f_exponential_survival),
             lapply(1:6, f_weibull_survival),
             lapply(1:6, f_llogis_survival),
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
        label = "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-logistic",
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
                 label = "What is the starting time of each time interval?",
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
                 label = "How many inner knots to use?",
                 value = 0,
                 min = 0, max = 10, step = 1),

               radioButtons(
                 "spline_scale",
                 label = "Which scale to model as a spline function?",
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
    "l", label = "Total number of drugs",
    choices = seq_len(12), selected = 2),

  fluidRow(
    column(6, selectInput(
      "l2",
      label = "Number of optional chemotherapy drugs",
      choices = c(0, 2:6), selected = 0),

      conditionalPanel(
        condition = "input.l2 > 0",
        selectInput(
          "m", label = "Number of optional chemotherapies",
          choices = 2:6, selected = 2))),

    column(6, conditionalPanel(
      condition = "input.l2 > 0",
      lapply(2:6, f_chemotherapy_occupancy)
    ))),

  lapply(1:12, f_drug_description),

  lapply(c(outer(1:6, 1:12, function(i,j) 1000*i+j)),
         function(k) f_treatment_by_drug(k %/% 1000, k %% 1000)),

  lapply(c(outer(1:6, 1:12, function(i,j) 1000*i+j)),
         function(k) f_chemotherapy_by_drug(k %/% 1000, k %% 1000)),

  lapply(1:12, f_dosing_schedule)
)


eventPredictPanel <- tabPanel(
  title = "Event Prediction",
  value = "event_prediction_panel",

  uiOutput("pred_date"),
  uiOutput("pred_plot"),

  downloadButton("downloadEventSummaryData", "Download summary data"),
  downloadButton("downloadEventSubjectData", "Download subject data")
)


dosingPredictPanel <- tabPanel(
  title = "Dosing Prediction",
  value = "dosing_prediction_panel",

  uiOutput("dosing_plot"),
  downloadButton("downloadDosingSummaryData", "Download summary data"),

  conditionalPanel(
    condition = "input.subject_dosing",
    downloadButton("downloadDosingSubjectData", "Download subject data")
  )
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
      paste0("width: ", validateCssUnit(width), ";",
             "padding-right: 5px; padding-bottom: 0px;
             display:inline-block;"),

    span(class = "btn btn-default btn-file",type="button",
         buttonLabel, inputTag,
         style=if (!is.null(width))
           paste0("width: ", validateCssUnit(width), ";",
                  "border-radius: 4px; padding-bottom:5px;"))
  )
}


# user interface ----------------
ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  prompter::use_prompt(),
  shinybusy::add_busy_spinner(),

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
        label = "Stage of the study",
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
          label = "What to predict?",
          choices = c("Enrollment only",
                      "Enrollment and event"),
          selected = "Enrollment and event",
          inline = FALSE)
      ),


      conditionalPanel(
        condition =
          "input.stage == 'Real-time after enrollment completion'",

        radioButtons(
          "to_predict2",
          label = "What to predict?",
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
                   label = "Target enrollment",
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
                   label = "Target events",
                   value = 200,
                   min = 1, max = 10000, step = 1)
               )
        )
      ),


      conditionalPanel(
        condition = "input.stage != 'Design stage'",

        fileInput(
          "file1",
          label = "Upload subject level data",
          accept = ".xlsx"
        )
      ),


      fluidRow(
        column(7, radioButtons(
          "pilevel",
          label = "Prediction interval",

          choices = c("95%" = "0.95", "90%" = "0.90", "80%" = "0.80"),
          selected = "0.95",
          inline = TRUE)
        ),

        column(5, numericInput(
          "nyears",
          label = "Years after cutoff",
          value = 1,
          min = 1, max = 10, step = 1)
        )
      ),


      conditionalPanel(
        condition = "input.to_predict == 'Enrollment and event' ||
        input.stage == 'Real-time after enrollment completion'",

        checkboxGroupInput(
          "to_show",
          label = "What to show on prediction plot?",
          choices = c("Enrollment", "Event", "Ongoing"),
          selected = c("Enrollment", "Event"),
          inline = TRUE
        )
      ),


      fluidRow(
        column(7, checkboxInput(
          "by_treatment", label = "By treatment?", value = TRUE),

          conditionalPanel(
            condition = "input.by_treatment &&
            (input.to_predict == 'Enrollment and event' ||
            input.stage == 'Real-time after enrollment completion')",

            checkboxInput(
              "predict_dosing", label = "Predict dosing?", value = TRUE),

            conditionalPanel(
              condition = "input.predict_dosing",
              checkboxInput(
                "subject_dosing", label = "Subject-level dosing?",
                value = FALSE)
            )
          )
        ),


        column(5, conditionalPanel(
          condition = "input.stage == 'Design stage' || input.by_treatment",

          selectInput(
            "k", label = "Treatments", choices = seq_len(6), selected = 2))
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
  # session$onSessionEnded(function() {
  #   stopApp()
  # })

  # whether to show or hide the observed data panel
  observeEvent(input$stage, {
    if (input$stage != 'Design stage') {
      showTab(inputId = "results", target = "observed_data_panel")
    } else {
      hideTab(inputId = "results", target = "observed_data_panel")
    }
  })


  # whether to allow the user to specify the number of treatments
  observeEvent(input$stage, {
    shinyjs::toggleState("k", input$stage == "Design stage")
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
        (to_predict() == 'Enrollment and event' ||
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


  # whether to generate subject-level dosing data
  subject_dosing <- reactive({
    if (input$predict_dosing) {
      input$subject_dosing
    } else {
      FALSE
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

    if (to_predict() == "Enrollment and event") {
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


  k <- reactive({
    if (!input$by_treatment && input$stage != "Design stage") {
      k = 1
    } else if (input$stage != "Design stage" && !is.null(df())) {
      k = length(table(df()$treatment))
      updateSelectInput(session, "k", selected=k)
    } else {
      k = as.numeric(input$k)
    }
    k
  })


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


  treatment_description <- reactive({
    req(k())
    if (k() > 1) {
      if (!input$by_treatment && input$stage != "Design stage") {
        a = "Overall"
      } else if (input$stage != "Design stage" && !is.null(df())) {
        treatment_mapping <- df() %>%
          dplyr::select(treatment, treatment_description) %>%
          dplyr::arrange(treatment, treatment_description) %>%
          dplyr::group_by(treatment, treatment_description) %>%
          dplyr::slice(dplyr::n())
        a = treatment_mapping$treatment_description
      } else {
        a = rownames(input[[paste0("treatment_allocation_", k())]])
      }
    } else {
      a = "Overall"
    }
    a
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

    matrix(c(t, lambda), ncol = 2,
           dimnames = list(paste("Interval", 1:length(t)),
                           c("Starting time", "Enrollment rate")))
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


  llogis_survival <- reactive({
    req(k())
    param = input[[paste0("llogis_survival_", k())]]
    locationlog = as.numeric(param[1,])
    scalelog = as.numeric(param[2,])

    valid = all(scalelog > 0)
    if (!valid) {
      showNotification(
        "Scale on the log scale must be positive"
      )
    }

    req(valid)

    matrix(c(locationlog, scalelog), nrow = 2, byrow = TRUE)
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


  observeEvent(treatment_description(), {
    if (input$stage == "Design stage") {
      updateMatrixInput(
        session, paste0("exponential_survival_", k()),
        value=matrix(exponential_survival(), ncol=k(),
                     dimnames = list(NULL, treatment_description())))
      updateMatrixInput(
        session, paste0("weibull_survival_", k()),
        value=matrix(weibull_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Shape", "Scale"),
                                     treatment_description())))
      updateMatrixInput(
        session, paste0("llogis_survival_", k()),
        value=matrix(llogis_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Location on log scale",
                                       "Scale on log scale"),
                                     treatment_description())))
      updateMatrixInput(
        session, paste0("lnorm_survival_", k()),
        value=matrix(lnorm_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Mean on log scale",
                                       "SD on log scale"),
                                     treatment_description())))

      npieces = nrow(piecewise_exponential_survival())
      updateMatrixInput(
        session, paste0("piecewise_exponential_survival_", k()),
        value=matrix(piecewise_exponential_survival(),
                     nrow=npieces, ncol=k()+1,
                     dimnames = list(
                       paste("Interval", seq_len(npieces)),
                       c("Starting time", treatment_description()))))
    } else if (input$by_treatment && !is.null(df())) {
      updateMatrixInput(
        session, paste0("treatment_allocation_", k()),
        value=matrix(treatment_allocation(), ncol = 1,
                     dimnames = list(treatment_description(), "Size")))
    }
  })



  l <- reactive(as.numeric(input$l))


  drug_name <- reactive({
    req(l())
    param = input[[paste0("drug_description_", l())]]
    as.character(param[,1])
  })


  dose_unit <- reactive({
    req(l())
    param = input[[paste0("drug_description_", l())]]
    as.character(param[,2])
  })


  dosing_schedule_df <- reactive({
    req(l())

    param <- dplyr::bind_rows(
      lapply(1:l(), function(i) {
        x <- input[[paste0("dosing_schedule_", l(), "_", i)]]
        dplyr::tibble(
          w = as.numeric(x[, "Days per Cycle"]),
          d = as.numeric(x[, "Dose per Cycle"]),
          N = as.numeric(x[, "Number of Cycles"]),
          drug = i,
          interval = rownames(x))
      }))

    valid1 = all(param$w > 0 & param$w == round(param$w))
    if (!valid1) {
      showNotification("Days per Cycle must be positive integers")
    }

    valid2 = all(param$N > 0 & param$N == round(param$N))
    if (!valid2) {
      showNotification("Number of Cycles must be positive integers")
    }

    valid3 = all(param$d >= 0)
    if (!valid3) {
      showNotification("Doses must be nonnegative")
    }

    req(valid1 && valid2 && valid3)

    param
  })


  observeEvent(input$l, {
    if (l() > 2) {
      choices <- c(0, 2:min((l()-1), 6))
    } else {
      choices <- 0
    }

    if (as.numeric(input$l2) < l()) {
      updateSelectInput(session, "l2", choices = choices, selected = input$l2)
    } else {
      updateSelectInput(session, "l2", choices = choices, selected = 0)
    }
  })


  l2 <- reactive({
    l2 = as.numeric(input$l2)
    ifelse(l2 >= l(), 0, l2)
  })


  l1 <- reactive(l() - l2())


  m <- reactive({
    req(l2())
    if (l2() > 0) {
      if (input$stage != "Design stage" && !is.null(df())) {
        m = length(table(df()$chemotherapy))
        updateSelectInput(session, "m", selected=m)
      } else {
        m = as.numeric(input$m)
      }
    } else {
      m = 1
    }
    m
  })


  chemotherapy_occupancy <- reactive({
    req(m())
    if (m() > 1 && input$stage == "Design stage") {
      d = input[[paste0("chemotherapy_occupancy_", m())]]
      d <- as.numeric(d)

      valid = all(d > 0) && (sum(d[-m()]) < 1)
      if (!valid) {
        showNotification(paste("Chemotherapy occupancy probabilities must be",
                               "positive and sum to 1"))
      }

      req(valid)
      d
    } else if (m() > 1 && !is.null(df())) {
      d = as.numeric(table(df()$chemotherapy))
      d/sum(d)
    } else if (m() > 1) {
      rep(1/m(), m())
    } else {
      0
    }
  })


  observeEvent(chemotherapy_occupancy(), {
    req(m())
    if (m() > 1 && input$stage == "Design stage") {
      a <- paste0("chemotherapy_occupancy_", m())
      prevValues <- as.numeric(input[[a]][-m()])
      lastValue <- 1 - sum(prevValues)
      b = matrix(c(prevValues, lastValue), nrow = m(), ncol = 1,
                 dimnames = dimnames(input[[a]]))
      updateMatrixInput(session, a, value=b)
    }
  })


  chemotherapy_description <- reactive({
    if (l2() > 0) {
      if (input$stage != "Design stage" && !is.null(df())) {
        chemotherapy_mapping <- df() %>%
          dplyr::select(chemotherapy, chemotherapy_description) %>%
          dplyr::arrange(chemotherapy, chemotherapy_description) %>%
          dplyr::group_by(chemotherapy, chemotherapy_description) %>%
          dplyr::slice(dplyr::n())
        chemotherapy_mapping$chemotherapy_description
      } else if (input$stage != "Design stage" && is.null(df())) {
        paste("Chemotherapy", 1:m())
      } else if (input$stage == "Design stage") {
        rownames(input[[paste0("chemotherapy_occupancy_", m())]])
      }
    } else {
      "None"
    }
  })



  observeEvent(list(treatment_description(), drug_name(), k(), l1()), {
    updateMatrixInput(
      session, paste0("treatment_by_drug_", k(), "_", l1()),
      value=matrix(treatment_by_drug(), ncol = l1(),
                   dimnames = dimnames(treatment_by_drug())))
  })


  observeEvent(list(chemotherapy_description(), drug_name(), m(), l2()), {
    if (l2() > 0) {
      updateMatrixInput(
        session, paste0("chemotherapy_by_drug_", m(), "_", l2()),
        value=matrix(chemotherapy_by_drug(), ncol = l2(),
                     dimnames = dimnames(chemotherapy_by_drug())))
    }
  })




  # whether to allow the user to specify the number of chemotherapies
  observeEvent(input$stage, {
    if (input$l2 > 0) {
      shinyjs::toggleState("m", input$stage == "Design stage")
    }
  })



  treatment_by_drug <- reactive({
    req(l())
    param = input[[paste0("treatment_by_drug_", k(), "_", l1())]]
    t = as.numeric(param)

    valid = all(t==1 | t==0)
    if (!valid) {
      showNotification(
        "Entries of treatment by drug matrix must be 1 or 0"
      )
    }

    req(valid)

    matrix(t, nrow = k(), ncol = l1(),
           dimnames = list(treatment_description(), drug_name()[1:l1()]))
  })


  treatment_by_drug_df <- reactive({
    req(l())
    param = input[[paste0("treatment_by_drug_", k(), "_", l1())]]
    t = as.numeric(param)

    valid = all(t==1 | t==0)
    if (!valid) {
      showNotification(
        "Entries of treatment by drug matrix must be 1 or 0"
      )
    }

    req(valid)

    dplyr::tibble(treatment = rep(1:k(), l1()),
                  drug = rep(1:l1(), each=k()),
                  drug_name = rep(drug_name()[1:l1()], each=k()),
                  dose_unit = rep(dose_unit()[1:l1()], each=k()),
                  included = as.logical(t)) %>%
      dplyr::filter(included) %>%
      dplyr::select(treatment, drug, drug_name, dose_unit)
  })



  chemotherapy_by_drug <- reactive({
    req(l2())
    if (l2() > 0) {
      param = input[[paste0("chemotherapy_by_drug_", m(), "_", l2())]]
      t = as.numeric(param)

      valid = all(t==1 | t==0)
      if (!valid) {
        showNotification(
          "Entries of chemotherapy by drug matrix must be 1 or 0"
        )
      }

      req(valid)

      matrix(t, nrow = m(), ncol = l2(),
             dimnames = list(chemotherapy_description(),
                             drug_name()[(l1()+1):l()]))
    } else {
      matrix(1, nrow = 1, ncol = 1, dimnames = list("None", "None"))
    }
  })


  chemotherapy_by_drug_df <- reactive({
    if (l2() > 0) {
      param = input[[paste0("chemotherapy_by_drug_", m(), "_", l2())]]
      t = as.numeric(param)

      valid = all(t==1 | t==0)
      if (!valid) {
        showNotification(
          "Entries of chemotherapy by drug matrix must be 1 or 0"
        )
      }

      req(valid)

      dplyr::tibble(chemotherapy = rep(1:m(), l2()),
                    drug = rep((l1()+1):l(), each=m()),
                    drug_name = rep(drug_name()[(l1()+1):l()], each=m()),
                    dose_unit = rep(dose_unit()[(l1()+1):l()], each=m()),
                    included = as.logical(t)) %>%
        dplyr::filter(included) %>%
        dplyr::select(chemotherapy, drug, drug_name, dose_unit)
    }
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

    if (input$l2 > 0) {
      required_columns <- c(required_columns, 'chemotherapy')
    }

    column_names <- colnames(df)

    shiny::validate(
      need(all(required_columns %in% column_names),
           "You don't have the right data"))

    if ('treatment' %in% column_names &&
        !('treatment_description' %in% column_names)) {
      df <- df %>% dplyr::mutate(
        treatment_description = paste0("Treatment ", treatment))
    }

    if ('chemotherapy' %in% column_names &&
        !('chemotherapy_description' %in% column_names)) {
      df <- df %>% dplyr::mutate(
        chemotherapy_description = paste0("Chemotherapy ", chemotherapy))
    }

    dplyr::tibble(df) %>%
      dplyr::mutate(trialsdt = as.Date(trialsdt),
                    randdt = as.Date(randdt),
                    cutoffdt = as.Date(cutoffdt))

  })


  # summarize observed data
  observed <- reactive({
    if (!is.null(df()))
      summarizeObserved(df(), to_predict(), showplot = FALSE,
                        input$by_treatment)
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
      fitEvent(df(), input$event_model, piecewiseSurvivalTime(),
               spline_k(), input$spline_scale, showplot = FALSE,
               input$by_treatment)
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
      w = treatment_allocation()/sum(treatment_allocation())

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
      if (to_predict() == "Enrollment and event") {
        model = input$event_prior
        event_prior <- list()

        for (i in 1:k()) {
          if (model == "Exponential") {
            theta = log(exponential_survival()[i])
          } else if (model == "Weibull") {
            theta = c(log(weibull_survival()[2,i]),
                      -log(weibull_survival()[1,i]))
          } else if (model == "Log-logistic") {
            theta = c(llogis_survival()[1,i], log(llogis_survival()[2,i]))
          } else if (model == "Log-normal") {
            theta = c(lnorm_survival()[1,i], log(lnorm_survival()[2,i]))
          } else if (model == "Piecewise exponential") {
            theta = log(piecewise_exponential_survival()[,i+1])
            piecewiseSurvivalTime = piecewise_exponential_survival()[,1]
          }

          if (model != "Piecewise exponential") {
            event_prior[[i]] <- list(
              model = model,
              theta = theta,
              vtheta = diag(length(theta))*1e-8,
              w = w[i])
          } else {
            event_prior[[i]] <- list(
              model = model,
              theta = theta,
              vtheta = diag(length(theta))*1e-8,
              piecewiseSurvivalTime = piecewiseSurvivalTime,
              w = w[i])
          }
        }

        if (k() == 1) event_prior <- event_prior[[1]]
      }

      # get prediction results based on what to predict
      if (to_predict() == "Enrollment only") {
        getPrediction(
          to_predict = to_predict(),
          target_n = target_n(),
          enroll_prior = enroll_prior,
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          ngroups = k(),
          alloc = treatment_allocation(),
          treatment_label = treatment_description())
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
          showplot = FALSE,
          by_treatment = input$by_treatment,
          ngroups = k(),
          alloc = treatment_allocation(),
          treatment_label = treatment_description())
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
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          alloc = treatment_allocation())
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
          showplot = FALSE,
          by_treatment = input$by_treatment,
          alloc = treatment_allocation())
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
          showplot = FALSE,
          by_treatment = input$by_treatment)
      }
    }
  })


  dosing <- reactive({
    if (predict_dosing()) {
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

      plower = (1 - pilevel())/2
      pupper = 1 - plower

      if (input$stage == 'Design stage') {
        newEvents <- pred()$event_pred$newEvents
        if (!("treatment" %in% colnames(newEvents))) {
          newEvents$treatment = 1
        }

        # generate chemotherapy for each subject
        if (l2() > 0) {
          chemoind <- t(rmultinom(nrow(newEvents), 1,
                                  chemotherapy_occupancy()))
          newEvents$chemotherapy = chemoind %*% seq(1,m())
          newEvents$chemotherapy_description =
            chemotherapy_description()[newEvents$chemotherapy]
        }

        # calculate cumulative doses
        t = pred()$event_pred$event_pred_df$t

        df2b = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents) %>%
          dplyr::filter(arrivalTime <= t)

        dfb <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df2b %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df2b %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t, draw) %>%
              dplyr::summarise(cum_dose = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N)),
                .groups = "drop_last")
          }))

        # predicted number of doses dispensed
        dosing_pred_df <- dfb %>%
          dplyr::group_by(drug, drug_name, dose_unit, t) %>%
          dplyr::summarise(n = quantile(cum_dose, probs = 0.5),
                           lower = quantile(cum_dose, probs = plower),
                           upper = quantile(cum_dose, probs = pupper),
                           mean = mean(cum_dose),
                           var = var(cum_dose),
                           .groups = "drop_last") %>%
          dplyr::ungroup()


        pred <- pred()
        pred$event_pred$newEvents <- newEvents
        pred$event_pred$dosing_pred_df <- dosing_pred_df

        if (subject_dosing()) {
          # subject-level dosing data
          dosing_df <- dplyr::bind_rows(
            lapply(1:l(), function(j) {
              schedule <- dosing_schedule_df() %>%
                dplyr::filter(drug == j)
              w = schedule$w
              d = schedule$d
              N = schedule$N

              if (j <= l1()) {
                dfx <- newEvents %>%
                  dplyr::inner_join(treatment_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "treatment")
              } else {
                dfx <- newEvents %>%
                  dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "chemotherapy")
              }

              dfx %>%
                dplyr::select(-c(event, dropout)) %>%
                dplyr::group_by(draw, usubjid, drug, drug_name, dose_unit) %>%
                dplyr::mutate(y = list(g_cum_dose(
                  min(time - 1, 365*nyears()), w, d, N))) %>%
                unnest(y) %>%
                dplyr::mutate(day_rel_trialsdt = arrivalTime +
                                day_rel_randdt - 1)
            })
          )

          pred$event_pred$dosing_df <- dosing_df
        }

        pred
      } else if (input$stage == 'Real-time before enrollment completion') {
        df0 <- dplyr::tibble(drug = 1:l(),
                             drug_name = drug_name(),
                             dose_unit = dose_unit(),
                             t = 1, n = 0, lower = NA, upper = NA,
                             mean = 0, var = 0)

        df <- df() %>%
          dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                        totalTime = arrivalTime + time - 1)

        trialsdt = pred()$observed$trialsdt
        cutoffdt = pred()$observed$cutoffdt
        t0 = pred()$observed$t0
        t1 = t0 + 365*nyears()

        # observed dosing before data cut
        t = unique(c(seq(1, t0, 30), t0))

        df1 <- dplyr::tibble(t = t) %>%
          dplyr::cross_join(df) %>%
          dplyr::filter(arrivalTime <= t)

        dosing_pred_obs_df <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df1 %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df1 %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t) %>%
              dplyr::summarise(n = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N)),
                .groups = "drop_last") %>%
              dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0)
          })) %>%
          dplyr::bind_rows(df0) %>%
          dplyr::group_by(drug, drug_name, dose_unit, t) %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::group_by(drug, drug_name, dose_unit)

        dosing_pred_t0 <- dosing_pred_obs_df %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::rename(cum_dose_t0 = n) %>%
          dplyr::select(drug, drug_name, dose_unit, cum_dose_t0)

        # predicted dosing after data cut
        newEvents <- pred()$event_pred$newEvents
        if (!("treatment" %in% colnames(newEvents))) {
          newEvents$treatment = 1
        }

        # obtain chemotherapy information
        if (l2() > 0) {
          # obtain chemotherapy info from observed data for ongoing patients
          ongoingSubjectsNewEvents <- newEvents %>%
            dplyr::inner_join(df %>% dplyr::select(usubjid, chemotherapy),
                              by = "usubjid")

          # obtain number of patients taking each chemotherapy in observed data
          s = as.numeric(table(factor(df$chemotherapy)))

          # draw chemotherapy occupancy probabilities from posterior
          p = rdirichlet(nreps(), s+1)

          # simulate chemotherapy for new patients
          n1 = target_n() - pred()$observed$n0
          usubjidNew <- paste0("Z-", 100000 + (1:n1))
          newSubjectsChemo <- rep(NA, nreps()*n1)
          for (i in 1:nreps()) {
            index = ((i-1)*n1+1):(i*n1)
            chemoind = t(rmultinom(n1, 1, p[i,]))
            newSubjectsChemo[index] <- chemoind %*% seq(1,m())
          }

          newSubjectsNewEvents <- newEvents %>%
            dplyr::right_join(dplyr::tibble(draw = rep(1:nreps(), each=n1),
                                            usubjid = rep(usubjidNew, nreps()),
                                            chemotherapy = newSubjectsChemo),
                              by = c("draw", "usubjid"))

          # combine ongoing and new patients and add chemotherapy description
          newEvents <- dplyr::bind_rows(ongoingSubjectsNewEvents,
                                        newSubjectsNewEvents) %>%
            dplyr::arrange(draw)

          newEvents$chemotherapy_description =
            chemotherapy_description()[newEvents$chemotherapy]
        }


        # calculate cumulative doses
        t = pred()$event_pred$event_pred_df$t
        t = unique(t[t >= t0])

        # from ongoing subjects
        df2a = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents) %>%
          dplyr::filter(arrivalTime <= t0 & totalTime > t0)

        dfa <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df2a %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df2a %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t, draw) %>%
              dplyr::summarise(cum_dose_a = sum(
                (f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N) -
                   f_cum_dose(t0 - arrivalTime, w, d, N))),
                .groups = "drop_last")
          }))


        # from new subjects
        df2b = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents) %>%
          dplyr::filter(arrivalTime > t0 & arrivalTime <= t)

        dfb <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df2b %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df2b %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t, draw) %>%
              dplyr::summarise(cum_dose_b = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N)),
                .groups = "drop_last")
          }))


        # dose dispense for time points after data cut
        dosing_pred_new_df <- dfa %>%
          dplyr::full_join(dfb, by = c("drug", "drug_name", "dose_unit",
                                       "t", "draw")) %>%
          dplyr::left_join(dosing_pred_t0, by = c('drug', 'drug_name',
                                                  'dose_unit')) %>%
          dplyr::mutate(cum_dose = cum_dose_t0 +
                          ifelse(is.na(cum_dose_a), 0, cum_dose_a) +
                          ifelse(is.na(cum_dose_b), 0, cum_dose_b)) %>%
          dplyr::group_by(drug, drug_name, dose_unit, t) %>%
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
          dplyr::arrange(drug, drug_name, dose_unit, t) %>%
          dplyr::ungroup()


        pred <- pred()
        pred$event_pred$newEvents <- newEvents
        pred$event_pred$dosing_pred_df <- dosing_pred_df

        if (subject_dosing()) {
          # subject-level dosing data
          discont_df <- df %>%
            dplyr::filter(event == 1) %>%
            dplyr::mutate(draw = 0)

          newEvents <- newEvents %>%
            dplyr::mutate(trialsdt = trialsdt,
                          cutoffdt = cutoffdt,
                          randdt = as.Date(arrivalTime - 1, origin = trialsdt))

          allsubjects <- dplyr::bind_rows(discont_df, newEvents) %>%
            dplyr::mutate(adt = as.Date(time - 1, origin = randdt))

          dosing_df <- dplyr::bind_rows(
            lapply(1:l(), function(j) {
              schedule <- dosing_schedule_df() %>%
                dplyr::filter(drug == j)
              w = schedule$w
              d = schedule$d
              N = schedule$N

              if (j <= l1()) {
                dfx <- allsubjects %>%
                  dplyr::inner_join(treatment_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "treatment")
              } else {
                dfx <- allsubjects %>%
                  dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "chemotherapy")
              }

              dfx %>%
                dplyr::select(-c(event, dropout)) %>%
                dplyr::group_by(draw, usubjid, drug, drug_name, dose_unit) %>%
                dplyr::mutate(y = list(g_cum_dose(
                  min(time - 1, t1 - arrivalTime), w, d, N))) %>%
                unnest(y) %>%
                dplyr::mutate(day_rel_trialsdt = arrivalTime +
                                day_rel_randdt - 1)
            })
          )

          dosing_df <- dosing_df %>%
            dplyr::mutate(date = as.Date(day_rel_trialsdt - 1,
                                         origin = trialsdt))

          pred$event_pred$dosing_df <- dosing_df
        }

        pred
      } else if (input$stage == "Real-time after enrollment completion") {
        df0 <- dplyr::tibble(drug = 1:l(),
                             drug_name = drug_name(),
                             dose_unit = dose_unit(),
                             t = 1, n = 0, lower = NA, upper = NA,
                             mean = 0, var = 0)

        df <- df() %>%
          dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                        totalTime = arrivalTime + time - 1)

        trialsdt = pred()$observed$trialsdt
        cutoffdt = pred()$observed$cutoffdt
        t0 = pred()$observed$t0
        t1 = t0 + 365*nyears()

        # observed dosing before data cut
        t = unique(c(seq(1, t0, 30), t0))

        df1 <- dplyr::tibble(t = t) %>%
          dplyr::cross_join(df) %>%
          dplyr::filter(arrivalTime <= t)

        dosing_pred_obs_df <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df1 %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df1 %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t) %>%
              dplyr::summarise(n = sum(
                f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N)),
                .groups = "drop_last") %>%
              dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0)
          })) %>%
          dplyr::bind_rows(df0) %>%
          dplyr::group_by(drug, drug_name, dose_unit, t) %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::group_by(drug, drug_name, dose_unit)

        dosing_pred_t0 <- dosing_pred_obs_df %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::rename(cum_dose_t0 = n) %>%
          dplyr::select(drug, drug_name, dose_unit, cum_dose_t0)

        # predicted dosing after data cut
        newEvents <- pred()$event_pred$newEvents
        if (!("treatment" %in% colnames(newEvents))) {
          newEvents$treatment = 1
        }

        # obtain chemotherapy info from observed data
        if (l2() > 0) {
          newEvents <- newEvents %>%
            dplyr::inner_join(df %>% dplyr::select(usubjid, chemotherapy,
                                                   chemotherapy_description),
                              by = "usubjid")
        }

        # calculate cumulative doses
        t = pred()$event_pred$event_pred_df$t
        t = unique(t[t >= t0])

        df2a = dplyr::tibble(t = t) %>%
          dplyr::cross_join(newEvents) %>%
          dplyr::filter(arrivalTime <= t0 & totalTime > t0)

        dfa <- dplyr::bind_rows(
          lapply(1:l(), function(j) {
            schedule <- dosing_schedule_df() %>%
              dplyr::filter(drug == j)
            w = schedule$w
            d = schedule$d
            N = schedule$N

            if (j <= l1()) {
              dfx <- df2a %>%
                dplyr::inner_join(treatment_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "treatment")
            } else {
              dfx <- df2a %>%
                dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                    dplyr::filter(drug == j),
                                  by = "chemotherapy")
            }

            dfx %>%
              dplyr::group_by(drug, drug_name, dose_unit, t, draw) %>%
              dplyr::summarise(cum_dose_a = sum(
                (f_cum_dose(pmin(totalTime, t) - arrivalTime, w, d, N) -
                   f_cum_dose(t0 - arrivalTime, w, d, N))),
                .groups = "drop_last")
          }))

        dosing_pred_new_df <- dfa %>%
          dplyr::left_join(dosing_pred_t0, by = c('drug', 'drug_name',
                                                  'dose_unit')) %>%
          dplyr::mutate(cum_dose = cum_dose_t0 +
                          ifelse(is.na(cum_dose_a), 0, cum_dose_a)) %>%
          dplyr::group_by(drug, drug_name, dose_unit, t) %>%
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
          dplyr::arrange(drug, drug_name, dose_unit, t) %>%
          dplyr::ungroup()


        pred <- pred()
        pred$event_pred$newEvents <- newEvents
        pred$event_pred$dosing_pred_df <- dosing_pred_df

        if (subject_dosing()) {
          # subject-level dosing data
          discont_df <- df %>%
            dplyr::filter(event == 1) %>%
            dplyr::mutate(draw = 0)

          newEvents <- newEvents %>%
            dplyr::mutate(trialsdt = trialsdt,
                          cutoffdt = cutoffdt,
                          randdt = as.Date(arrivalTime - 1, origin = trialsdt))

          allsubjects <- dplyr::bind_rows(discont_df, newEvents) %>%
            dplyr::mutate(adt = as.Date(time - 1, origin = randdt))

          dosing_df <- dplyr::bind_rows(
            lapply(1:l(), function(j) {
              schedule <- dosing_schedule_df() %>%
                dplyr::filter(drug == j)
              w = schedule$w
              d = schedule$d
              N = schedule$N

              if (j <= l1()) {
                dfx <- allsubjects %>%
                  dplyr::inner_join(treatment_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "treatment")
              } else {
                dfx <- allsubjects %>%
                  dplyr::inner_join(chemotherapy_by_drug_df() %>%
                                      dplyr::filter(drug == j),
                                    by = "chemotherapy")
              }

              dfx %>%
                dplyr::select(-c(event, dropout)) %>%
                dplyr::group_by(draw, usubjid, drug, drug_name, dose_unit) %>%
                dplyr::mutate(y = list(g_cum_dose(
                  min(time - 1, t1 - arrivalTime), w, d, N))) %>%
                unnest(y) %>%
                dplyr::mutate(day_rel_trialsdt = arrivalTime +
                                day_rel_randdt - 1)
            })
          )

          dosing_df <- dosing_df %>%
            dplyr::mutate(date = as.Date(day_rel_trialsdt - 1,
                                         origin = trialsdt))

          pred$event_pred$dosing_df <- dosing_df
        }

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

      if (input$by_treatment && k() > 1) {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          sum_by_trt <- df() %>%
            dplyr::bind_rows(df() %>% dplyr::mutate(
              treatment = 9999, treatment_description = "Overall")) %>%
            dplyr::group_by(treatment, treatment_description) %>%
            dplyr::summarise(n0 = n(),
                             d0 = sum(event),
                             r0 = sum(!event),
                             rp = sum((time < as.numeric(
                               cutoffdt - randdt + 1)) & !event),
                             .groups = "drop")

          if (any(sum_by_trt$rp) > 0) {
            table <- t(sum_by_trt %>% dplyr::select(n0, d0, r0, rp))
            colnames(table) <- sum_by_trt$treatment_description
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects",
                                 "  With ongoing date before cutoff")
          } else {
            table <- t(sum_by_trt %>% dplyr::select(n0, d0, r0))
            colnames(table) <- sum_by_trt$treatment_description
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects")
          }
        } else {
          sum_by_trt <- df() %>%
            dplyr::bind_rows(df() %>% dplyr::mutate(
              treatment = 9999, treatment_description = "Overall")) %>%
            dplyr::group_by(treatment, treatment_description) %>%
            dplyr::summarise(n0 = n(), .groups = "drop")

          table <- t(sum_by_trt %>% dplyr::select(n0))
          colnames(table) <- sum_by_trt$treatment_description
          rownames(table) <- c("Current number of subjects")
        }
      } else {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          sum_overall <- dplyr::tibble(n0 = observed()$n0,
                                       d0 = observed()$d0,
                                       r0 = observed()$r0,
                                       rp = observed()$rp)

          if (sum_overall$rp > 0) {
            table <- t(sum_overall %>% dplyr::select(n0, d0, r0, rp))
            colnames(table) <- "Overall"
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects",
                                 "  With ongoing date before cutoff")
          } else {
            table <- t(sum_overall %>% dplyr::select(n0, d0, r0))
            colnames(table) <- "Overall"
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects")
          }
        } else {
          table <- t(dplyr::tibble(n0 = observed()$n0))
          colnames(table) <- "Overall"
          rownames(table) <- c("Current number of subjects")
        }
      }

      print(table, quote=FALSE)
    }
  })


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


  output$input_df <- renderDataTable(
    df(), options = list(pageLength = 10)
  )


  output$enroll_fit <- renderPlotly({
    if (!is.null(enroll_fit())) enroll_fit()$enroll_fit_plot
  })


  output$event_fit1 <- renderPlotly({
    if (!is.null(event_fit())) event_fit()$event_fit_plot
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
      enroll_pred_df <- pred()$enroll_pred$enroll_pred_df
      if ((!input$by_treatment || k() == 1) ||
          ((input$by_treatment || input$stage == 'Design stage') &&
           k() > 1 &&
           length(table(enroll_pred_df$treatment)) == k() + 1)) {
        g1 <- enroll_pred_plot
      } else {
        g1 <- NULL
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

      if (showEnrollment())
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$enroll_pred_df)

      if (showEvent())
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$event_pred_df)

      if (showOngoing())
        dfs <- dfs %>% dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)


      dfs$parameter <- factor(dfs$parameter, levels = c(
        "Enrollment", "Event", "Ongoing"))


      if ((!input$by_treatment || k() == 1) &&
          !("treatment" %in% names(dfs))) { # overall
        if (input$stage != 'Design stage') {
          dfa <- dfs %>% dplyr::filter(is.na(lower))
          dfb <- dfs %>% dplyr::filter(!is.na(lower))

          g1 <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", fillcolor = ~parameter,
              line = list(width=0)) %>%
            plotly::add_lines(
              data = dfb, x = ~date, y = ~n, color = ~parameter,
              line = list(width=2)) %>%
            plotly::add_lines(
              data = dfa, x = ~date, y = ~n, color = ~parameter,
              line = list(shape="hv", width=2)) %>%
            plotly::add_lines(
              x = rep(observed()$cutoffdt, 2), y = range(dfs$n),
              name = "cutoff", line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              annotations = list(
                x = observed()$cutoffdt, y = 0, text = 'cutoff',
                xanchor = "left", yanchor = "bottom", font = list(size = 12),
                showarrow = FALSE),
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'))

          if (observed()$tp < observed()$t0) {
            g1 <- g1 %>%
              plotly::add_lines(
                x = rep(observed()$cutofftpdt, 2), y = range(dfs$n),
                name = "prediction start",
                line = list(dash="dash", color="grey"),
                showlegend = FALSE) %>%
              plotly::layout(
                annotations = list(
                  x = observed()$cutofftpdt, y = 0,
                  text = 'prediction start',
                  xanchor = "left", yanchor = "bottom",
                  font = list(size=12), showarrow = FALSE))
          }

          if (showEvent()) {
            g1 <- g1 %>%
              plotly::add_lines(
                x = range(dfs$date), y = rep(target_d(), 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
              plotly::layout(
                annotations = list(
                  x = 0.95, xref = "paper", y = target_d(),
                  text = 'target events', xanchor = "right",
                  yanchor = "bottom", font = list(size = 12),
                  showarrow = FALSE))
          }
        } else { # Design stage
          g1 <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfs, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", fillcolor = ~parameter,
              line = list(width=0)) %>%
            plotly::add_lines(
              data = dfs, x = ~t, y = ~n, color = ~parameter,
              line = list(width=2)) %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'))

          if (showEvent()) {
            g1 <- g1 %>%
              plotly::add_lines(
                x = range(dfs$t), y = rep(target_d(), 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
              plotly::layout(
                annotations = list(
                  x = 0.95, xref = "paper", y = target_d(),
                  text = 'target events', xanchor = "right",
                  yanchor = "bottom", font = list(size = 12),
                  showarrow = FALSE))
          }
        }
      } else if (((input$by_treatment || input$stage == 'Design stage') &&
                  k() > 1) &&
                 ("treatment" %in% names(dfs)) &&
                 (length(table(dfs$treatment)) == k() + 1)) { # by treatment

        if (input$stage != 'Design stage') {
          dfa <- dfs %>% dplyr::filter(is.na(lower))
          dfb <- dfs %>% dplyr::filter(!is.na(lower))

          g <- list()
          for (i in c(9999, 1:k())) {
            dfsi <- dfs %>% dplyr::filter(treatment == i)
            dfbi <- dfb %>% dplyr::filter(treatment == i)
            dfai <- dfa %>% dplyr::filter(treatment == i)

            g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
              plotly::add_ribbons(
                data = dfbi, x = ~date, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", fillcolor = ~parameter,
                line = list(width=0)) %>%
              plotly::add_lines(
                data = dfbi, x = ~date, y = ~n, color = ~parameter,
                line = list(width=2)) %>%
              plotly::add_lines(
                data = dfai, x = ~date, y = ~n, color = ~parameter,
                line = list(shape="hv", width=2)) %>%
              plotly::add_lines(
                x = rep(observed()$cutoffdt, 2), y = range(dfsi$n),
                name = "cutoff", line = list(dash="dash"),
                showlegend = FALSE) %>%
              plotly::layout(
                xaxis = list(title = "", zeroline = FALSE),
                yaxis = list(zeroline = FALSE),
                legend = list(x = 0, y = 1.05, yanchor = "bottom",
                              orientation = 'h')) %>%
              plotly::layout(
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))


            if (observed()$tp < observed()$t0) {
              g[[(i+1) %% 9999]] <- g[[(i+1) %% 9999]] %>%
                plotly::add_lines(
                  x = rep(observed()$cutofftpdt, 2), y = range(dfsi$n),
                  name = "prediction start",
                  line = list(dash="dash", color="grey"),
                  showlegend = FALSE)
            }


            if (i == 9999) {
              g[[1]] <- g[[1]] %>%
                plotly::layout(
                  annotations = list(
                    x = observed()$cutoffdt, y = 0, text = 'cutoff',
                    xanchor = "left", yanchor = "bottom",
                    font = list(size = 12), showarrow = FALSE))

              if (observed()$tp < observed()$t0) {
                g[[1]] <- g[[1]] %>%
                  plotly::layout(
                    annotations = list(
                      x = observed()$cutofftpdt, y = 0,
                      text = 'prediction start',
                      xanchor = "left", yanchor = "bottom",
                      font = list(size=12), showarrow = FALSE))
              }

              if (showEvent()) {
                g[[1]] <- g[[1]] %>%
                  plotly::add_lines(
                    x = range(dfsi$date), y = rep(target_d(), 2),
                    name = 'target events', showlegend = FALSE,
                    line = list(dash="dot",
                                color="rgba(128, 128, 128, 0.5")) %>%
                  plotly::layout(
                    annotations = list(
                      x = 0.95, xref = "paper", y = target_d(),
                      text = 'target events', xanchor = "right",
                      yanchor = "bottom", font = list(size = 12),
                      showarrow = FALSE))
              }
            }
          }

        } else {  # Design stage
          g <- list()

          for (i in c(9999, 1:k())) {
            dfsi <- dfs %>% dplyr::filter(treatment == i)

            g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
              plotly::add_ribbons(
                data = dfsi, x = ~t, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", fillcolor = ~parameter,
                line = list(width=0)) %>%
              plotly::add_lines(
                data = dfsi, x = ~t, y = ~n, color = ~parameter,
                line = list(width=2)) %>%
              plotly::layout(
                xaxis = list(title = "Days since trial start",
                             zeroline = FALSE),
                yaxis = list(zeroline = FALSE),
                legend = list(x = 0, y = 1.05, yanchor = "bottom",
                              orientation = 'h')) %>%
              plotly::layout(
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))


            if (i == 9999) {
              if (showEvent()) {
                g[[1]] <- g[[1]] %>%
                  plotly::add_lines(
                    x = range(dfsi$t), y = rep(target_d(), 2),
                    name = 'target events', showlegend = FALSE,
                    line = list(dash="dot",
                                color="rgba(128, 128, 128, 0.5")) %>%
                  plotly::layout(
                    annotations = list(
                      x = 0.95, xref = "paper", y = target_d(),
                      text = 'target events', xanchor = "right",
                      yanchor = "bottom", font = list(size = 12),
                      showarrow = FALSE))
              }
            }
          }
        }

        g1 <- plotly::subplot(g, nrows = k() + 1, margin = 0.05)
      } else {
        g1 <- NULL
      }

    }

    g1

  })


  output$pred_plot <- renderUI({
    if (input$by_treatment && k() > 1) {
      plotlyOutput("pred_plot1", height=250*(k()+1))
    } else {
      plotlyOutput("pred_plot1")
    }
  })


  output$downloadEventSummaryData <- downloadHandler(
    filename = function() {
      paste0("event_summary_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      if (to_predict() == "Enrollment only") {
        eventsummarydata <- pred()$enroll_pred$enroll_pred_df
      } else {
        eventsummarydata <- pred()$event_pred$enroll_pred_df %>%
          dplyr::bind_rows(pred()$event_pred$event_pred_df) %>%
          dplyr::bind_rows(pred()$event_pred$ongoing_pred_df)
      }
      writexl::write_xlsx(eventsummarydata, file)
    }
  )


  output$downloadEventSubjectData <- downloadHandler(
    filename = function() {
      paste0("event_subject_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      if (to_predict() == "Enrollment only") {
        eventsubjectdata <- pred()$enroll_pred$newSubjects
        if (input$stage != 'Design stage') {
          df <- df() %>%
            dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                          totalTime = arrivalTime + time - 1)

          if (input$by_treatment) {
            eventsubjectdata <- df %>%
              dplyr::mutate(draw = 0) %>%
              dplyr::select(draw, usubjid, arrivalTime,
                            treatment, treatment_description) %>%
              dplyr::bind_rows(eventsubjectdata)
          } else {
            eventsubjectdata <- df %>%
              dplyr::mutate(draw = 0) %>%
              dplyr::select(draw, usubjid, arrivalTime) %>%
              dplyr::bind_rows(eventsubjectdata)
          }
        }
      } else {
        eventsubjectdata <- pred()$event_pred$newEvents
        if (input$stage != 'Design stage') {
          df <- df() %>%
            dplyr::mutate(arrivalTime = as.numeric(randdt - trialsdt + 1),
                          totalTime = arrivalTime + time - 1)

          if (input$by_treatment) {
            eventsubjectdata <- df %>%
              dplyr::filter(event == 1 | dropout == 1) %>%
              dplyr::mutate(draw = 0) %>%
              dplyr::select(draw, usubjid, arrivalTime,
                            treatment, treatment_description,
                            time, event, dropout, totalTime) %>%
              dplyr::bind_rows(eventsubjectdata)
          } else {
            eventsubjectdata <- df %>%
              dplyr::filter(event == 1 | dropout == 1) %>%
              dplyr::mutate(draw = 0) %>%
              dplyr::select(draw, usubjid, arrivalTime,
                            time, event, dropout, totalTime) %>%
              dplyr::bind_rows(eventsubjectdata)
          }
        }
      }
      writexl::write_xlsx(eventsubjectdata, file)
    }
  )



  lapply(1:6, function(i) {
    pwexp <- paste0("piecewise_exponential_survival_", i)
    observeEvent(input[[paste0("add_piecewise_exponential_survival_", i)]], {
      a = matrix(as.numeric(input[[pwexp]]), ncol=ncol(input[[pwexp]]))
      b = matrix(a[nrow(a),], nrow=1)
      b[1,1] = b[1,1] + 1
      c = rbind(a, b)
      rownames(c) = paste("Interval", seq(1,nrow(c)))
      colnames(c) = colnames(input[[pwexp]])
      updateMatrixInput(session, pwexp, c)
    })
  })


  output$dosing_plot1 <- renderPlotly({
    if (predict_dosing()) {
      req(dosing()$event_pred)
      req(dosing()$stage == input$stage &&
            dosing()$to_predict == to_predict())

      dosing_pred_df <- dosing()$event_pred$dosing_pred_df

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

        g <- list()
        for (j in 1:l()) {
          dfs <- dplyr::filter(dosing_pred_df, drug == j)
          dfa <- dplyr::filter(dosing_pred_df, drug == j & is.na(lower))
          dfb <- dplyr::filter(dosing_pred_df, drug == j & !is.na(lower))

          g[[j]] <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval") %>%
            plotly::add_lines(
              data = dfb, x = ~date, y = ~n, name = "median prediction",
              line = list(width=2)) %>%
            plotly::add_lines(
              data = dfa, x = ~date, y = ~n, name = "observed",
              line = list(width=2)) %>%
            plotly::add_lines(
              x = rep(observed()$cutoffdt, 2), y = range(dfs$n),
              name = "cutoff", line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(title = paste0("Doses to dispense ",
                                          "(", dfs$dose_unit[1], ")"),
                           zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'),
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", dfs$drug_name[1], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))

          if (observed()$tp < observed()$t0) {
            g[[j]] <- g[[j]] %>%
              plotly::add_lines(
                x = rep(observed()$cutofftpdt, 2), y = range(dfs$n),
                name = "prediction start",
                line = list(dash="dash", color="grey"),
                showlegend = FALSE)
          }

          if (j==1) {
            g[[j]] <- g[[j]] %>%
              plotly::layout(
                annotations = list(
                  x = observed()$cutoffdt, y = 0, text = 'cutoff',
                  xanchor = "left", yanchor = "bottom",
                  font = list(size = 12), showarrow = FALSE))

            if (observed()$tp < observed()$t0) {
              g[[j]] <- g[[j]] %>%
                plotly::layout(
                  annotations = list(
                    x = observed()$cutofftpdt, y = 0,
                    text = 'prediction start',
                    xanchor = "left", yanchor = "bottom",
                    font = list(size=12), showarrow = FALSE))
            }
          }
        }
      } else {
        g <- list()
        for (j in 1:l()) {
          dfs <- dplyr::filter(dosing_pred_df, drug == j)

          g[[j]] <- plotly::plot_ly() %>%
            plotly::add_ribbons(
              data = dfs, x = ~t, ymin = ~lower, ymax = ~upper,
              name = "prediction interval",
              fill = "tonexty", line = list(width=0)) %>%
            plotly::add_lines(
              data = dfs, x = ~t, y = ~n, name = "median prediction",
              line = list(width=2)) %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(title = paste0("Doses to dispense ",
                                          "(", dfs$dose_unit[1], ")"),
                           zeroline = FALSE),
              legend = list(x = 0, y = 1.05, yanchor = "bottom",
                            orientation = 'h'),
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", dfs$drug_name[1], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))
        }
      }

      g1 <- plotly::subplot(g, nrows = l(),
                            titleX = TRUE, titleY = TRUE, margin = 0.1)
    } else {
      g1 <- NULL
    }

    g1
  })


  output$dosing_plot <- renderUI({
    if (l() > 1) {
      plotlyOutput("dosing_plot1", height=250*l())
    } else {
      plotlyOutput("dosing_plot1")
    }
  })


  output$downloadDosingSummaryData <- downloadHandler(
    filename = function() {
      paste0("dosing_summary_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      dosingsummarydata <- dosing()$event_pred$dosing_pred_df
      writexl::write_xlsx(dosingsummarydata, file)
    }
  )

  output$downloadDosingSubjectData <- downloadHandler(
    filename = function() {
      paste0("dosing_subject_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      dosingsubjectdata <- dosing()$event_pred$dosing_df
      write.csv(dosingsubjectdata, file)
    }
  )

  observeEvent(input$add_accrualTime, {
    a = matrix(as.numeric(input$accrualTime),
               ncol=ncol(input$accrualTime))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$accrualTime)
    updateMatrixInput(session, "accrualTime", c)
  })


  observeEvent(input$del_accrualTime, {
    if (nrow(input$accrualTime) >= 2) {
      a = matrix(as.numeric(input$accrualTime),
                 ncol=ncol(input$accrualTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
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
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$piecewise_poisson_rate)
    updateMatrixInput(session, "piecewise_poisson_rate", c)
  })


  observeEvent(input$del_piecewise_poisson_rate, {
    if (nrow(input$piecewise_poisson_rate) >= 2) {
      a = matrix(as.numeric(input$piecewise_poisson_rate),
                 ncol=ncol(input$piecewise_poisson_rate))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
      colnames(b) = colnames(input$piecewise_poisson_rate)
      updateMatrixInput(session, "piecewise_poisson_rate", b)
    }
  })


  observeEvent(input$add_piecewiseSurvivalTime, {
    a = matrix(as.numeric(input$piecewiseSurvivalTime),
               ncol=ncol(input$piecewiseSurvivalTime))
    b = matrix(a[nrow(a),] + 1, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$piecewiseSurvivalTime)
    updateMatrixInput(session, "piecewiseSurvivalTime", c)
  })


  observeEvent(input$del_piecewiseSurvivalTime, {
    if (nrow(input$piecewiseSurvivalTime) >= 2) {
      a = matrix(as.numeric(input$piecewiseSurvivalTime),
                 ncol=ncol(input$piecewiseSurvivalTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
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
      rownames(c) = paste("Interval", seq(1,nrow(c)))
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
        rownames(b) = paste("Interval", seq(1,nrow(b)))
        colnames(b) = colnames(input[[pwexp]])
        updateMatrixInput(session, pwexp, b)
      }
    })
  })


  observeEvent(l(), {
    lapply(1:l(), function(i) {
      output[[paste0("label_dosing_schedule_", l(), "_", i)]] <-
        renderText(paste("<b>Dosing schedule for", drug_name()[i], "</b>"))
    })
  })


  lapply(1:12, function(j) {
    lapply(1:j, function(i) {
      dosing <- paste0("dosing_schedule_", j, "_", i)
      observeEvent(input[[paste0("add_dosing_schedule_", j, "_", i)]], {
        a = matrix(as.numeric(input[[dosing]]), ncol=ncol(input[[dosing]]))
        b = matrix(a[nrow(a),], nrow=1)
        c = rbind(a, b)
        rownames(c) = c(rownames(input[[dosing]]), paste("Interval", nrow(c)))
        colnames(c) = colnames(input[[dosing]])
        updateMatrixInput(session, dosing, c)
      })
    })
  })


  lapply(1:12, function(j) {
    lapply(1:j, function(i) {
      dosing <- paste0("dosing_schedule_", j, "_", i)
      observeEvent(input[[paste0("del_dosing_schedule_", j, "_", i)]], {
        if (nrow(input[[dosing]]) >= 2) {
          a = matrix(as.numeric(input[[dosing]]), ncol=ncol(input[[dosing]]))
          b = matrix(a[-nrow(a),], ncol=ncol(a))
          rownames(b) = rownames(input[[dosing]])[1:nrow(b)]
          colnames(b) = colnames(input[[dosing]])
          updateMatrixInput(session, dosing, b)
        }
      })
    })
  })


  # save inputs
  output$saveInputs <- downloadHandler(
    filename = function() {
      paste0("inputs_", Sys.Date(), "_drug_supply.rds")
    },

    content = function(file) {
      x <- list(
        treatment_allocation = matrix(
          treatment_allocation(), ncol=1,
          dimnames = list(treatment_description(), "Size")),
        exponential_survival = matrix(
          exponential_survival(), nrow = 1,
          dimnames = list(NULL, treatment_description())),
        weibull_survival = matrix(
          weibull_survival(), nrow = 2,
          dimnames = list(c("Shape", "Scale"), treatment_description())),
        llogis_survival = matrix(
          llogis_survival(), nrow = 2,
          dimnames = list(c("Location on log scale", "Scale on log scale"),
                          treatment_description())),
        lnorm_survival = matrix(
          lnorm_survival(), nrow = 2,
          dimnames = list(c("Mean on log scale", "SD on log scale"),
                          treatment_description())),
        piecewise_exponential_survival = matrix(
          piecewise_exponential_survival(), ncol = k()+1,
          dimnames = list(paste("Interval",
                                1:nrow(piecewise_exponential_survival())),
                          c("Starting time", treatment_description()))),
        chemotherapy_occupancy = matrix(
          chemotherapy_occupancy(), ncol = 1,
          dimnames = list(chemotherapy_description(), "Probability")),
        drug_description = matrix(
          as.character(input[[paste0("drug_description_", l())]]), ncol = 2,
          dimnames = list(NULL, c("Drug Name", "Dose Unit"))),
        treatment_by_drug = treatment_by_drug(),
        chemotherapy_by_drug = chemotherapy_by_drug(),
        enroll_prior = input$enroll_prior,
        poisson_rate = poisson_rate(),
        mu = mu(),
        delta = delta(),
        piecewise_poisson_rate = piecewise_poisson_rate(),
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
          dimnames = list(paste("Interval",
                                1:length(piecewiseSurvivalTime())),
                          "Starting time")),
        spline_k = spline_k(),
        spline_scale = input$spline_scale,
        l = l(),
        l2 = l2(),
        m = m(),
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

      for (i in 1:l()) {
        a = paste0("dosing_schedule_", l(), "_", i)
        x[[a]] = input[[a]]
      }

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

    if ((x$stage == 'Design stage' ||
         (x$by_treatment &&
          x$stage != 'Real-time after enrollment completion')) && x$k > 1) {
      updateMatrixInput(
        session, paste0("treatment_allocation_", x$k),
        value=x$treatment_allocation)
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Exponential') {
      updateMatrixInput(
        session, paste0("exponential_survival_", x$k),
        value=x$exponential_survival)
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Weibull') {
      updateMatrixInput(
        session, paste0("weibull_survival_", x$k),
        value=x$weibull_survival)
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Log-logistic') {
      updateMatrixInput(
        session, paste0("llogis_survival_", x$k),
        value=x$llogis_survival)
    }

    if (x$stage == 'Design stage' && x$event_prior == 'Log-normal') {
      updateMatrixInput(
        session, paste0("lnorm_survival_", x$k),
        value=x$lnorm_survival)
    }

    if (x$stage == 'Design stage' &&
        x$event_prior == 'Piecewise exponential') {
      updateMatrixInput(
        session, paste0("piecewise_exponential_survival_", x$k),
        value=x$piecewise_exponential_survival)
    }

    if (x$predict_dosing) {
      updateSelectInput(session, "l", selected=x$l)

      if (x$l > 2) {
        choices <- c(0, 2:min((x$l-1), 6))
      } else {
        choices <- 0
      }

      updateSelectInput(session, "l2", choices = choices, selected = x$l2)

      if (x$l2 > 0) {
        updateSelectInput(session, "m", selected=x$m)

        if (x$stage == 'Design stage') {
          updateMatrixInput(
            session, paste0("chemotherapy_occupancy_", x$m),
            value=x$chemotherapy_occupancy)
        }
      }

      updateMatrixInput(
        session, paste0("drug_description_", x$l),
        value=x$drug_description)
      updateMatrixInput(
        session, paste0("treatment_by_drug_", x$k, "_", x$l - x$l2),
        value=x$treatment_by_drug)

      if (x$l2 > 0) {
        updateMatrixInput(
          session, paste0("chemotherapy_by_drug_", x$m, "_", x$l2),
          value=x$chemotherapy_by_drug)
      }

      lapply(1:x$l, function(i) {
        a = paste0("dosing_schedule_", x$l, "_", i)
        updateMatrixInput(
          session, a, value=x[[a]])
      })
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
          session, "piecewise_poisson_rate", value=x$piecewise_poisson_rate)
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
            session, "accrualTime", value=x$accrualTime)
        }
      }

      if ((x$stage == 'Real-time before enrollment completion' &&
           x$to_predict == 'Enrollment and event') ||
          x$stage == 'Real-time after enrollment completion') {

        updateRadioButtons(session, "event_model", selected=x$event_model)

        if (x$event_model == "Piecewise exponential") {
          updateMatrixInput(
            session, "piecewiseSurvivalTime", value=x$piecewiseSurvivalTime)
        } else if (x$event_model == "Spline") {
          updateNumericInput(session, "spline_k", value=x$spline_k)
          updateRadioButtons(session, "spline_scale", selected=x$spline_scale)
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


    if (x$stage == 'Design stage' || x$by_treatment) {
      updateSelectInput(session, "k", selected=x$k)
    }

    updateNumericInput(session, "nreps", value=x$nreps)
    updateNumericInput(session, "seed", value=x$seed)

  })
}

# Run the application
shinyApp(ui = ui, server = server)

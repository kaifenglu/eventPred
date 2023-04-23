library(magrittr)

uniroot(function(t) lrstat::accrual(
  time = t, accrualTime = seq(0, 9),
  accrualIntensity = c(26/9*seq(1, 9), 26),
  accrualDuration = 22) - 300,
  c(1, 22))

lrstat::lrsamplesize(
  beta = 0.2, kMax = 1,
  accrualTime = seq(0, 9),
  accrualIntensity = c(26/9*seq(1, 9), 26),
  piecewiseSurvivalTime = c(0, 6),
  lambda1 = c(0.0533, 0.0309),
  lambda2 = c(0.0533, 0.0533),
  gamma1 = -log(1-0.05)/12,
  gamma2 = -log(1-0.05)/12,
  accrualDuration = 15.54,
  followupTime = NA, fixedFollowup = FALSE)

## code to prepare `observedData` dataset goes here
set.seed(123)
n = 300
trialsdt = as.Date("2018-03-01")

# piecewise accrual with 8-month ramp-up and 26 patients per month thereafter
accrualTime = seq(0, 8)*30.4375
accrualIntensity = 1/30.4375*26/9*seq(1, 9)

# general arrival time
rhs = cumsum(-log(runif(n)))
J = length(accrualIntensity)
psum = c(0, cumsum(accrualIntensity[1:(J-1)] * diff(accrualTime)))
j1 = findInterval(rhs, psum)
arrivalTime = ceiling(accrualTime[j1] +
                        (rhs - psum[j1])/accrualIntensity[j1])

# underlying time to event and time to dropout
k = 2
alloc = c(2,2)
blocksize = sum(alloc)
nblocks = ceiling(n/blocksize)
treats = rep(1:k, alloc)
treatment = c(replicate(nblocks, sample(treats)))[1:n]

n1 = sum(treatment == 1)
n2 = sum(treatment == 2)

# piecewise exponential distributions for time to event
u = c(0, 6)*30.4375 # left end points of time intervals
lambda1 = c(0.0533, 0.0309)/30.4375 # hazard rates for the treatment group
lambda0 = c(0.0533, 0.0533)/30.4375 # hazard rates for the control group

# function to generate the survival time for each treatment group
fsurvtime <- function(n, u, lambda) {
  J = length(u)
  psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
  rhs = rexp(n)
  j1 = findInterval(rhs, psum)
  survivalTime = u[j1] + (rhs - psum[j1])/lambda[j1]
  survivalTime
}

survivalTime = rep(NA, n)
survivalTime[treatment == 1] = fsurvtime(n1, u, lambda1)
survivalTime[treatment == 2] = fsurvtime(n2, u, lambda0)
survivalTime = ceiling(survivalTime)

# exponential distribution for time to dropout
dropoutTime = ceiling(rexp(n, rate = -log(1-0.05)/365))

# observed survival time and event indicator
time = pmin(survivalTime, dropoutTime)
event = 1*(time == survivalTime)
dropout = 1*(dropoutTime < survivalTime)

# complete data without administrative censoring
dfcomplete <- dplyr::tibble(trialsdt, arrivalTime, treatment,
                            time, event, dropout) %>%
  dplyr::mutate(totalTime = arrivalTime + time) %>%
  dplyr::arrange(totalTime) %>%
  dplyr::mutate(n = cumsum(event))

# end the study when the target number of events is reached
cutoff = floor(max((dfcomplete %>%
                      dplyr::filter(n == 244))$totalTime))


# complete data with administrative censoring
finalData <- dfcomplete %>%
  dplyr::filter(arrivalTime <= cutoff) %>%
  dplyr::mutate(cutoff = cutoff,
                followupTime = cutoff - arrivalTime,
                event = ifelse(time <= followupTime, event, 0),
                dropout = ifelse(time <= followupTime, dropout, 0),
                time = pmax(pmin(time, followupTime), 1e-8),
                randdt = as.Date(arrivalTime - 1, origin = trialsdt),
                cutoffdt = as.Date(cutoff - 1, origin = trialsdt)) %>%
  dplyr::select(trialsdt, randdt, cutoffdt, treatment,
                time, event, dropout) %>%
  dplyr::arrange(randdt)

trialsdt = min(finalData$randdt)

# partial data before enrollment completion

# time when 75% subjects are enrolled
cutoffdt = (finalData %>%
              dplyr::filter(dplyr::row_number() == 225))$randdt


interimData1 <- finalData %>%
  dplyr::select(-cutoffdt) %>%
  dplyr::filter(randdt <= cutoffdt) %>%
  dplyr::mutate(cutoffdt = cutoffdt,
                cutoff = as.numeric(cutoffdt - trialsdt + 1),
                arrivalTime = as.numeric(randdt - trialsdt + 1),
                followupTime = cutoff - arrivalTime,
                event = ifelse(time <= followupTime, event, 0),
                dropout = ifelse(time <= followupTime, dropout, 0),
                time = pmax(pmin(time, followupTime), 1e-8)) %>%
  dplyr::select(trialsdt, randdt, cutoffdt, treatment, time, event, dropout)


# partial data after enrollment completion

# time when about 85% events are observed
cutoffdt = max((finalData %>%
                  dplyr::mutate(adt = as.Date(time, origin = randdt)) %>%
                  dplyr::arrange(adt) %>%
                  dplyr::mutate(n = cumsum(event)) %>%
                  dplyr::filter(n == 183))$adt)
cutoff = as.numeric(cutoffdt - trialsdt + 1)

interimData2 <- finalData %>%
  dplyr::select(-cutoffdt) %>%
  dplyr::filter(randdt <= cutoffdt) %>%
  dplyr::mutate(cutoffdt = cutoffdt,
                cutoff = as.numeric(cutoffdt - trialsdt + 1),
                arrivalTime = as.numeric(randdt - trialsdt + 1),
                followupTime = cutoff - arrivalTime,
                event = ifelse(time <= followupTime, event, 0),
                dropout = ifelse(time <= followupTime, dropout, 0),
                time = pmax(pmin(time, followupTime), 1e-8)) %>%
  dplyr::select(trialsdt, randdt, cutoffdt, treatment, time, event, dropout)


# save to data/ folder
usethis::use_data(interimData1, overwrite = TRUE)
usethis::use_data(interimData2, overwrite = TRUE)
usethis::use_data(finalData, overwrite = TRUE)

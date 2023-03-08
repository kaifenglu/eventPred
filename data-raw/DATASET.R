## code to prepare `observedData` dataset goes here
n0 <- 300
set.seed(123)
trialsdt = as.Date("2021-03-01")

# piecewise accrual with 8-month ramp-up and 26 patients per month thereafter
accrualTime = seq(0, 8)*30.4375
accrualIntensity = 1/30.4375*26/9*seq(1, 9)

# general arrival time
rhs = cumsum(-log(runif(n0)))
J = length(accrualIntensity)
psum = c(0, cumsum(accrualIntensity[1:(J-1)] * diff(accrualTime)))
j1 = findInterval(rhs, psum)
arrivalTime = ceiling(accrualTime[j1] +
                        (rhs - psum[j1])/accrualIntensity[j1])

# time of data cut
analysisTime <- max(arrivalTime) + 1

# potential follow-up time for each subject
followupTime <- analysisTime - arrivalTime

# underlying time to event and time to dropout
treatmentGroup <- rbinom(n0, 1, 0.5)
rate <- (0.0309*(treatmentGroup == 1) +
           0.0533*(treatmentGroup == 0))/30.4375
survivalTime <- ceiling(rexp(n0, rate=rate))
dropoutTime <- ceiling(rexp(n0, rate=-log(1-0.05)/365))

# observed time to event or dropout, and event and dropout indicators
time <- pmin(survivalTime, dropoutTime, followupTime)
event <- 1*(time == survivalTime)
dropout <- 1*(time == dropoutTime)

# observed data
observedData <- data.frame(
  RANDDT = as.Date(arrivalTime - 1, origin = trialsdt),
  CUTOFFDT = as.Date(analysisTime - 1, origin = trialsdt),
  ADT = as.Date(arrivalTime + time - 1, origin = trialsdt),
  event = event, dropout = dropout)

# save to data/ folder
usethis::use_data(observedData, overwrite = TRUE)

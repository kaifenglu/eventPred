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

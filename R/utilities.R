# determine the placement of major ticks of x-axis
fbw <- function(n_months) {
  if (n_months <= 12) {
    bw = '1 month'
  } else if (n_months <= 24) {
    bw = '2 months'
  } else if (n_months <= 36) {
    bw = '3 months'
  } else if (n_months <= 72) {
    bw = '6 months'
  } else {
    bw = '1 year'
  }

  bw
}


# generate year labels
flabel <- function(df, trialsdt) {
  # set up labels for year part of the dates
  yrlim = range(df$year)

  # left end point of each year is Jan 1 except for the first year
  # right end point of each year is Dec 31 except for the last year
  label_range <- df %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarise(xmin = min(.data$date), xmax = max(.data$date)) %>%
    dplyr::mutate(xmin = as.Date(ifelse(.data$year == yrlim[1], .data$xmin,
                                        as.Date(paste0(.data$year, "-01-01"))),
                                 origin = "1970-01-01"),
                  xmax = as.Date(ifelse(.data$year == yrlim[2], .data$xmax,
                                        as.Date(paste0(.data$year, "-12-31"))),
                                 "1970-01-01"),
                  ymin = 0,
                  ymax = 1)

  # center the year label within each rectangle
  label_range <- label_range %>%
    dplyr::mutate(x1 = as.numeric(.data$xmin - trialsdt),
                  x2 = as.numeric(.data$xmax - trialsdt),
                  x = as.Date((.data$x1 + .data$x2)/2 + trialsdt))

  # year labels
  g2 <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = label_range, fill = "lightcoral",
                       color = "#f2f2f2",
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    group = .data$year)) +
    ggplot2::geom_text(data = label_range,
                       ggplot2::aes(x = .data$x, y = 0.5,
                                    group = .data$year, label = .data$year)) +
    ggplot2::scale_x_date(expand = c(0.01, 0.01)) +
    ggplot2::theme_void()

  g2
}


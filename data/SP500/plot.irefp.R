plot.irefp <- function(x, functional = "max", type = c("CUSUM", "MOSUM")) {
  type <- match.arg(type)
  x <- c(x$efpprocess, x$process)

  if(type == "CUSUM") {
    bound <- cus.bound
    x <- irts(stock.prices$time, x)
    y.at <- as.POSIXct(c("2001-08-01", "2001-09-01", "2001-10-01", "2001-11-01", "2001-12-01", "2002-01-01"))
  } else {
    bound <- mos.bound
    x <- irts(stock.prices$time[-(1:14)], x)
    y.at <- as.POSIXct(c("2001-09-01", "2001-10-01", "2001-11-01", "2001-12-01", "2002-01-01"))
  }

  if(is.null(functional)) {
    ylim <- c(min(-bound$value, x$value), max(bound$value, x$value))
  } else {
    x$value <- abs(x$value)
    ylim <- c(0, max(bound$value, x$value))
  }

  plot(x, ylim = ylim, ylab = "empirical fluctuation process", xlab = "Time", xaxt = "n", main = "")
  axis.POSIXct(1, at = y.at)
  lines(bound, col = 2)
  if(is.null(functional)) {
    bound$value <- -bound$value
    lines(bound, col = 2)
  }
  abline(h = 0)
  abline(v = as.POSIXct("2001-09-10"), lty = 2)
}

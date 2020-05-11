scatterplot.density.out <- function (x, y, zlim, xylim, num.bins = 64, col = kristen.colors(32), 
                                     xlab, ylab, main, density.in.percent = TRUE, col.regression.line = 1, 
                                     col.one.to.one.line = grey(0.4), col.bar.legend = TRUE, plt.beyond.zlim = FALSE, ...) {
  if (!is.numeric(x)) 
    stop("x must be numeric")
  if (!is.numeric(y)) 
    stop("y must be numeric")
  x.data <- as.vector(x)
  y.data <- as.vector(y)
  if (length(x.data) != length(y.data)) 
    stop("x and y must have the same length.")
  if (is.list(num.bins)) {
    num.bins.x <- num.bins$num.bins.x
    num.bins.y <- num.bins$num.bins.y
  }
  else {
    num.bins.x <- num.bins
    num.bins.y <- num.bins
  }
  if (missing(xylim)) {
    xlim <- range(x.data)
    ylim <- range(y.data)
  }
  else {
    if (is.list(xylim)) {
      xlim <- xylim$xlim
      ylim <- xylim$ylim
    }
    else {
      xlim <- xylim
      ylim <- xylim
    }
  }
  data.bins.x <- seq(xlim[1], xlim[2], length = num.bins.x)
  bin.x.length <- data.bins.x[2] - data.bins.x[1]
  plot.seq.x <- seq(xlim[1] - (bin.x.length/2), xlim[2] + (bin.x.length/2), 
                    length = num.bins.x + 1)
  data.bins.y <- seq(ylim[1], ylim[2], length = num.bins.y)
  bin.y.length <- data.bins.y[2] - data.bins.y[1]
  plot.seq.y <- seq(ylim[1] - (bin.y.length/2), ylim[2] + (bin.y.length/2), 
                    length = num.bins.y + 1)
  x.cut <- cut(x.data, plot.seq.x)
  y.cut <- cut(y.data, plot.seq.y)
  tab.x.y <- table(x.cut, y.cut)
  if (density.in.percent) 
    tab.x.y <- tab.x.y/length(x.data) * 100
  tab.x.y[tab.x.y == 0] <- NA
  if (missing(xlab)) 
    xlab <- deparse(substitute(x))
  if (missing(ylab)) 
    ylab <- deparse(substitute(y))
  if (missing(main)) {
    if (density.in.percent) 
      main <- "Data Density Plot (%)"
    else main <- "Data Frequency Plot (counts)"
  }
  if (plt.beyond.zlim) {
    if (missing(zlim)) {
      warning("plt.beyond.zlim=TRUE is not a valid option if zlim argument is not provided, changing to plt.beyond.zlim=FALSE")
      plt.beyond.zlim <- FALSE
    }
    else {
      tab.x.y[tab.x.y < zlim[1]] <- zlim[1]
      tab.x.y[tab.x.y > zlim[2]] <- zlim[2]
    }
  }
  if (missing(zlim)) 
    zlim <- range(tab.x.y, na.rm = T)
  image(x = plot.seq.x, y = plot.seq.y, z = tab.x.y, zlim = zlim, 
        col = col, xlab = xlab, ylab = ylab, main = main, ...)
  if (!is.null(col.one.to.one.line)) 
    abline(0, 1, col = col.one.to.one.line, lty = 3)
  if (!is.null(col.regression.line)) {
    y.x.lm <- lm(y.data ~ x.data)$coeff
    abline(y.x.lm, col = col.regression.line)
    legend.txt <- bquote(paste(hat(y), "=", .(signif(y.x.lm[1], 
                                                     2)), "+", .(signif(y.x.lm[2], 2)), "x"))
    legend("topleft", legend = do.call("expression", c(legend.txt)), 
           bty = "n", text.col = col.regression.line)
  }
  if (col.bar.legend) 
    vertical.image.legend(col = col, zlim = zlim)
  return(list(x = plot.seq.x, y = plot.seq.y, z = tab.x.y))
}



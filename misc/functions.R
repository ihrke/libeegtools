headline <- function(h){
  cat("+-----------------------------------------------\n");
  cat("| ",h,"\n");
  cat("+-----------------------------------------------\n");
}

subheadline <- function(h){
  cat(" ",h,"\n");
  cat("+-----------------------------------------------\n");
}
printf <- function(fmt,...){
  cat(sprintf(fmt, ...));
}

errbarplot <- function (x, y = NULL, uiw, liw = uiw, add=F, ..., sfrac = 0.01)  {
  if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (is.null(y)) {
    if (is.null(x)) 
      stop("both x and y NULL")
    y <- as.numeric(x)
    x <- seq(along = x)
  }
  ui <- y + uiw;
  li <- y - liw;
  if(add)
    points(x, y, ylim = range(c(y, ui, li)), ...)
  else
    plot(x, y,...)# ylim = range(c(y, ui, li)), ...)

  smidge <- diff(par("usr")[1:2]) * sfrac
  segments(x, li, x, ui,...)
  x2 <- c(x, x)
  ul <- c(li, ui)
  segments(x2 - smidge, ul, x2 + smidge, ul, ...)
  invisible(list(x = x, y = y))
}

sem <- function(y){
  return(sd(y)/(sqrt(length(y))));
}

cohensd <- function(v1, v2){
  # Effektstaerke als d = M1-M2/s_pooled mit
  # s_pooled = sqrt(1/2 * (s_1^2+s_2^2))
  d = abs(mean(v1)-mean(v2));
  d = d/sqrt((sd(v1)^2 + sd(v2)^2)/2);
  return(d);
};

cohensd_abs<- function(v1, v2){
  # Effektstaerke als d = M1-M2/s_pooled mit
  # s_pooled = sqrt(1/2 * (s_1^2+s_2^2))
  d = mean(abs(v1-v2));
  d = d/sqrt((sd(v1)^2 + sd(v2)^2)/2);
  return(d);
};


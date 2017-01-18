# These functions compute confidence intervals around the d'^ corresponding to
# the observed data according to Gourevitch and Galanter's (1967) approximation.
# There is no adjustment for perfect hit or false alarm rates.  The functions
# compute 1 - a confidence intervals, where a is specified in an additional
# argument; the default is 0.05.  The functions return a numeric vector of the
# lower and upper bounds of the confidence interval.

G.int.count <- function(hits, misses, false_alarms, correct_rejections, a = 0.05) {
  # get hits, 'signal' trials, and hit rate
  h = hits
  S = hits + misses
  ph = h / S
  # get false alarms, 'noise' trials, and false alarm rate
  f = false_alarms
  N = false_alarms + correct_rejections
  pf = f / N
  
  # compute d'
  dprime = qnorm(ph) - qnorm(pf)
  # Gourevitch and Galanter's equations
  sd = (((ph * (1 - ph)) / (S * (dnorm(qnorm(ph)) ^ 2))) + ((pf * (1 - pf)) / (N * (dnorm(qnorm(pf)) ^ 2)))) ^ .5
  int = c(dprime + (qnorm(a / 2) * sd), dprime + (qnorm(1 - (a / 2)) * sd))
  return(int)
}

G.int.list <- function(responses, a = 0.05) {
  # get hits, 'signal' trials, and hit rate
  h = length(responses[responses == "H"])
  S = length(responses[responses == "H" | responses == "M"])
  ph = h / S
  # get false alarms, 'noise' trials, and false alarm rate
  f = length(responses[responses == "F"])
  N = length(responses[responses == "F" | responses == "R"])
  pf = f / N
  
  # compute d'
  dprime = qnorm(ph) - qnorm(pf)
  # Gourevitch and Galanter's equations
  sd = (((ph * (1 - ph)) / (S * (dnorm(qnorm(ph)) ^ 2))) + ((pf * (1 - pf)) / (N * (dnorm(qnorm(pf)) ^ 2)))) ^ .5
  int = c(dprime + (qnorm(a / 2) * sd), dprime + (qnorm(1 - (a / 2)) * sd))
  return(int)
}

# These functions compute confidence intervals around the d'^ corresponding to
# the observed data according to Miller's (1996) method.  There is no adjustment
# for perfect hit or false alarm rates.  The functions compute 1 - a confidence
# intervals, where a is specified in an additional argument; the default is
# 0.05.  The functions return a numeric vector of the lower and upper bounds of
# the confidence interval.

pmf.int.count <- function(hits, misses, false_alarms, correct_rejections, a = 0.05) {
  # get number of 'signal' and 'noise' trials
  S = hits + misses
  N = false_alarms + correct_rejections
  
  # set bounds of confidence interval at observed d'
  observed_dprime = dprime.count(hits, misses, false_alarms, correct_rejections)
  lower_bound = observed_dprime
  upper_bound = observed_dprime
  
  # assume that the upper bound of the confidence interval is between the
  # observed d' and 10
  max_upper_bound = observed_dprime + 10
  min_upper_bound = observed_dprime
  current_upper_bound = observed_dprime + 5
  # if we have not identified the upper bound to an accuracy of at least 5
  # decimal places,
  while((round(current_upper_bound, 5) != round(max_upper_bound, 5)) & (round(current_upper_bound, 5) != round(min_upper_bound, 5))) {
    # find ph and pf that yield the current upper bound where ph = 1 - pf
    ph = pnorm(current_upper_bound / 2)
    pf = 1 - ph
    # calculate pmf of current upper bound above observed d'
    pmf = dprime.pmf(ph, S, pf, N, observed_dprime, "max")
    # probability of observed d' is equal to sum of pmf of current upper bound
    # above observed d'
    percentile = sum(pmf[,2])
    # if current upper bound is too low,
    if(percentile > (a / 2)) {
      # new upper bound to test: halfway between current upper bound and maximum
      # upper bound
      min_upper_bound = current_upper_bound
      # narrow search space appopriately
      current_upper_bound = min_upper_bound + ((max_upper_bound - min_upper_bound) / 2)
    }
    # if current upper bound is too high,
    else {
      # new upper bound to test: halfway between current upper bound and minimum
      # upper bound
      max_upper_bound = current_upper_bound
      # narrow search space appropriately
      current_upper_bound = min_upper_bound + ((max_upper_bound - min_upper_bound) / 2)
    }
  }
  # last identified uppper bound must be correct
  upper_bound = current_upper_bound
  
  # analogous procedure to identify lower bound
  max_lower_bound = observed_dprime
  min_lower_bound = observed_dprime - 10
  current_lower_bound = observed_dprime - 5
  while((round(current_lower_bound, 5) != round(max_lower_bound, 5)) & (round(current_lower_bound, 5) != round(min_lower_bound, 5))) {
    ph = pnorm(current_lower_bound / 2)
    pf = 1 - ph
    pmf = dprime.pmf(ph, S, pf, N, observed_dprime, "min")
    percentile = sum(pmf[,2])
    if(percentile < (a / 2)) {
      min_lower_bound = current_lower_bound
      current_lower_bound = min_lower_bound + ((max_lower_bound - min_lower_bound) / 2)
    }
    else {
      max_lower_bound = current_lower_bound
      current_lower_bound = min_lower_bound + ((max_lower_bound - min_lower_bound) / 2)
    }
  }
  lower_bound = current_lower_bound 
  
  return(c(lower_bound, upper_bound))
}

pmf.int.list <- function(responses, a = 0.05) {
  # get number of 'signal' and 'noise' trials
  S = length(responses[responses == "H" || responses == "M"])
  N = length(responses[responses == "F" || responses == "R"])
  
  # set bounds of confidence interval at observed d'
  observed_dprime = dprime.list(responses)
  lower_bound = observed_dprime
  upper_bound = observed_dprime
  
  # assume that the upper bound of the confidence interval is between the
  # observed d' and 10
  max_upper_bound = observed_dprime + 10
  min_upper_bound = observed_dprime
  current_upper_bound = observed_dprime + 5
  # if we have not identified the upper bound to an accuracy of at least 5
  # decimal places,
  while((round(current_upper_bound, 5) != round(max_upper_bound, 5)) & (round(current_upper_bound, 5) != round(min_upper_bound, 5))) {
    # find ph and pf that yield the current upper bound where ph = 1 - pf
    ph = pnorm(current_upper_bound / 2)
    pf = 1 - ph
    # calculate pmf of current upper bound above observed d'
    pmf = dprime.pmf(ph, S, pf, N, observed_dprime, "max")
    # probability of observed d' is equal to sum of pmf of current upper bound
    # above observed d'
    percentile = sum(pmf[,2])
    # if current upper bound is too low,
    if(percentile > (a / 2)) {
      # new upper bound to test: halfway between current upper bound and maximum
      # upper bound
      min_upper_bound = current_upper_bound
      # narrow search space appopriately
      current_upper_bound = min_upper_bound + ((max_upper_bound - min_upper_bound) / 2)
    }
    # if current upper bound is too high,
    else {
      # new upper bound to test: halfway between current upper bound and minimum
      # upper bound
      max_upper_bound = current_upper_bound
      # narrow search space appropriately
      current_upper_bound = min_upper_bound + ((max_upper_bound - min_upper_bound) / 2)
    }
  }
  # last identified uppper bound must be correct
  upper_bound = current_upper_bound
  
  # analogous procedure to identify lower bound
  max_lower_bound = observed_dprime
  min_lower_bound = observed_dprime - 10
  current_lower_bound = observed_dprime - 5
  while((round(current_lower_bound, 5) != round(max_lower_bound, 5)) & (round(current_lower_bound, 5) != round(min_lower_bound, 5))) {
    ph = pnorm(current_lower_bound / 2)
    pf = 1 - ph
    pmf = dprime.pmf(ph, S, pf, N, observed_dprime, "min")
    percentile = sum(pmf[,2])
    if(percentile < (a / 2)) {
      min_lower_bound = current_lower_bound
      current_lower_bound = min_lower_bound + ((max_lower_bound - min_lower_bound) / 2)
    }
    else {
      max_lower_bound = current_lower_bound
      current_lower_bound = min_lower_bound + ((max_lower_bound - min_lower_bound) / 2)
    }
  }
  lower_bound = current_lower_bound 
  
  return(c(lower_bound, upper_bound))
}

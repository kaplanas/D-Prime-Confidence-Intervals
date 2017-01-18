# This function computes the sampling distribution of a given d' (that is, its
# pmf); it takes as arguments values for ph, nS, pf, and nN.  It returns a two-
# column matrix where the first column contains an ordered list of possible d'^s
# and the second column their corresponding probabilities.
#
# If the last argument is specified as "max" or "min", the function returns
# values only for d'^s above (or below) a certain value.  The limit is specified# in the second-to-last argument.  Using this option makes Miller's method
# slightly more efficient.

dprime.pmf <- function(p_hit, signal_trials, p_false_alarm, noise_trials, limit = 0, max_or_min = "neither") {
  # rename arguments for compactness
  ph = p_hit
  S = signal_trials
  pf = p_false_alarm
  N = noise_trials
  
  # store possible observed hit and false alarm rates
  Hs = 0:S
  Fs = 0:N
  # create variables to hold calculated values based on possible hit and false
  # alarm rates
  pmfHs = rep(0, length(Hs))
  pmfFs = rep(0, length(Fs))
  qnormHs = rep(0, length(Hs))
  qnormFs = rep(0, length(Fs))
  ps = rep(0, length(Hs) * length(Fs))
  dim(ps) = c(length(Hs), length(Fs))
  dprimes = rep(0, length(Hs) * length(Fs))
  dim(dprimes) = c(length(Hs), length(Fs))
  
  # for each possible hit rate,
  for(i in 1:length(Hs)) {
    # compute the probability of observing that hit rate
    pmfHs[i] = choose(S, Hs[i]) * (ph ^ Hs[i]) * ((1 - ph) ^ (S - Hs[i]))
    # adjust hit rate if necessary
    adjusted_H = Hs[i]
    if(adjusted_H == 0) { adjusted_H = .5 }
    else if(adjusted_H == S) { adjusted_H = S - .5 }
    # compute inverse-normal transform of hit rate
    qnormHs[i] = qnorm(adjusted_H / S)
  }
  # analogous calculations for each possible false alarm rate
  for(j in 1:length(Fs)) {
    pmfFs[j] = choose(N, Fs[j]) * (pf ^ Fs[j]) * ((1 - pf) ^ (N - Fs[j]))
    adjusted_F = Fs[j]
    if(adjusted_F == 0) { adjusted_F = .5 }
    else if(adjusted_F == N) { adjusted_F = N - .5 }
    qnormFs[j] = qnorm(adjusted_F / N)
  }
  # for each possible combination of hit and false alarm rates,
  for(i in 1:length(Hs)) {
    for(j in 1:length(Fs)) {
      # compute the probability of observing that combination
      ps[i,j] = pmfHs[i] * pmfFs[j]
      # compute the d' for that combination
      dprimes[i,j] = qnormHs[i] - qnormFs[j]
    }
  }
  
  # sort possible d's and eliminate duplicates
  sorted_dprimes = unique(sort(dprimes))
  # eliminate d's above maximum or below minimum, if necessary
  if(max_or_min == "max") {
    sorted_dprimes = sorted_dprimes[sorted_dprimes < limit]
  }
  else if(max_or_min == "min") {
    sorted_dprimes = sorted_dprimes[sorted_dprimes > limit]
  }
  pmf = rep(0, 2 * length(sorted_dprimes))
  dim(pmf) = c(length(sorted_dprimes), 2)
  pmf[,1] = sorted_dprimes
  # for each possible d',
  for(i in 1:length(sorted_dprimes)) {
    # the probability of observing that d' is the sum of the probabilities of
    # observing each combination of hit rate and false alarm rate that yields
    # that d'
    pmf[i,2] = sum(ps[dprimes == sorted_dprimes[i]])
  }
  
  return(pmf)
}

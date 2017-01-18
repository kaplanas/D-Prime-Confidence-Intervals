# These functions compute confidence intervals around the d'^ corresponding to
# the observed data according to the maximum likelihood estimation method.
# Perfect hit and false alarm rates are adjusted upward or downward by 0.5.  The
# functions compute 1 - a confidence intervals, where a is specified in an
# additional argument; the default is 0.05. The difference between discrete
# values of ph and pf is given in an additional argument; the default is 0.001.
# The functions return a numeric vector of the lower and upper bounds of the
# confidence interval.

MLE.int.count <- function(hits, misses, false_alarms, correct_rejections, a = 0.05, res = 0.001) {
  # get hits and 'signal' trials; store 'choose' computation
  h = hits
  S = hits + misses
  choose_S_h = choose(S, h)
  # get false alarms and 'noise' trials; store 'choose' computation
  f = false_alarms
  N = false_alarms + correct_rejections
  choose_N_f = choose(N, f)
  
  # store possible observed hit and false alarm rates
  phs = seq(0, 1, res)
  pfs = seq(0, 1, res)
  # create variables to hold calculated values based on possible hit and false
  # alarm rates
  pmfHs = rep(0, length(phs))
  pmfFs = rep(0, length(pfs))
  qnormHs = rep(0, length(phs))
  qnormFs = rep(0, length(pfs))
  Ls = rep(0, length(phs) * length(pfs))
  dim(Ls) = c(length(phs), length(pfs))
  # create variable to store d's in confidence set
  dprimes = rep(0, length(phs) * length(pfs))
  dim(dprimes) = c(length(phs), length(pfs))
  
  # for each possible hit rate,
  for(i in 1:length(phs)) {
    # calculate the probability of observing that hit rate
    pmfHs[i] = choose_S_h * (phs[i] ^ h) * ((1 - phs[i]) ^ (S - h))
    # adjust hit rate if necessary
    adjusted_ph = phs[i]
    if(adjusted_ph == 0) { adjusted_ph = .5 / S }
    else if(adjusted_ph == 1) { adjusted_ph = (S - .5) / S }
    # compute inverse-normal transform of hit rate
    qnormHs[i] = qnorm(adjusted_ph)
  }
  # analogous calculations for each possible false alarm rate
  for(j in 1:length(pfs)) {
    pmfFs[j] = choose_N_f * (pfs[j] ^ f) * ((1 - pfs[j]) ^ (N - f))
    adjusted_pf = pfs[j]
    if(adjusted_pf == 0) { adjusted_pf = .5 / N }
    else if(adjusted_pf == 1) { adjusted_pf = (N - .5) / N }
    qnormFs[j] = qnorm(adjusted_pf)
  }
  # for each possible combination of hit and false alarm rates,
  for(i in 1:length(phs)) {
    for(j in 1:length(pfs)) {
      # compute the likelihood of observed hit/false alarm rates if this is the
      # true hit/false alarm rate
      Ls[i,j] = pmfHs[i] * pmfFs[j]
      # compute the d' for that combination
      dprimes[i,j] = qnormHs[i] - qnormFs[j]
    }
  }
  
  # get maximum likelihood
  maxL = max(Ls)
  # get LRT statistic for each combination of hit and false alarm rates
  ps = Ls / maxL
  # interval is range of dprimes where we can't reject the null hypothesis that
  # the true d' equals the observed d'^
  int = range(dprimes[ps > a])
  return(int)
}

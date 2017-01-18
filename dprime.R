# These functions compute d' from a given number of hits, misses, false alarms,
# and correct rejections.  Perfect hit and false alarm rates are adjusted
# upward or downward by 0.5.

dprime.count <- function(hits, misses, false_alarms, correct_rejections) {
  # get hits and 'signal' trials
  h = hits
  S = hits + misses
  # adjust perfect hit rate, if necessary
  if(h == 0) {
    h = .5
  }
  if(h == S) {
    h = h - .5
  }
  # get false alarms and 'noise' trials
  f = false_alarms
  N = false_alarms + correct_rejections
  # adjust perfect false alarm rate, if necessary
  if(f == 0) {
    f = .5
  }
  if(f == N) {
    f = f - .5
  }
  # compute d'
  dprime = qnorm(h / S) - qnorm(f / N)
  return(dprime)
}

dprime.list <- function(responses) {
  # get hits and 'signal' trials
  h = length(responses[responses == "H"])
  S = length(responses[responses == "H" | responses == "M"])
  # adjust perfect hit rate, if necessary
  if(h == 0) {
    h = .5
  }
  if(h == S) {
    h = h - .5
  }
  # get false alarms and 'noise' trials
  f = length(responses[responses == "F"])
  N = length(responses[responses == "F" | responses == "R"])
  # adjust perfect false alarm rate, if necessary
  if(f == 0) {
    f = .5
  }
  if(f == N) {
    f = f - .5
  }
  # compute d'
  dprime = qnorm(h / S) - qnorm(f / N)
  return(dprime)
}

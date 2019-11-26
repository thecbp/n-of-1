# Function for creating backbone for data of an N-of-1 trial
create.backbone = function(order = c("A", "B"), s.freq = 1, 
                           p.length = 1, n.blocks = 1) {
  # Parameters
  # - order: the order that the treatments appear in
  # - s.freq: sampling frequency, how many times a patient will be measured in a 
  #           given treatment period
  # - p.length: period length
  # - n.blocks: how many treatment blocks there are
  
  # Working Assumptions
  # - Assumes that we will only be working with two treatments 
  
  # Build out a block skeleton detailing the treatment schedule for the subject
  single.block = tibble(
    treatment = rep(order, each = s.freq * p.length)
  )
  
  # Append this data structure to itself repeatedly
  all.blocks = NULL
  for (i in 1:n.blocks) {
    all.blocks = bind_rows(all.blocks, single.block)
  }
  
  # Correctly label day of trial and treatment blocks
  all.blocks = all.blocks %>% 
    mutate(
      block = rep(1:n.blocks, each = s.freq * p.length * length(order)),
      day = rep(1:(n.blocks * p.length * length(order)), each = s.freq)
    )
  
  return(all.blocks)
}

expdecay = function(start, target, tau, delta_t) {
  # Parameters
  # - start: the value we want to start at
  # - target: the value we want to end up at after all delta_t
  # - tau: scalar factor for how fast the decay happens (related to carryover, run-in)
  # - delta_t: vector of numbers to calculate the transition
  
  # Returns
  # - a vector of values going from a start value to a target value in a decay fashion
  
  # Notes:
  # if the decay profile doesn't seem correct, just tune it to your liking
  
  target + (start - target) * exp(-delta_t / tau)
}

collect.intervals = function(trt, trt.vec, type = "run") {
  # Parameters
  # - trt: the given treatment you want to find the treatment intervals for
  # - trt.vec: the vector of treatments
  # - match: do we want to see when the treatment is being used (TRUE) or not being used (FALSE)
  
  n = length(trt.vec)
  all.indices = c()
  matched = ifelse(type == "run", TRUE, FALSE)
  indices = integer()
  
  for (i in 1:n) {
    match = (trt.vec[i] == trt)
    # The character matches the current index in the treatment vector
    # Add the index
    if (match == matched) {
      indices = c(indices, i)
    }
    
    # Edge case for capturing treatments at the end
    if (match == matched & i == n) {
      indices = c(indices, i)
      all.indices = rbind(all.indices, c(min(indices), max(indices)))
    }
    
    # if the character didn't match and we have a bunch of indices
    # This signifies that there was a change in treatment
    # We should reset everything in preparation for seeing the treatment again
    if (match == !matched & length(indices) > 0) {
      all.indices = rbind(all.indices, c(min(indices), max(indices)))
      indices = integer()
    }
    
  }
  return(all.indices)
}

model.outcome = function(bb, 
                         trts,
                         mu.b = 100,
                         sd.b = 5,
                         sd.o = 1) {
  # Parameters
  # - bb: a treatment schedule backbone (create.backbone)
  # - trts: list of treatments to get the effect for
  # - mu.b: mean of the baseline value (what level should outcome take)
  # - sd.b: standard deviation of the baseline
  # - sd.o: standard deviation of observation noise
  # - trt.type: indicates if the treatments add/subtract from outcome ("-" for reducing effect)
  
  # Working Assumptions:
  # - just working with a continuous outcome for now
  # - baseline will just have regular (normal) noise
  # - assume that the carryover and run-in of any treatment always will be smaller than the 
  #     treatment period
  # - treatments will just be sequentially named as letters starting at A, B, ... 
  # - both treatments used will have the same effect direction (both reducing, etc)
  
  # Extra:
  # - expand to other outcome types?
  # - allow for interaction between treatments?
  
  # Create columns for the effects of all potential treatments
  effects = tibble(
    idx = 1:nrow(bb),
    A.eff = 0,
    B.eff = 0,
    C.eff = 0,
    D.eff = 0
  )
  
  eff.names = c("A.eff", "B.eff", "C.eff", "D.eff")
  
  for (i in 1:length(all.trts$testing)) {
    cur.trt = all.trts$testing[i]
    run.idx = collect.intervals(cur.trt, bb$treatment, type = "run")
    carry.idx = collect.intervals(cur.trt, bb$treatment, type = "carry")
    
    # Need to trim the carry for all the treatments that don't go first
    first.trt = bb$treatment[1]
    if (first.trt != cur.trt) {
      carry.idx = carry.idx[-1,]
    }
    
    # fill in the treatment effects using the run-in
    for (j in 1:nrow(run.idx)) {
      start = run.idx[j,1]
      end = run.idx[j,2]
      distance = end - start + 1
      effects[[eff.names[i]]][start:end] = expdecay(0, 
                                                    all.trts[[cur.trt]]$effect, 
                                                    all.trts[[cur.trt]]$run,
                                                    1:distance) +
        rnorm(distance, 0, all.trts[[cur.trt]]$sd) # process noise of the treatments
    }
    
    # carryover
    for (j in 1:nrow(carry.idx)) {
      start = carry.idx[j,1]
      end = carry.idx[j,2]
      distance = end - start + 1
      effects[[eff.names[i]]][start:end] = expdecay(all.trts[[cur.trt]]$effect, 
                                                    0, 
                                                    all.trts[[cur.trt]]$carry, 
                                                    1:distance) +
        rnorm(distance, 0, all.trts[[cur.trt]]$sd) # process noise of the treatments
    }
    
    
  }
  
  # Recombine the backbone with the effects
  final.bb = bind_cols(bb, effects)
  
  # Calculate the observed effect on the baseline 
  final.bb = final.bb %>% 
    mutate(
      baseline = rnorm(nrow(bb), mu.b, sd.b),
      obs = baseline + A.eff + B.eff + C.eff + D.eff + rnorm(nrow(bb), 0, sd.o)
    ) %>% 
    select(-idx)
  
  return(final.bb)
}

# Putting the two functions together to make simulations easy
simulate.trial = function(trts, mu.b = 100, sd.b = 5, sd.o = 1,
                          order = c("A", "B"), s.freq = 1, p.length = 1, n.blocks = 1) {
  
  # Create the backbone
  bb = create.backbone(order, s.freq, p.length, n.blocks)
  
  # Model the treatment effects
  full.data = model.outcome(bb, trts, mu.b = mu.b, sd.b = sd.b, sd.o = sd.o)
  
  return(full.data)
}

# USE CASE:

# trt.A = list(
#   name = "A",
#   effect = -30,
#   sd = 1,
#   run = 2,
#   carry = 3
# )
# 
# trt.B = list(
#   name = "B",
#   effect = -50,
#   sd = 1,
#   run = 4,
#   carry = 2
# )
#
# bb = create.backbone(order = c("A", "B", "B", "A"), 
#                      s.freq = 1, p.length = 10, n.blocks = 5)
# test = model.outcome(bb, all.trts, 100, 2) %>% 
#   pivot_longer(
#     ., 
#     baseline:obs,
#     names_to = "effect",
#     values_to = "value"
#   )
# 
# test %>% 
#   ggplot(aes(x = day, y = value, color = effect)) +
#   geom_line() +
#   labs(
#     title = "Visualization of simulated baseline, treatment effects, & observation",
#     x = "Day of treatment", 
#     y = "Outcome value"
#   ) + 
#   theme(legend.position = "bottom")

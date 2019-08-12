#################################################
# Arvid Edenheim - arved490                     #
# Linköpings University                         #
# Advanced Machine Learning - TDDE15            #
# Lab 2                                         #
# Last edited: 2019-08-12                       #
#################################################

### Setup ###

if(!require(HMM)) {
  install.packages("HMM")
}

library('HMM')

if(!require(entropy)) {
  install.packages("entropy")
}

library('entropy')

set.seed(1234)
tProb <- 0.5
nStates <- 10
states = symbols <- 1:10


############### 1 ###############

getNoiseInterval <- function(index, nStates, noise) {
  start <- index - noise
  end <- index + noise
  interval <- start:end %% nStates
  interval[interval==0] <- nStates
  return(interval)
}

# Initialize matrix with transition probabilities. P(X->X) = 0.5, P(X -> X+1) = 0.5
transProbs <- tProb * diag(nStates)
for(i in 1:nStates) {
  transProbs[i,(i %% nStates + 1)] <- tProb
}

# Initialize emission probabilities. P(X-2) = P(X-1) = P(X) = P(X+1) = P(X+2) = 0.2
eProb <- 0.2
emissionProbs <- matrix(0, nrow = nStates, ncol = nStates)
noise <- 2
for(i in 1:nStates) {
  emissionProbs[i,getNoiseInterval(i, nStates, noise)] <- eProb
}

HMM <- initHMM(States = states, Symbols = symbols, transProbs = transProbs, emissionProbs = emissionProbs)


############### 2 ###############

nSteps <- 100
sims <- simHMM(HMM, nSteps)


############### 3 ###############

# Exponential due to forward/backward returning log values
alpha <- exp(forward(HMM, sims$observation))
beta <- exp(backward(HMM, sims$observation))

# Normalizing at each time step to prevent underflow (?) and to, you know, normalize to obtain a prob. distr.
# Basically performs col[,i] / sum(col[,i]) to get the normalized distr.
filter <- prop.table(alpha, margin = 2)
smoother <- prop.table(alpha*beta, margin = 2)

# Evaluate most probable path using viterbi. Would be interesting to compare with Minimum Bayes Risk
viterbiPath = viterbi(HMM, sims$observation)
filterPath <- apply(filter, MARGIN = 2, which.max)
smootherPath <- apply(smoother, MARGIN = 2, which.max)

# Plot paths
plot(sims$states, type = 'l', main = 'Viterbi')
lines(viterbiPath, col='red')
lines(filterPath, col='green')
lines(smootherPath, col='blue')
legend(x='topright', fill = c('black', 'red', 'green', 'blue'), legend = c('States', 'Viterbi', 'Filter', 'Smoother'))


############### 4 ###############

accuracy <- function(path, states) {
  return(sum(diag(table(path, states)) / length(path)))
}

# 0.49
viterbiAccuracy = accuracy(viterbiPath, sims$states)

# 0.63
filterAccuracy =  accuracy(filterPath, sims$states)

# 0.75
smootherAccuracy =  accuracy(smootherPath, sims$states)


############### 5 ###############

sims <- simHMM(HMM, nSteps)

alpha <- exp(forward(HMM, sims$observation))
beta <- exp(backward(HMM, sims$observation))

filter <- prop.table(alpha, margin = 2)
smoother <- prop.table(alpha*beta, margin = 2)

viterbiPath = viterbi(HMM, sims$observation)
filterPath <- apply(filter, MARGIN = 2, which.max)
smootherPath <- apply(smoother, MARGIN = 2, which.max)

# 0.55
viterbiAccuracy = accuracy(viterbiPath, sims$states)

# 0.58
filterAccuracy =  accuracy(filterPath, sims$states)

# 0.73
smootherAccuracy =  accuracy(smootherPath, sims$states)


# Smoothing appears to perform better at evaluating the hidden states

############### 6 ###############

# With 300 time steps
nSteps <- 300
sims300 <- simHMM(HMM, nSteps)
alpha300 <- exp(forward(HMM, sims300$observation))
filter300 <- prop.table(alpha300, margin = 2)

entropy300 <- apply(filter300, MARGIN = 2, entropy.empirical)

# Evaluate rolling mean from step 1 to current step i
rollingMeanEntropy <- rep(0, nSteps)
for(i in 1:nSteps) {
  rollingMeanEntropy[i] <- mean(entropy300[1:i])
}

plot(rollingMeanEntropy, type='l', main="Rolling mean of entropy", ylab = 'Entropy', xlab='Step')

# Comment: Converges after approx 70-100 steps. More observations does not appear to provide further information


############### 7 ###############

hiddenStates101 <- colSums(filter[,100] * transProbs)

# P(X = 1) = 0.33, P(X = 2) = 0.50, P(X = 3) = 0.17, given filter[,100] = [0.65, 0.35, 0, ...]

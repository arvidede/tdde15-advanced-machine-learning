#################################################
# Arvid Edenheim - arved490                     #
# Linköpings University                         #
# Advanced Machine Learning - TDDE15            #
# Lab 1                                         #
# Last edited: 2019-xx-xx                       #
#################################################

##### Setup #####

if(!require(bnlearn)) {
  install.packages("bnlearn")
}

if(!require(gRain)) {
  install.packages("gRain")
}

library(bnlearn)
library(gRain)

###### Assignment 1 ######
# Show non-equivalent BN structures for multiple runs of the HC algorithm

data('asia')

par(mfrow=c(1,2))
# With log-likelihood as score function
set.seed(1)
hc_loglik <- hc(x = asia, score='loglik')
hc_loglik <- cpdag(hc_loglik)
plot(hc_loglik) # Generates complete network due to lack of parameter punishment

# With bayesian information criterion as score function
set.seed(1)
hc_bic <- hc(x = asia, score = 'bic')
hc_bic <- cpdag(hc_bic)
plot(hc_bic)

# Check equality
print(all.equal(hc_loglik, hc_bic))



###### Assignment 2 ######
# Learn a BN from 80% of the data and classify on the rest, compare with true BN

confusion_matrix <- function(BN, data, obs_var, pred_var) {
  predictions <- rep(0, nrow(data))
  for(i in 1:nrow(data)) {
    X <- NULL
    for(j in obs_var) {
      X[j] <- if(data[i,j] == 'yes') 'yes' else 'no'
    }
    find <- setEvidence(object=BN, nodes=obs_var, states=X) # Recommended over setFinding
    dist <- querygrain(object=find, nodes=pred_var)[[pred_var]]
    predictions[i] <- if(dist['yes'] >= 0.5) 'yes' else 'no'
  }
  return (table(data[,pred_var], predictions, dnn=c("TRUE", "PRED")))
}

set.seed(123)
samples <- sample(1:nrow(asia), size = floor(0.8*nrow(asia)), replace = F)
train <- asia[samples, ]
test  <- asia[-samples, ]

BN_train <- hc(train, restart=3, score = 'bic')
BN_true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

# Convert to grain and run Lauritzen-Spiegelhalter
BN_train <- compile(as.grain(bn.fit(BN_train, train)))
BN_true <- compile(as.grain(bn.fit(BN_true, train)))

conf_matrix <- confusion_matrix(BN_train, data=test, obs_var=c('A', 'D', 'X', 'E', 'B', 'L', 'T'), pred_var='S')
conf_matrix_true <- confusion_matrix(BN_true, data=test, obs_var=c('A', 'D', 'X', 'E', 'B', 'L', 'T'), pred_var='S')

# How can they be identical? Am I a pro?

###### Assignment 3 ######

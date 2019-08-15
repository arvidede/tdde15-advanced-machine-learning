#################################################
# Arvid Edenheim - arved490                     #
# Linköping University                          #
# Advanced Machine Learning - TDDE15            #
# Lab 4                                         #
# Last edited: 2019-08-15                       #
#################################################

####### Functions #######



####### 2.1 - Implementing GP Regression #######

####### 2.1.1 #######

####### 2.1.2 #######

####### 2.1.3 #######

####### 2.1.4 #######

####### 2.1.5 #######



####### 2.2 - GP Regression with kernlab #######

# Setup
read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv",
         header=TRUE,
         sep=";")

####### 2.2.1 #######

####### 2.2.2 #######

####### 2.2.3 #######

####### 2.2.4 #######

####### 2.2.5 #######



####### 2.3 - GP Classification with kernlab #######

# Setup
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv",
                 header=FALSE,
                 sep=",")

names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])

set.seed(111)
SelectTraining <- sample(1:dim(data)[1],
                         size = 1000,
                         replace = FALSE)

####### 2.3.1 #######

####### 2.3.2 #######

####### 2.3.3 #######

####### 2.3.4 #######

####### 2.3.5 #######

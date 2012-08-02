# reshaping bruce's data from wide to long
# I want each row to be a spot check 
# (spot check = counts from one time period)
# 8 Feb 2012

# this is the first script to run.
# it preps the data for "drift GLMM by time period.R"

data = read.csv(
"/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$N1 = data$N0 - data$t1
data$N2 = data$N1 - data$t2
data$N3 = data$N2 - data$t3
data$N4 = data$N3 - data$t4
data$N0t = data$N0


# first get rid of block, 
# because i don't much care about it
# or do we?
data = data[,-c(2)]


# using function reshape

ldata = reshape(data, varying=list(c('t1','t2','t3','t4'),
	c('N0t', 'N1', 'N2', 'N3'), c('N1', 'N2', 'N3', 'N4')),
	v.names=c('drift', 'initial', 'stay'), timevar='time', 
	direction='long', idvar='rep')
ldata$foodperinitial = ldata$food / ldata$initial
# each row is the count from one time one channel was checked
# initial is the number in the channel at the beginning of
# that time period.
# drift is number drifted
# stay is number stayed
# probability models for drift simulation and fitting
# CANONPY + DENSITY
# 16 Jan 2012

# run this before running 'CANOPY drift simulation by time period.R'
# or before 'MLE for drift data by time period.R'


# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

# logit transform:
logit = function(x){
	log( x / (1 - x))
}




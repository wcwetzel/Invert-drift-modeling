# probability models for drift simulation and fitting
# 14 Dec 2011

# run this before running 'drift simulation by time period.R'
# or before 'MLE for drift data by time period.R'


# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

# logit transform:
logit = function(x){
	log( x / (1 - x))
}


###############################
## models with no transform: ##

# constant rate
pm1 = function(m) { m } 

# linear decrease with food
pm2 = function(m, b) { m - b * data$food } 

# rational decrease with food
pm3 = function(m, c) { m / (c + data$food) } 

# linear increase with density
pm4 = function(m, a, N) { m + a * N } 

# linear decrease with food, increase with density
pm5 = function(m, a, b, N) { m + a * N - b * data$food }
 
# rational decrease with food, linear increase with density
pm6 = function(m, a, c, N) { (m + a * N) / (c + data$food) } 
###############################



###############################
## models with logisitic(y): ##

# constant rate
pm1 = function(m) {logistic( m )} 

# linear decrease with food
pm2 = function(m, b) {logistic( m - b * data$food )} 

# rational decrease with food
pm3 = function(m, c) {logistic( m / (c + data$food) )} 

# linear increase with density
pm4 = function(m, a, N) {logistic( m + a * N )} 

# linear decrease with food, increase with density
pm5 = function(m, a, b, N) {logistic( m + a * N - b * data$food )}
 
# rational decrease with food, linear increase with density
pm6 = function(m, a, c, N) {logistic( (m + a * N) / (c + data$food) )} 
###############################


###############################
## models with logisitic(m): ##

# constant rate
pm1 = function(m) {logistic( m )} 

# linear decrease with food
pm2 = function(m, b) {logistic( m ) - b * data$food } 

# rational decrease with food
pm3 = function(m, c) {logistic( m ) / (c + data$food) } 

# linear increase with density
pm4 = function(m, a, N) {logistic( m ) + a * N } 

# linear decrease with food, increase with density
pm5 = function(m, a, b, N) {logistic( m ) + a * N - b * data$food }
 
# rational decrease with food, linear increase with density
pm6 = function(m, a, c, N) { (logistic(m) + a * N) / (c + data$food) }

# rational decrease with food, linear increase with density
pm6n = function(a, c, N) { (a * N) / (c + data$food) }

###############################


###############################
## models with logisitic(m and b and a): ##

# constant rate
pm1 = function(m) {logistic( m )} 

# linear decrease with food
pm2 = function(m, b) {logistic( m ) - logistic(b) * data$food } 

# rational decrease with food
pm3 = function(m, c) {logistic( m ) / (c + data$food) } 

# linear increase with density
pm4 = function(m, a, N) {logistic( m ) + logistic(a) * N } 

# linear decrease with food, increase with density
pm5 = function(m, a, b, N) {logistic( m ) + logistic(a) * N - logistic(b) * data$food }
 
# rational decrease with food, linear increase with density
pm6 = function(m, a, c, N) { (logistic(m) + logistic(a) * N) / (c + data$food) } 
###############################





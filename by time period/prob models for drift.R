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
## models with logisitic(m): ##

# constant rate
pm1 = function(m) {logistic( m )} 

# linear decrease with food
pm2 = function(m, b) {logistic( m ) - b * data$food } 

# rational decrease with food
pm3 = function(m, c) {logistic( m ) / (c + data$food) } 

# rational decrease with food WITH k!!!
pm3k = function(m, k, c) {logistic( m ) + {logistic(k) / (c + data$food)} } 

# linear increase with density
pm4 = function(m, a, N) {logistic( m ) + a * N } 

# linear decrease with food, increase with density
pm5 = function(m, a, b, N) {logistic( m ) + a * N - b * data$food }
 
# rational decrease with food, linear increase with density
pm6 = function(m, a, c, N) { (logistic(m) + a * N) / (c + data$food) }

# rational decrease with food, linear increase with density WITH k!!!
pm6k = function(m, k, a, c, N) { logistic(m) + {(logistic(k) + a * N) / (c + data$food)} }

# rational decrease with food, linear increase with density WITH k!!! and no intercept
pm6kn = function(m, a, c, N) { logistic(m) + {(a * N) / (c + data$food)} }

# rational decrease with food, linear increase with density, NO intercept
pm6n = function(a, c, N) { (a * N) / (c + data$food) }

# rational decrease with food, nonlinear increase with density
pm7 = function(m, a, c, alpha, N) { (logistic(m) + {a * N}^alpha) / (c + data$food) }

# rational decrease with food, nonlinear increase with density
pm7k = function(m, k, a, c, alpha, N) { logistic(m) + (logistic(k) + {a * N}^alpha) / (c + data$food) }

# linear decrease with food, nonlinear increase with density
pm8 = function(m, a, b, alpha, N) {logistic( m ) + {a * N}^alpha - b * data$food }

# nonlinear increase with density
pm9 = function(m, a, alpha, N) {logistic( m ) + {a * N}^alpha } 

# constant rate WITH VELOCITY
pm1v = function(m, m1) {logistic( m + m1 * data$velocity )} 

# linear decrease with food WITH VELOCITY
pm2v = function(m, m1, b) {logistic( m + m1 * data$velocity ) - b * data$food } 

# rational decrease with food WITH VELOCITY
pm3v = function(m, m1, c) {logistic( m + m1 * data$velocity ) / (c + data$food) } 

# linear increase with density WITH VELOCITY
pm4v = function(m, m1, a, N) {logistic( m + m1 * data$velocity ) + a * N } 

# linear decrease with food, increase with density WITH VELOCITY
pm5v = function(m, m1, a, b, N) {logistic( m + m1 * data$velocity ) + a * N - b * data$food }
 
# rational decrease with food, linear increase with density WITH VELOCITY
pm6v = function(m, m1, a, c, N) { (logistic( m + m1 * data$velocity ) + a * N) / (c + data$food) }

# rational decrease with food, linear increase with density, WITH k!!!, WITH VELOCITY
pm6kv = function(m, m1, k, a, c, N) { logistic( m + m1 * data$velocity ) + {(logistic(k) + a * N) / (c + data$food)} }



pm3negexp = function(m, a, b) {logistic( m ) + a * exp(-b * data$food) } 


########################
#### with predation ####
# rational decrease with food, linear increase with density
pm6p = function(p, m, a, c, N) { logistic(p * data$preds) + (logistic(m) + a * N) / (c + data$food) }
pm6p(0, 0.1, 0.0001, 1, 113)


# just predation
pmp1 = function(p, m) {logistic(p * datacp$preds + m)}

# just canopy
pmp2 = function(c, m) {logistic(c * datacp$canopy + m)}

# predation and canopy: different intercepts, same slope
pmp3 = function(p, c, m) {logistic(p * datacp$preds + c * datacp$canopy + m)}

# predation and canopy: different intercepts, different slopes
pmp4 = function(p, c, m, cp) {
	logistic(m + p * datacp$preds + c * datacp$canopy + cp * datacp$preds * datacp$canopy)}

# rational decrease with canopy
pmp5 = function(p, m, c) { (p * datacp$preds + logistic(m)) / (c + datacp$revcanopy) } 

###############################




# ###############################
# ## models with no transform: ##

# # constant rate
# pm1 = function(m) { m } 

# # linear decrease with food
# pm2 = function(m, b) { m - b * data$food } 

# # rational decrease with food
# pm3 = function(m, c) { m / (c + data$food) } 

# # linear increase with density
# pm4 = function(m, a, N) { m + a * N } 

# # linear decrease with food, increase with density
# pm5 = function(m, a, b, N) { m + a * N - b * data$food }
 
# # rational decrease with food, linear increase with density
# pm6 = function(m, a, c, N) { (m + a * N) / (c + data$food) } 

# # rational decrease with food, nonlinear increase with density
# pm7 = function(m, a, c, alpha, N) { (m + a * {N^alpha}) / (c + data$food) }

# ###############################



# ###############################
# ## models with logisitic(y): ##

# # constant rate
# pm1 = function(m) {logistic( m )} 

# # linear decrease with food
# pm2 = function(m, b) {logistic( m - b * data$food )} 

# # rational decrease with food
# pm3 = function(m, c) {logistic( m / (c + data$food) )} 

# # linear increase with density
# pm4 = function(m, a, N) {logistic( m + a * N )} 

# # linear decrease with food, increase with density
# pm5 = function(m, a, b, N) {logistic( m + a * N - b * data$food )}
 
# # rational decrease with food, linear increase with density
# pm6 = function(m, a, c, N) {logistic( (m + a * N) / (c + data$food) )} 
# ###############################



# ###############################
# ## models with logisitic(m and b and a): ##

# # constant rate
# pm1 = function(m) {logistic( m )} 

# # linear decrease with food
# pm2 = function(m, b) {logistic( m ) - logistic(b) * data$food } 

# # rational decrease with food
# pm3 = function(m, c) {logistic( m ) / (c + data$food) } 

# # linear increase with density
# pm4 = function(m, a, N) {logistic( m ) + logistic(a) * N } 

# # linear decrease with food, increase with density
# pm5 = function(m, a, b, N) {logistic( m ) + logistic(a) * N - logistic(b) * data$food }
 
# # rational decrease with food, linear increase with density
# pm6 = function(m, a, c, N) { (logistic(m) + logistic(a) * N) / (c + data$food) } 
# ###############################





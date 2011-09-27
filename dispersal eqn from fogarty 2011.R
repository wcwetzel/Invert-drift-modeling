f = 10
N = 0:113

# my model 6 parms
b = 1
m = 0.1 # drift rate when N=0 and f--> infinity

# forgarty parms, it's just a logistic
a = 1


# dispersal probabilities/rates
D1 =  m * N / b + f # my model 6, linear relationship between per cap drift
# and density

D2 = 1 / (1 + exp(  f / (a * N) + log(((1 - 0.9) / 0.9))))


plot(D1 ~ N, type='l')
plot(D2 ~ N, type='l')


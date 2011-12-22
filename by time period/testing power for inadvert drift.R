library(bbmle)
data = data.frame(food=food, total=total, N0=213)

logistic = function(x){1/(1 + exp(-x))}
p3sk = function(m, k, c){logistic(m) + k/(c+data$food)}
p3s = function(k, c){k/(c+data$food)}

model3sk = mle2(total ~ dbinom(size=data$N0, prob=p3sk(m, k, c)),
	start=list(m = -10, k=0.01, c = 1), data=data)

model3s = mle2(total ~ dbinom(size=data$N0, prob=p3s(k, c)),
	start=list(k=0.01, c = 1), data=data)

AICtab(model3sk, model3s)

anova(model3sk, model3s)

model3sk

logistic(-3.2481687)

prom3sk = profile(model3sk)
ci3sk = confint(prom3sk)
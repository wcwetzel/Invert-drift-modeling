#### plotting Hept and Baet differences ####
## 13 April 2012 ##

## first run "drift data reshape.R"
## second run "drift data reshape by species.R"

par(mfrow=c(2,2))
plot(drift ~ time, data=ldata)
plot(drifth ~ time, data=ldataspp)
plot(driftb ~ time, data=ldataspp)

par(mfrow=c(2,2))
plot(drift ~ food, data=ldata[ldata$experiment=='excess_food',])
plot(drifth ~ food, data=ldataspp[ldataspp$experiment=='excess_food',])
plot(driftb ~ food, data=ldataspp[ldataspp$experiment=='excess_food',])

par(mfrow=c(2,2))
plot(drift ~ N0, data=ldata[ldata$experiment=='density',])
plot(drifth ~ N0, data=ldataspp[ldataspp$experiment=='density',])
plot(driftb ~ N0, data=ldataspp[ldataspp$experiment=='density',])

par(mfrow=c(2,2))
plot(drift ~ preds, data=ldata[ldata$experiment=='pred_canopy',])
plot(drifth ~ preds, data=ldataspp[ldataspp$experiment=='pred_canopy',])
plot(driftb ~ preds, data=ldataspp[ldataspp$experiment=='pred_canopy',])

par(mfrow=c(2,2))
plot(drift ~ canopy, data=ldata[ldata$experiment=='pred_canopy',])
plot(drifth ~ canopy, data=ldataspp[ldataspp$experiment=='pred_canopy',])
plot(driftb ~ canopy, data=ldataspp[ldataspp$experiment=='pred_canopy',])


# h and b look very similar



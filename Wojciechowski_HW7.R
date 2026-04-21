#Homework 7

##Tutorial 

# read in greenhouse gas data from reservoirs
ghg <- read.csv("/cloud/project/activity07/Deemer_GHG_Data.csv")

install.packages(c("dplyr", "ggplot2","olsrr", "PerformanceAnalytics"))

library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)

ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)

# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)

# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)


# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe

summary(mod.full)

res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)

plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 

# check full model
full.step$model

# plot AIC over time
plot(full.step )

# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")
# look at prediction with 95% confidence interval of the mean

predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")


#Question 3

#log transformation depth to reduce skew 

ghg$log.depth <- log(ghg$mean.depth + 1)

mod.improved <- lm(log.ch4 ~ airTemp+
                 log.age+
                  log.depth+
                 log.DIP+
                 log.precip+ 
                  BorealV +
                  TropicalV , data=ghg)

summary(mod.improved)

table_out <- summary(mod.improved)$coefficients

write.csv(table_out, "/cloud/project/mod.csv", row.names =TRUE)

##Assumptions Check for improved model

res.imp <- rstandard(mod.improved)
fit.imp <- fitted.values(mod.improved)


qqnorm(res.imp, pch = 19, col = "grey50", main = "Q-Q Plot: Improved Model")
qqline(res.imp)
shapiro.test(res.imp)

plot(fit.imp, res.imp, pch = 19, col = "grey50",
     xlab = "Fitted Values", ylab = "Standardized Residuals",
     main = "Residuals vs Fitted: Improved Model")
abline(h = 0)

ols_vif_tol(mod.improved)

imp.step <- ols_step_forward_aic(mod.improved)
imp.step
plot(imp.step)
summary(imp.step$model)

predict.lm(imp.step$model,
           data.frame(airTemp    = 25,
                      log.age    = log(2),
                      log.depth  = log(16),
                      log.DIP    = log(51),
                      log.precip = log(1500),
                      BorealV    = 0,
                      TropicalV  = 1),
           interval = "prediction")

##start of chapter 7 tutorial 

ETdat <- read.csv("/cloud/project/activity07/ETdata.csv")

unique(ETdat$crop)

install.packages("lubridate")
library(lubridate)
library(ggplot2)
library(forecast)
library(dplyr)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

almondTrend <- almond_dec$trend
almondSeason <- almond_dec$seasonal

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")



##start of homework 

## Question 1

##transformation for co2
ghg$co2.transformation <- 1/ (ghg$co2 + 1000)

## model implementation 
model.co2 <- lm(co2.transformation ~ log.age + log.depth + BorealV +
                  log.DIP + log.precip, ghg)
summary(model.co2)

##check assumptions 

res.co2 <- rstandard(model.co2)
fit.co2 <- fitted.values(model.co2)
qqnorm(res.co2, pch = 19, col = "grey50", main = "Q-Q Plot: co2")
qqline(res.co2)
shapiro.test(res.co2)

plot(fit.co2, res.co2, pch = 19, col = "grey50",
     xlab = "Fitted Values", ylab = "Standardized Residuals",
     main = "Residuals vs Fitted: co2")
abline(h = 0)

ols_vif_tol(model.co2)

co2.step <- ols_step_forward_aic(model.co2)
co2.step
plot(co2.step)
summary(co2.step$model)

## print out table for summary 
table_out <- summary(model.co2)$coefficients

write.csv(table_out, "/cloud/project/model.csv", row.names =TRUE)

nrow(model.co2$model)

## Question 3.

##decompose almonds

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

## decompose pistachios


pistachios <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(pistachios, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# ET time series
pistachios_ts <- ts(pistachios$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose ET time series
pistachios_dec <- decompose(pistachios_ts)
# plot decomposition
plot(pistachios_dec)

##decompose fallow/idle fields


fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(fallow, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# ET time series
fallow_ts <- ts(fallow$ET.in, # data
                    start = c(2016,1), #start year 2016, month 1
                    #first number is unit of time and second is observations within a unit
                    frequency= 12) # frequency of observations in a unit

# decompose ET time series
fallow_dec <- decompose(fallow_ts)
# plot decomposition
plot(fallow_dec)

##decompose corn time series 

corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(corn, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

#  ET time series
corn_ts <- ts(corn$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose ET time series
corn_dec <- decompose(corn_ts)
# plot decomposition
plot(corn_dec)

## decompose table grapes time series


grapes <- ETdat %>% # ET data
  filter(crop == "Grapes (Table/Raisin)") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(grapes, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

#  ET time series
grapes_ts <- ts(grapes$ET.in, # data
              start = c(2016,1), #start year 2016, month 1
              #first number is unit of time and second is observations within a unit
              frequency= 12) # frequency of observations in a unit

# decompose ET time series
grapes_dec <- decompose(grapes_ts)
# plot decomposition
plot(grapes_dec)


## Question 4 

##pistachio ar model 

pistachiosTrend <- pistachios_dec$trend
pistachiosSeason <- pistachios_dec$seasonal

acf(na.omit(pistachios_ts),
    lag.max = 24)

pacf.plot <- pacf(na.omit(pistachios_ts))

pistachios_y <- na.omit(pistachios_ts)

model1 <- arima(pistachios_y,
                order = c(1,0,0))
model1

model4 <- arima(pistachios_y,
                order = c(4,0,0))
model4

AR_fit1 <- pistachios_y - residuals(model1)
AR_fit4 <- pistachios_y - residuals(model4)

plot(pistachios_y)
points(AR_fit1, type = "l", col = "tomato3",        lty = 2, lwd = 2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd = 2)
legend("topleft", c("data","AR1","AR4"),
       lty = c(1,2,2), lwd = c(1,2,2),
       col = c("black","tomato3","darkgoldenrod4"),
       bty = "n")

newPistachio <- forecast(model4)
newPistachio

newPistachioF <- data.frame(newPistachio)

years <- c(rep(2021,4), rep(2022,12), rep(2023,8))
month <- c(seq(9,12), seq(1,12), seq(1,8))
newPistachioF$dateF <- ymd(paste(years,"/",month,"/",1))

ggplot() +
  geom_line(data = pistachios, aes(x = ymd(date), y = ET.in)) +
  xlim(ymd(pistachios$date[1]), newPistachioF$dateF[24]) +
  geom_line(data = newPistachioF, aes(x = dateF, y = Point.Forecast),
            col = "red") +
  geom_ribbon(data = newPistachioF,
              aes(x = dateF, ymin = Lo.95, ymax = Hi.95),
              fill = rgb(0.5,0.5,0.5,0.5)) +
  theme_classic() +
  labs(x = "year", y = "Evapotranspiration (in)", main = "Pistachois")

## fallow ar model 

fallowTrend <- fallow_dec$trend
fallowSeason <- fallow_dec$seasonal

acf(na.omit(fallow_ts),
    lag.max = 24)

pacf.plot <- pacf(na.omit(fallow_ts))

fallow_y <- na.omit(fallow_ts)

model1 <- arima(fallow_y,
                order = c(1,0,0))
model1

model4 <- arima(fallow_y,
                order = c(4,0,0))
model4

AR_fit1 <- fallow_y - residuals(model1)
AR_fit4 <- fallow_y - residuals(model4)

plot(fallow_y)
points(AR_fit1, type = "l", col = "tomato3",        lty = 2, lwd = 2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd = 2)
legend("topleft", c("data","AR1","AR4"),
       lty = c(1,2,2), lwd = c(1,2,2),
       col = c("black","tomato3","darkgoldenrod4"),
       bty = "n")

newFallow <- forecast(model4)
newFallow

newFallowF <- data.frame(newFallow)

years <- c(rep(2021,4), rep(2022,12), rep(2023,8))
month <- c(seq(9,12), seq(1,12), seq(1,8))
newFallowF$dateF <- ymd(paste(years,"/",month,"/",1))

ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in)) +
  xlim(ymd(fallow$date[1]), newFallowF$dateF[24]) +
  geom_line(data = newFallowF, aes(x = dateF, y = Point.Forecast),
            col = "red") +
  geom_ribbon(data = newFallowF,
              aes(x = dateF, ymin = Lo.95, ymax = Hi.95),
              fill = rgb(0.5,0.5,0.5,0.5)) +
  theme_classic() +
  labs(x = "year", y = "Evapotranspiration (in)")

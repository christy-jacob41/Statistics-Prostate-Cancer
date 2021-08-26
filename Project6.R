# Christy Jacob
# Mini Project 6

# read in file along with the variables
data <- read.csv("prostate_cancer.csv")
psa <- data[,2]
cancervol <- data[,3]
weight <- data[,4]
age <- data[,5]
benpros <- data[,6]
vesinv <- data[,7]
capspen <- data[,8]
gleason <- data[,9]

# plotting histogram of psa
hist(psa, main="psa")

# plotting boxplot of psa
boxplot(psa)

# plotting boxplot of psa that is natural log transformed
boxplot(log(psa))

# plotting qq plot of psa
qqnorm(psa)
qqline(psa)

# histograms of the variables that aren't the response variable
hist(cancervol, main="cancervol")
hist(weight, main="weight")
hist(age, main="age")
hist(benpros, main="benpros")
hist(vesinv, main="vesinv")
hist(capspen, main="capspen")
hist(gleason, main="gleason")

# plotting scatterplots of the variables without natural log transformation
pairs(~psa + cancervol + weight + age + benpros + capspen + gleason, data = data)

# getting correlation of all the variables
correl <- cor(data[,2:9])
correl

# plotting boxplot of vesinv with psa
boxplot(psa~vesinv)

# plotting scatterplots of the variables with natural log transformation
pairs(~log(psa) + cancervol + capspen + gleason, data = data)

# plotting the natural log transformation of psa with gleason, capspen, and cancervol
plot(cancervol, log(psa))
abline(lm(log(psa)~cancervol))
plot(capspen, log(psa))
abline(lm(log(psa)~capspen))
plot(gleason, log(psa))
abline(lm(log(psa)~gleason))

# replacing psa with log(psa)
data[,2] <- log(psa)

# getting correlation of the variables with the natural log transformation
correl.log <- cor(data[,2:9],log(psa))
correl.log

# plotting boxplot of vesinv with natural log transformation of psa
boxplot(log(psa)~vesinv)

# building  model using 4 predictors
model1 <- lm(log(psa)~cancervol+capspen+gleason+vesinv)
model1
summary(model1)

# building model using 3 predictors wihout capspen
model2 <- lm(log(psa)~cancervol+gleason+vesinv)
model2
summary(model2)

# partial F test to see if capspen is significant
anova(model1, model2)

# building model using 2 predictors, without gleason
model3 <- lm(log(psa)~cancervol+vesinv)
model3
summary(model3)

# building model using 2 predictors, without cancervol
model4 <- lm(log(psa)~gleason+vesinv)
model4
summary(model3)

# building model using 2 predictors, without vesinv
model5 <- lm(log(psa)~cancervol+gleason)
model5
summary(model3)

# partial F test to see if gleason is significant
anova(model2, model3)

# partial F test to see if cancervol is significant
anova(model2, model4)

# partial F test to see if vesinv is significant
anova(model2, model5)

nat.log.psa <- log(psa)

# AIC based models
# forward selection AIC
forward <- step(lm(nat.log.psa ~ 1, data = data), scope = list(upper = ~cancervol+capspen+gleason+vesinv),direction = "forward")
# backward elimination AIC
backward <- step(lm(nat.log.psa~cancervol+capspen+gleason+vesinv, data = data), scope = list(lower = ~1), direction = "backward")
# hybrid AIC of forward and backward
hybrid <-  step(lm(nat.log.psa ~ 1, data = data), scope = list(lower = ~1, upper = ~cancervol+capspen+gleason+vesinv),direction = "both")

# summary analysis of preliminary model
summary(model2)
# residual plot of preliminary model
plot(fitted(model2), resid(model2))
abline(h=0)
# absolute residual plot of preliminary model
plot(fitted(model2), abs(resid(model2)))
# normal qq plot of preliminary model
qqnorm(resid(model2))
qqline(resid(model2))

# calculating psa level for patients whose quantitative predictors are at the sample means of the variable & qualitative predictors are at the most frequency category 
mean_cancervol <- mean(cancervol)
mean_gleason <- mean(gleason)
psa_predicted <- -0.72120 + (0.05981 * mean_cancervol) + (0.38491 * mean_gleason) + (0.62117 * 0)
psa_predicted
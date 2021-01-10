# VIF removal, outlier removal

train <- read.csv(file="protein-train.csv")

library(leaps)
library(MASS)
attach(train)

#Declaring these as null

train$scArgN_bbC_medshort = NULL
train$scArgN_bbO_short = NULL

m1 <- lm(accuracy ~ ., data = train)

vifval <- max(vif(m1))
removed <- vector()

while(vifval > 10)
{
  removal <- names(which.max(vif(m1)))
  m1 <- update(m1, paste0(".~. - ", removal))
  vifval <- max(vif(m1))
  removed <- c(removed, removal)
}

train = train[,!(names(train) %in% removed)]


m1 <- lm(accuracy ~ ., data = train)

outliers = which(cooks.distance(m1) > 4/nrow(train))
train <- train[-outliers,]


#Seed remains the same for all

N <- nrow(train)
set.seed(123423)
trainInd <- sample(1:N, round(N*0.8), replace=F)
trainSet <- train[trainInd,]
validSet <- train[-trainInd,]

#Declaring the models

full <- lm(accuracy ~ ., data = trainSet)
empty <- lm(accuracy ~ 1, data = trainSet)

#Model 1 - BIC with standard penalty

m_k_bic <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
                   direction = "both", trace = 0, k = log(nrow(trainSet)))
pred1 <- predict(m_k_bic, newdata = validSet)
RMSE1 <- sqrt(mean((validSet$accuracy - pred1)^2))
RMSE1_t <- (mean(m_k_bic$residuals^2)) #Rmse on train

#Model 2 - BIC with double penalty

m_k_bic_2 <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
                     direction = "both", trace = 0, k = 2 * log(nrow(trainSet)))
pred2 <- predict(m_k_bic_2, newdata = validSet)
RMSE2 <- sqrt(mean((validSet$accuracy - pred2)^2))
RMSE2_t <- (mean(m_k_bic_2$residuals^2)) #Rmse on train

#Model 3 - BIC with quad penalty

m_k_bic_3 <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
                     direction = "both", trace = 0, k = 4 * log(nrow(trainSet)))
pred3 <- predict(m_k_bic_3, newdata = validSet)
RMSE3 <- sqrt(mean((validSet$accuracy - pred3)^2))
RMSE3_t <- (mean(m_k_bic_3$residuals^2)) #Rmse on train



#Model 5 - ICM with BIC penalty

pen <- log(nrow(trainSet))
varlist = c()
varnames = names(trainSet)
n = nrow(trainSet)
varorder <- sample(1:ncol(trainSet))
minCrit = Inf
noChange = F
while (!noChange) {
  noChange = T
  for (i in varorder) { 
    if (i == 1)
      next
    
    if (i %in% varlist & length(varlist) > 1) {
      index = c(1, varlist[varlist != i]) 
      trainVars = trainSet[, index]
      
      fit = lm(accuracy ~ ., data = trainVars)
      
      if (AIC(fit, k = pen) < minCrit) {
        minCrit = AIC(fit, k = pen)
        varlist = varlist[varlist != i]
        print(paste0("Criterion: ", round(minCrit, 1), ", variables: ", paste0(varnames[varlist], collapse = " ")))
        best.model = fit
        noChange = F
      }
      
    } else if (!i %in% varlist) {
      index = c(1, varlist, i) 
      trainVars = trainSet[, index]
      
      fit = lm(accuracy ~ ., data = trainVars)
      
      if (AIC(fit, k = pen) < minCrit) {
        minCrit = AIC(fit, k = pen)
        varlist = c(varlist, i)
        print(paste0("Criterion: ", round(minCrit, 1), ", variables: ", paste0(varnames[varlist], collapse = " ")))
        best.model = fit
        noChange = F
      }      
    }
  }
}

predICM <- predict(best.model, newdata = validSet)
RMSE5 = sqrt(mean((validSet$accuracy - predICM)^2)) # RMSE on validation
RMSE5_t = sqrt(mean(best.model$residuals^2)) # RMSE on train

#Model 6: full model

m6 <- lm(accuracy ~ ., data = trainSet)

pred6 <- predict(m6, newdata = validSet)
RMSE6 <- sqrt(mean((validSet$accuracy - pred6)^2))
RMSE6_t <- (mean(m6$residuals^2)) #Rmse on train

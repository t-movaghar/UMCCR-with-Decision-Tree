#install.packages("e1071")
library(e1071)

# Breast Cancer = Risk 1 (primary risk)
# Other Cancers = Risk 2 
# Other Causes = Risk 3
setwd("C:/Users/taram/Desktop/Dr Pal/")

brst.data = read.table("Brstdata.txt",header=T)
head(brst.data)

for(i in 1:length(brst.data$T)){
  if(brst.data$Cause[i]=="Other_Cancers"){
    brst.data$D[i] = 2
  }else if(brst.data$Cause[i]=="Other_Causes"){
    brst.data$D[i] = 3
  }
}

#head(brst.data)
#stat.desc(brst.data)

# finding the initial value of the cure rate

library(survival)
library(survminer)
library(rpart)


J = rep(NA,length(brst.data$T))
for(j in 1:length(brst.data$T)){
  if(brst.data$Cause[j]=="Breast_Cancer"){
    J[j] = 1
  }else{
    J[j] = 0
  }
}

head(brst.data)
J

data.init = data.frame(Y=brst.data$T/12,J)
KM.init = survfit(Surv(Y,J)~1,data=data.init)
ggsurvplot(KM.init,data=data.init,xlab="Time (years)")

# survival probability at the largest uncensored time = 0.289, which is an initial guess for cure rate

# finding the initial values for the Weibull parameters for each risk type

library(BB)

# Risk 1:

data1 = brst.data[brst.data$D==1,2]
data1.yr = data1/12
mean(data1.yr) # 2.41
var(data1.yr) # 3.56

# Risk 2:

data2 = brst.data[brst.data$D==2,2]
data2.yr = data2/12
mean(data2.yr) # 2.11
var(data2.yr) # 4.29


# Risk 3:

data3 = brst.data[brst.data$D==3,2]
data3.yr = data3/12
mean(data3.yr) # 2.66
var(data3.yr) # 4.51


f1 = function(v=c(alp,lam)){
  f = rep(NA,2)
  f[1] =  ((1/(v[2]^(1/v[1])))*gamma(1+(1/v[1]))) - 2.66
  f[2] =  ((1/(v[2]^(1/v[1])))^2)*((gamma(1+(2/v[1])))-((gamma(1+(1/v[1])))^2)) - 4.51
  f
}

BBsolve(par=c(0.5,0.5),fn=f1)

# alp1.init = 1.29; lam1.init = 0.29
# alp2.init = 1.02; lam2.init = 0.46
# alp3.init = 1.26; lam3.init = 0.26

# calculation of observed log-likelihod function

log.lik.fn = function(param=c(cure,alp1,lam1,alp2,lam2,alp3,lam3),y0,y1,y2,y3){
  
  h1.0 = param[2]*param[3]*(y0^(param[2]-1))
  h1.1 = param[2]*param[3]*(y1^(param[2]-1)) 
  h1.2 = param[2]*param[3]*(y2^(param[2]-1))
  h1.3 = param[2]*param[3]*(y3^(param[2]-1))
  
  h2.0 = param[4]*param[5]*(y0^(param[4]-1))
  h2.1 = param[4]*param[5]*(y1^(param[4]-1)) 
  h2.2 = param[4]*param[5]*(y2^(param[4]-1))
  h2.3 = param[4]*param[5]*(y3^(param[4]-1))
  
  h3.0 = param[6]*param[7]*(y0^(param[6]-1))
  h3.1 = param[6]*param[7]*(y1^(param[6]-1)) 
  h3.2 = param[6]*param[7]*(y2^(param[6]-1))
  h3.3 = param[6]*param[7]*(y3^(param[6]-1))
  
  H1.0 = param[3]*(y0^param[2])
  H1.1 = param[3]*(y1^param[2])
  H1.2 = param[3]*(y2^param[2])
  H1.3 = param[3]*(y3^param[2])
  
  H2.0 = param[5]*(y0^param[4])
  H2.1 = param[5]*(y1^param[4])
  H2.2 = param[5]*(y2^param[4])
  H2.3 = param[5]*(y3^param[4])
  
  H3.0 = param[7]*(y0^param[6])
  H3.1 = param[7]*(y1^param[6])
  H3.2 = param[7]*(y2^param[6])
  H3.3 = param[7]*(y3^param[6])
  
  l.obs = sum(log(1-param[1])+log(h1.1)-H1.1-H2.1-H3.1)+sum(log((param[1]*exp(-H2.0-H3.0))+((1-param[1])*exp(-H1.0-H2.0-H3.0))))+
    sum(log(h2.2))+sum(log((param[1]*exp(-H2.2-H3.2))+((1-param[1])*exp(-H1.2-H2.2-H3.2))))+
    sum(log(h3.3))+sum(log((param[1]*exp(-H2.3-H3.3))+((1-param[1])*exp(-H1.3-H2.3-H3.3))))
  return(l.obs)
}


library(rpart)  
library(numDeriv)  

EM.MCCR3 <- function(data, train_ratio = 0.7, tol = 1e-5, c = 0.5, maxit = 100,
                     alpha1 = 1, lambda1 = 1, alpha2 = 1, lambda2 = 1, alpha3 = 1, lambda3 = 1) {
  
  library(e1071)
  
  set.seed(123)
  train_indices <- sample(1:nrow(data), size = floor(train_ratio * nrow(data)))
  training_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Split training into train/validation
  val_indices <- sample(1:nrow(training_data), size = floor(0.2 * nrow(training_data)))
  valid_split <- training_data[val_indices, ]
  training_data <- training_data[-val_indices, ]
  
  training_data$d <- factor(training_data$d, levels = c("0", "1", "2", "3"))
  test_data$d <- factor(test_data$d, levels = c("0", "1", "2", "3"))
  
  # Fit best SVM ONCE using training/validation data
  cost_vals <- c(0.1, 1, 10)
  gamma_vals <- c(0.01, 0.1, 1)
  
  best_accuracy <- 0
  best_model <- NULL
  best_params <- list()
  
  for (cost in cost_vals) {
    for (gamma in gamma_vals) {
      model <- svm(d ~ ., data = training_data, kernel = "radial",
                   cost = cost, gamma = gamma, probability = TRUE)
      
      pred <- predict(model, newdata = valid_split)
      acc <- mean(pred == valid_split$d)
      
      if (acc > best_accuracy) {
        best_accuracy <- acc
        best_model <- model
        best_params <- list(cost = cost, gamma = gamma)
      }
    }
  }
  
  cat("Best SVM Params:\n")
  print(best_params)
  cat("Validation Accuracy:", best_accuracy, "\n")
  
  pred_prob <- attr(predict(best_model, newdata = test_data, probability = TRUE), "probabilities")
  cured_prob <- if ("0" %in% colnames(pred_prob)) pred_prob[, "0"] else rep(0.5, nrow(test_data))
  
  # Split by event type
  data0 <- training_data[training_data$d == 0, ]
  data1 <- training_data[training_data$d == 1, ]
  data2 <- training_data[training_data$d == 2, ]
  data3 <- training_data[training_data$d == 3, ]
  
  y0 <- data0$y
  y1 <- data1$y
  y2 <- data2$y
  y3 <- data3$y
  
  p.old <- matrix(c(c, alpha1, lambda1, alpha2, lambda2, alpha3, lambda3), ncol = 1)
  continue <- TRUE
  iter <- 1
  
  while (continue) {
    # Current parameters
    a1 <- p.old[2, 1]; l1 <- p.old[3, 1]
    a2 <- p.old[4, 1]; l2 <- p.old[5, 1]
    a3 <- p.old[6, 1]; l3 <- p.old[7, 1]
    
    # Cumulative hazards
    H1.0 <- l1 * y0^a1; H1.2 <- l1 * y2^a1; H1.3 <- l1 * y3^a1
    H2.0 <- l2 * y0^a2; H2.2 <- l2 * y2^a2; H2.3 <- l2 * y3^a2
    H3.0 <- l3 * y0^a3; H3.2 <- l3 * y2^a3; H3.3 <- l3 * y3^a3
    
    # Update weights
    numer <- (1 - cured_prob) * exp(-H1.0 - H2.0 - H3.0)
    denom <- cured_prob * exp(-H2.0 - H3.0) + numer
    w0 <- numer / denom
    
    w2 <- ((1 - p.old[1, 1]) * exp(-H1.2 - H2.2 - H3.2)) /
      ((p.old[1, 1] * exp(-H2.2 - H3.2)) + ((1 - p.old[1, 1]) * exp(-H1.2 - H2.2 - H3.2)))
    w3 <- ((1 - p.old[1, 1]) * exp(-H1.3 - H2.3 - H3.3)) /
      ((p.old[1, 1] * exp(-H2.3 - H3.3)) + ((1 - p.old[1, 1]) * exp(-H1.3 - H2.3 - H3.3)))
    
    # Q1
    Q1 <- function(par) {
      a <- par[1]; l <- par[2]
      res <- sum(log(a * l * y1^(a - 1))) - sum(l * y1^a) - sum(w0 * l * y0^a)
      return(-res)
    }
    risk1.new <- optim(par = c(a1, l1), fn = Q1, method = "Nelder-Mead")$par
    
    # Q2
    Q2 <- function(par) {
      a <- par[1]; l <- par[2]
      res <- sum(log(a * l * y2^(a - 1))) - sum(l * y2^a) - sum(w0 * l * y0^a)
      return(-res)
    }
    risk2.new <- optim(par = c(a2, l2), fn = Q2, method = "Nelder-Mead")$par
    
    # Q3
    Q3 <- function(par) {
      a <- par[1]; l <- par[2]
      res <- sum(log(a * l * y3^(a - 1))) - sum(l * y3^a) - sum(w0 * l * y0^a)
      return(-res)
    }
    risk3.new <- optim(par = c(a3, l3), fn = Q3, method = "Nelder-Mead")$par
    
    # Update
    p.new <- matrix(c(mean(cured_prob), risk1.new, risk2.new, risk3.new), ncol = 1)
    
    continue <- (max(abs(p.new - p.old)) > tol) && (iter < maxit)
    p.old <- p.new
    iter <- iter + 1
  }
  
  # Final test predictions
  pred <- predict(best_model, newdata = test_data)
  cm <- table(Predicted = pred, Actual = test_data$d)
  accuracy <- sum(diag(cm)) / sum(cm)
  
  precision <- diag(cm) / rowSums(cm)
  recall <- diag(cm) / colSums(cm)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  cat("\nFinal SVM Model Metrics:\n")
  print(cm)
  cat("Accuracy:", accuracy, "\n")
  for (i in 1:4) {
    cat(sprintf("Class %d Precision: %.4f\n", i - 1, precision[i]))
    cat(sprintf("Class %d Recall: %.4f\n", i - 1, recall[i]))
    cat(sprintf("Class %d F1-Score: %.4f\n", i - 1, f1_score[i]))
  }
  
  list(parameters = p.new, iterations = iter, cured_fraction = cured_prob)
  print(iter)
  print(p.new)
}




# calling the functions

library(numDeriv)

data = data.frame(y=brst.data$T/12,d=brst.data$D)
#data = data.frame(brst.data)

brst.data.new <- brst.data[, !(names(brst.data) %in% c("D", "T"))]  # Drop columns D and T

# Add new columns y and d
brst.data.new$y <- data$y
brst.data.new$d <- data$d

brst.data.new <- brst.data.new[, !(names(brst.data.new) %in% c("No", "Cause", "Race", "Type"))]
#data

brst.data.new$race_white <- ifelse(brst.data.new$race_black == 0 & brst.data.new$race_asian == 0, 1, 0)
brst.data.new$luminal_a <- ifelse(brst.data.new$luminal_b == 0 & brst.data.new$her2_enrch == 0 & brst.data.new$unkn == 0 & brst.data.new$tripl_neg == 0, 1, 0)


head(brst.data.new)


c.init = 0.289
alp1.init = 1.29
lam1.init = 0.29
alp2.init = 1.02
lam2.init = 0.46
alp3.init = 1.26
lam3.init = 0.26

EM.MCCR3(brst.data.new,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
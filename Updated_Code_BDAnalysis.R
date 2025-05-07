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

EM.MCCR3 <- function(data, train_ratio = 0.7, tol, c, maxit, alpha1, lambda1, alpha2, lambda2, alpha3, lambda3) {
  
  set.seed(123)  # For reproducibility
  train_indices <- sample(1:nrow(data), size = floor(train_ratio * nrow(data)))
  training_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Split training data by event type
  data0 = training_data[training_data$d == 0, ]  # Alive
  data1 = training_data[training_data$d == 1, ]
  data2 = training_data[training_data$d == 2, ]
  data3 = training_data[training_data$d == 3, ]
  
  y0 = data0$y
  y1 = data1$y
  y2 = data2$y
  y3 = data3$y
  
  p.new = matrix(0, ncol = 1, nrow = 7)
  p.old = matrix(0, ncol = 1, nrow = 7)
  
  p.old[1, 1] = c 
  p.old[2, 1] = alpha1
  p.old[3, 1] = lambda1
  p.old[4, 1] = alpha2
  p.old[5, 1] = lambda2
  p.old[6, 1] = alpha3
  p.old[7, 1] = lambda3
  
  continue = TRUE
  iter = 1
  
  while (continue) {
    
    # E-step: Creating hazard functions
    h1.0 <- p.old[2,1] * p.old[3,1] * (y0^(p.old[2,1] - 1))
    h1.1 <- p.old[2,1] * p.old[3,1] * (y1^(p.old[2,1] - 1))
    h1.2 <- p.old[2,1] * p.old[3,1] * (y2^(p.old[2,1] - 1))
    h1.3 <- p.old[2,1] * p.old[3,1] * (y3^(p.old[2,1] - 1))
    
    h2.0 <- p.old[4,1] * p.old[5,1] * (y0^(p.old[4,1] - 1))
    h2.1 <- p.old[4,1] * p.old[5,1] * (y1^(p.old[4,1] - 1))
    h2.2 <- p.old[4,1] * p.old[5,1] * (y2^(p.old[4,1] - 1))
    h2.3 <- p.old[4,1] * p.old[5,1] * (y3^(p.old[4,1] - 1))
    
    h3.0 <- p.old[6,1] * p.old[7,1] * (y0^(p.old[6,1] - 1))
    h3.1 <- p.old[6,1] * p.old[7,1] * (y1^(p.old[6,1] - 1))
    h3.2 <- p.old[6,1] * p.old[7,1] * (y2^(p.old[6,1] - 1))
    h3.3 <- p.old[6,1] * p.old[7,1] * (y3^(p.old[6,1] - 1))
    
    # Cumulative hazard functions
    H1.0 <- p.old[3,1] * (y0^p.old[2,1])
    H1.1 <- p.old[3,1] * (y1^p.old[2,1])
    H1.2 <- p.old[3,1] * (y2^p.old[2,1])
    H1.3 <- p.old[3,1] * (y3^p.old[2,1])
    
    H2.0 <- p.old[5,1] * (y0^p.old[4,1])
    H2.1 <- p.old[5,1] * (y1^p.old[4,1])
    H2.2 <- p.old[5,1] * (y2^p.old[4,1])
    H2.3 <- p.old[5,1] * (y3^p.old[4,1])
    
    H3.0 <- p.old[7,1] * (y0^p.old[6,1])
    H3.1 <- p.old[7,1] * (y1^p.old[6,1])
    H3.2 <- p.old[7,1] * (y2^p.old[6,1])
    H3.3 <- p.old[7,1] * (y3^p.old[6,1])
    

    # M-Step
    
    dt_model <- rpart(d ~ ., data = training_data, method = "class")
    
    
    cured_prob <- predict(dt_model, newdata = test_data, type = "prob")[, "0"]
    
    
    w0 <- ((1 - cured_prob) * exp(-H1.0 - H2.0 - H3.0)) / 
      ((cured_prob * exp(-H2.0 - H3.0)) + ((1 - cured_prob) * exp(-H1.0 - H2.0 - H3.0)))
    
    w2 <- ((1 - p.old[1,1]) * exp(-H1.2 - H2.2 - H3.2)) /
      ((p.old[1,1] * exp(-H2.2 - H3.2)) + ((1 - p.old[1,1]) * exp(-H1.2 - H2.2 - H3.2)))
    w3 <- ((1 - p.old[1,1]) * exp(-H1.3 - H2.3 - H3.3)) /
      ((p.old[1,1] * exp(-H2.3 - H3.3)) + ((1 - p.old[1,1]) * exp(-H1.3 - H2.3 - H3.3)))
    
    
    Q1 <- function(par1 = c(a1, l1)) {
      haz1.1 = par1[1] * par1[2] * (y1^(par1[1] - 1))
      chaz1.0 = par1[2] * (y0^par1[1])
      chaz1.1 = par1[2] * (y1^par1[1])
      chaz1.2 = par1[2] * (y2^par1[1])
      chaz1.3 = par1[2] * (y3^par1[1])
      res1 = sum(log(haz1.1) - chaz1.1) - sum(w0 * chaz1.0)
      return(-res1)
    }
    risk1.new <- optim(par = c(p.old[2, 1], p.old[3, 1]), fn = Q1, method = "Nelder-Mead")$par
    
    Q2 <- function(par2 = c(a2, l2)) {
      haz2.2 = par2[1] * par2[2] * (y2^(par2[1] - 1))
      chaz2.0 = par2[2] * (y0^par2[1])
      chaz2.1 = par2[2] * (y1^par2[1])
      chaz2.2 = par2[2] * (y2^par2[1])
      chaz2.3 = par2[2] * (y3^par2[1])
      res2 = sum(log(haz2.2)) - sum(chaz2.1) - sum(w0 * chaz2.0)
      return(-res2)
    }
    risk2.new <- optim(par = c(p.old[4, 1], p.old[5, 1]), fn = Q2, method = "Nelder-Mead")$par
    
    Q3 <- function(par3 = c(a3, l3)) {
      haz3.3 = par3[1] * par3[2] * (y3^(par3[1] - 1))
      chaz3.0 = par3[2] * (y0^par3[1])
      chaz3.1 = par3[2] * (y1^par3[1])
      chaz3.2 = par3[2] * (y2^par3[1])
      chaz3.3 = par3[2] * (y3^par3[1])
      res3 = sum(log(haz3.3)) - sum(chaz3.1) - sum(w0 * chaz3.0)
      return(-res3)
    }
    risk3.new <- optim(par = c(p.old[6, 1], p.old[7, 1]), fn = Q3, method = "Nelder-Mead")$par
    
    p.new = matrix(c(mean(cured_prob), risk1.new, risk2.new, risk3.new))
    
    iter = iter + 1
    continue = (max(abs(p.new - p.old)) > tol) && (iter < maxit)
    p.old <- p.new
  }

  
  list(parameters = p.new, iterations = iter, cured_fraction = cured_prob)
  print(iter)
  print(p.new)
  #print(se.est)
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

#[1] 6 iterations
#[,1]
#[1,] 0.28183867
#[2,] 0.89473888
#[3,] 0.22694818
#[4,] 0.64982165
#[5,] 0.01293865
#[6,] 0.93945217
#[7,] 0.02185218


# Subgroup analysis with respect to race 

brst.data = read.table("Brstdata.txt",header=T)

for(i in 1:length(brst.data$T)){
  if(brst.data$Cause[i]=="Other_Cancers"){
    brst.data$D[i] = 2
  }else if(brst.data$Cause[i]=="Other_Causes"){
    brst.data$D[i] = 3
  }
}

J = rep(NA,length(brst.data$T))
for(j in 1:length(brst.data$T)){
  if(brst.data$Cause[j]=="Breast_Cancer"){
    J[j] = 1
  }else{
    J[j] = 0
  }
}

brst.data$Race = rep(NA,length(brst.data$T))
for(j in 1:length(brst.data$T)){
  if(brst.data$race_asian[j]==1 & brst.data$race_black[j]==0){
    brst.data$Race[j] = "Asian"
  }else if(brst.data$race_asian[j]==0 & brst.data$race_black[j]==1){
    brst.data$Race[j] = "Black"
  }else{
    brst.data$Race[j] = "White"
  }
}

brst.data$Type = rep(NA,length(brst.data$T))
for(j in 1:length(brst.data$T)){
  if(brst.data$luminal_b[j]==1){
    brst.data$Type[j] = "LUM-B"
  }else if(brst.data$her2_enrch[j]==1){
    brst.data$Type[j] = "HER2-ER"
  }else if(brst.data$unkn[j]==1){
    brst.data$Type[j] = "UN"
  }else if(brst.data$tripl_neg[j]==1){
    brst.data$Type[j] = "TRI-NEG"
  }else{
    brst.data$Type[j] = "LUM-A"
  }
}

data.race = data.frame(Y=brst.data$T/12,J,Race=brst.data$Race,d=brst.data$D)

#head(data.race)

KM.race = survfit(Surv(Y,J)~Race,data=data.race)
ggsurvplot(KM.race,data=data.race,xlab="Time (years)")

data.black = data.race[data.race$Race=="Black",]
data.asian = data.race[data.race$Race=="Asian",]
data.white = data.race[data.race$Race=="White",]

KM.black = survfit(Surv(Y,J)~1,data=data.black)
KM.asian = survfit(Surv(Y,J)~1,data=data.asian)
KM.white = survfit(Surv(Y,J)~1,data=data.white)

# c.black.init = 0.135; c.white.init = 0.314; c.asian.init = 0.396

library(BB)

# Risk 1:

data.asian1 = data.asian[data.asian$D==1,1]
mean(data.asian1) # 2.52
var(data.asian1) # 3.61

# Risk 2:

data.asian2 = data.asian[data.asian$D==2,1]
mean(data.asian2) # 3.97
var(data.asian2) # 4.01


# Risk 3:

data.asian3 = data.asian[data.asian$D==3,1]
mean(data.asian3) # 3.24
var(data.asian3) # 6.18


f1 = function(v=c(alp,lam)){
  f = rep(NA,2)
  f[1] =  ((1/(v[2]^(1/v[1])))*gamma(1+(1/v[1]))) - 3.24
  f[2] =  ((1/(v[2]^(1/v[1])))^2)*((gamma(1+(2/v[1])))-((gamma(1+(1/v[1])))^2)) - 6.18
  f
}

BBsolve(par=c(0.5,0.5),fn=f1)

##############################################################################################################################################

# Race Black

# alp1.init = 1.14; lam1.init = 0.42
# alp2.init = 0.79; lam2.init = 0.76
# alp3.init = 0.97; lam3.init = 0.49

#data = data.frame(y=data.black$Y,d=data.black$D)

#head(brst.data)

#CONSIDERING ALL COVARIATES

#####################################################################################################################################################

#Race Black
c.init = 0.135
alp1.init = 1.14
lam1.init = 0.42
alp2.init = 0.79
lam2.init = 0.76
alp3.init = 0.97
lam3.init = 0.49

brst.data.black <- subset(brst.data.new, race_black == 1)

EM.MCCR3(brst.data.black,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] 6
#[,1]
#[1,] 0.11283215
#[2,] 0.81059240
#[3,] 0.27397663
#[4,] 0.74008845
#[5,] 0.01161256
#[6,] 0.84830386
#[7,] 0.01978104

###########################################################################################################################################
# Race White

# alp1.init = 1.33; lam1.init = 0.26
# alp2.init = 0.99; lam2.init = 0.50
# alp3.init = 1.36; lam3.init = 0.23

#data = data.frame(y=data.white$Y,d=data.white$D)

c.init = 0.314
alp1.init = 1.33
lam1.init = 0.26
alp2.init = 0.99
lam2.init = 0.50
alp3.init = 1.36
lam3.init = 0.23

brst.data.white <- subset(brst.data.new, race_white == 1)

EM.MCCR3(brst.data.white,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] 6
#[,1]
#[1,] 0.25560349
#[2,] 0.87427706
#[3,] 0.21790935
#[4,] 0.63086710
#[5,] 0.01187222
#[6,] 0.92533822
#[7,] 0.01691911

#############################################################################################################################

# Race Asian

# alp1.init = 1.34; lam1.init = 0.26
# alp2.init = 2.08; lam2.init = 0.04 
# alp3.init = 1.31; lam3.init = 0.19


#data = data.frame(y=data.asian$Y,d=data.asian$D)


c.init = 0.396
alp1.init = 1.34
lam1.init = 0.26
alp2.init = 2.08
lam2.init = 0.04
alp3.init = 1.31
lam3.init = 0.19

brst.data.asian <- subset(brst.data.new, race_asian == 1)

EM.MCCR3(brst.data.asian,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results:
#[1] 6
#[,1]
#[1,] 0.257575758
#[2,] 0.890349402
#[3,] 0.192416958
#[4,] 1.620762834
#[5,] 0.003888169
#[6,] 1.249929009
#[7,] 0.014994453

########################################


# Subgroup analysis with respect to Type of breast cancer

data.type = data.frame(Y=brst.data$T/12,J,Type=brst.data$Type,D=brst.data$D)
KM.type = survfit(Surv(Y,J)~Type,data=data.type)
ggsurvplot(KM.type,data=data.type,xlab="Time (years)")

data.luma = data.type[data.type$Type=="LUM-A",]
data.lumb = data.type[data.type$Type=="LUM-B",]
data.her = data.type[data.type$Type=="HER2-ER",]
data.tri = data.type[data.type$Type=="TRI-NEG",]
data.un = data.type[data.type$Type=="UN",]


KM.luma = survfit(Surv(Y,J)~1,data=data.luma)
KM.lumb = survfit(Surv(Y,J)~1,data=data.lumb)
KM.her = survfit(Surv(Y,J)~1,data=data.her)
KM.tri = survfit(Surv(Y,J)~1,data=data.tri)
KM.un = survfit(Surv(Y,J)~1,data=data.un)

#c.luma=0.281; c.lumb=0.431; c.her=0.380; c.tri=0.226; c.un=0.157

library(BB)

# Risk 1:
data.un1 = data.un[data.un$D==1,1]
mean(data.un1) # 2.02
var(data.un1) #  3.26

# Risk 2:
data.un2 = data.un[data.un$D==2,1]
mean(data.un2) # 1.12
var(data.un2) #  4.26

# Risk 3:
data.un3 = data.un[data.un$D==3,1]
mean(data.un3) # 3.07
var(data.un3) # 4.89

f1 = function(v=c(alp,lam)){
  f = rep(NA,2)
  f[1] =  ((1/(v[2]^(1/v[1])))*gamma(1+(1/v[1]))) - 3.07
  f[2] =  ((1/(v[2]^(1/v[1])))^2)*((gamma(1+(2/v[1])))-((gamma(1+(1/v[1])))^2)) - 4.89
  f
}
BBsolve(par=c(0.5,0.5),fn=f1)

####################################################################################################################################################

# LUM-A

brst.data.luma <- subset(brst.data.new, luminal_a == 1)

c.init = 0.281
alp1.init = 1.49
lam1.init = 0.17
alp2.init = 1.59
lam2.init = 0.14
alp3.init = 1.37
lam3.init = 0.22

EM.MCCR3(brst.data.luma,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] 6
#[,1]
#[1,] 0.254309018
#[2,] 0.978947912
#[3,] 0.184105689
#[4,] 1.390431102
#[5,] 0.003239407
#[6,] 1.015996346
#[7,] 0.017234042

##################################################################################################################################

# LUM-B

c.init = 0.431
alp1.init = 1.43
lam1.init = 0.22
alp2.init = 3.09
lam2.init = 0.03
alp3.init = 1.25
lam3.init = 0.27

brst.data.lumb <- subset(brst.data.new, luminal_b == 1)

EM.MCCR3(brst.data.lumb,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] 7
#[,1]
#[1,] 0.432988363
#[2,] 0.901797300
#[3,] 0.229996582
#[4,] 1.330212358
#[5,] 0.005025655
#[6,] 0.606073636
#[7,] 0.023643003


################################################################################################################

# HER2-ER

c.init = 0.380
alp1.init = 1.44
lam1.init = 0.27
alp2.init = 0.78
lam2.init = 1.00
alp3.init = 0.98
lam3.init = 0.35

brst.data.her2 <- subset(brst.data.new, her2_enrch == 1)

EM.MCCR3(brst.data.her2,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] 6
#[,1]
#[1,] 0.380995787
#[2,] 0.939429489
#[3,] 0.211528249
#[4,] 1.017305892
#[5,] 0.009686423
#[6,] 0.892428937
#[7,] 0.027429371

####################################################################################################################################

# TRI-NEG

c.init = 0.226
alp1.init = 1.01
lam1.init = 0.67
alp2.init = 1.04
lam2.init = 0.42
alp3.init = 0.85
lam3.init = 0.63

brst.data.trig <- subset(brst.data.new, tripl_neg == 1)

EM.MCCR3(brst.data.trig,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] 7
#[,1]
#[1,] 0.23154762
#[2,] 0.68888979
#[3,] 0.37159333
#[4,] 0.76465161
#[5,] 0.01513079
#[6,] 0.64417738
#[7,] 0.02665233

##############################################################################################################

# UN

c.init = 0.157
alp1.init = 1.12
lam1.init = 0.43
alp2.init = 0.57
lam2.init = 1.22
alp3.init = 1.41
lam3.init = 0.18

brst.data.unkn <- subset(brst.data.new, unkn == 1)

EM.MCCR3(brst.data.unkn,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] 5
#[,1]
#[1,] 0.14814815
#[2,] 0.74805626
#[3,] 0.29436863
#[4,] 0.43048614
#[5,] 0.04588792
#[6,] 0.97632903
#[7,] 0.01327434

####################################################################################################################################################

# Finding the CIF

library(calculus)


c = 0.366
alp1 = 1.301
lam1 = 0.287
alp2 = 0.455 
lam2 = 0.024
alp3 = 0.844   
lam3 = 0.016

# risk 1

f1 = function(t){
  h1 = alp1*lam1*t^{alp1-1}
  H1 = lam1*t^alp1
  H2 = lam2*t^alp2
  H3 = lam3*t^alp3
  res1 = h1*exp(-H1-H2-H3)
  return(res1)
}

# risk 2

f22 = function(t){
  h2 = alp2*lam2*t^{alp2-1}
  H2 = lam2*t^alp2
  H3 = lam3*t^alp3
  res22 = h2*exp(-H2-H3)
  return(res22)
}

f2 = function(t){
  h2 = alp2*lam2*t^{alp2-1}
  H1 = lam1*t^alp1
  H2 = lam2*t^alp2
  H3 = lam3*t^alp3
  res2 = h2*exp(-H1-H2-H3)
  return(res2)
}

# risk 3

f3 = function(t){
  h3 = alp3*lam3*t^{alp3-1}
  H1 = lam1*t^alp1
  H2 = lam2*t^alp2
  H3 = lam3*t^alp3
  res3 = h3*exp(-H1-H2-H3)
  return(res3)
}

f33 = function(t){
  h3 = alp3*lam3*t^{alp3-1}
  H2 = lam2*t^alp2
  H3 = lam3*t^alp3
  res33 = h3*exp(-H2-H3)
  return(res33)
}


x = seq(0.1,8,by=0.1)
i1 = rep(NA,length(x))
i2 = rep(NA,length(x))
i3 = rep(NA,length(x))
i22 = rep(NA,length(x))
i33 = rep(NA,length(x))

for(i in 1:length(x)){
  i1[i] = integral(f1, bounds = list(t = c(0,x[i])))$value
  i2[i] = integral(f2, bounds = list(t = c(0,x[i])))$value
  i3[i] = integral(f3, bounds = list(t = c(0,x[i])))$value
  i22[i] = integral(f22, bounds = list(t = c(0,x[i])))$value
  i33[i] = integral(f33, bounds = list(t = c(0,x[i])))$value
}



plot(x,(1-c)*i1,xlab="Time (in years)",ylab="Cumulative probability of death",
     ylim=c(0,1),type="l",lty=1,col="red",main="HER2-ER")
lines(x,(c*i22)+((1-c)*i2),lty=2,col="purple")
lines(x,(c*i33)+((1-c)*i3),lty=3,col="blue")
legend("topleft",legend=c("Breast Cancer","Other Cancers","Other Causes"),lty=c(1,2,3),col=c("red","purple","blue"))


# Plots of overall survival stratified by race


t = seq(0,8,by=0.1)

c1 = 0.088
alp11 = 0.995
lam11 = 0.377
alp21 = 0.666 
lam21 = 0.021
alp31 = 0.866   
lam31 = 0.034 

H11 = lam11*t^alp11
H21 = lam21*t^alp21
H31 = lam31*t^alp31

c2 = 0.219
alp12 = 1.050
lam12 = 0.251
alp22 = 0.671 
lam22 = 0.012
alp32 = 0.996   
lam32 = 0.014 

H12 = lam12*t^alp12
H22 = lam22*t^alp22
H32 = lam32*t^alp32

c3 = 0.326
alp13 = 1.104
lam13 = 0.231
alp23 = 1.820  
lam23 = 0.003
alp33 = 1.017    
lam33 = 0.022

H13 = lam13*t^alp13
H23 = lam23*t^alp23
H33 = lam33*t^alp33


S1 = (c1*exp(-H21-H31)) + ((1-c1)*exp(-H11-H21-H31))
S2 = (c2*exp(-H22-H32)) + ((1-c2)*exp(-H12-H22-H32))
S3 = (c3*exp(-H23-H33)) + ((1-c3)*exp(-H13-H23-H33))

plot(t,S1,xlab="Time (in years)",ylab="Overall survival probability"
     ,type="l",lty=1,col="red")
lines(t,S2,lty=2,col="blue")
lines(t,S3,lty=3,col="purple")

legend("topright",legend=c("Black","White","Asian"),lty=c(1,2,3),
       col=c("red","blue","purple"))



# Plots of overall survival stratified by cancer type


t = seq(0,8,by=0.1)

c1 = 0.104
alp11 = 1.122
lam11 = 0.165
alp21 = 1.254 
lam21 = 0.004
alp31 = 1.029   
lam31 = 0.018 

H11 = lam11*t^alp11
H21 = lam21*t^alp21
H31 = lam31*t^alp31

c2 = 0.105
alp12 = 0.962
lam12 = 0.388
alp22 = 0.506 
lam22 = 0.050
alp32 = 1.190   
lam32 = 0.014 

H12 = lam12*t^alp12
H22 = lam22*t^alp22
H32 = lam32*t^alp32


c3 = 0.225 
alp13 = 1.067
lam13 = 0.588
alp23 = 0.959 
lam23 = 0.012
alp33 = 0.809   
lam33 = 0.025

H13 = lam13*t^alp13
H23 = lam23*t^alp23
H33 = lam33*t^alp33


c4 = 0.309 
alp14 = 1.010
lam14 = 0.229
alp24 = 1.211  
lam24 = 0.003
alp34 = 0.763   
lam34 = 0.014 

H14 = lam14*t^alp14
H24 = lam24*t^alp24
H34 = lam34*t^alp34


c5 = 0.366 
alp15 = 1.301
lam15 = 0.287
alp25 = 0.455  
lam25 = 0.024
alp35 = 0.844  
lam35 = 0.016 

H15 = lam15*t^alp15
H25 = lam25*t^alp25
H35 = lam35*t^alp35

S1 = (c1*exp(-H21-H31)) + ((1-c1)*exp(-H11-H21-H31))
S2 = (c2*exp(-H22-H32)) + ((1-c2)*exp(-H12-H22-H32))
S3 = (c3*exp(-H23-H33)) + ((1-c3)*exp(-H13-H23-H33))
S4 = (c4*exp(-H24-H34)) + ((1-c4)*exp(-H14-H24-H34))
S5 = (c5*exp(-H25-H35)) + ((1-c5)*exp(-H15-H25-H35))

plot(t,S1,xlab="Time (in years)",ylab="Overall survival probability"
     ,type="l",lty=1,col="red",ylim=c(0,1))
lines(t,S2,lty=2,col="blue")
lines(t,S3,lty=3,col="purple")
lines(t,S4,lty=4,col="green")
lines(t,S5,lty=5,col="orange")

legend("topright",legend=c("LUM-A","UN","TRI-NEG","LUM-B","HER2-ER"),lty=c(1,2,3,4,5),
       col=c("red","blue","purple","green","orange"))



####################################################################################################
#PLOTTING MY RESULTS

# Time vector
t <- seq(0, 20, by = 0.1)

# Black group parameters
c1 <- 0.11283215
alp11 <- 0.81059240
lam11 <- 0.27397663
alp21 <- 0.74008845
lam21 <- 0.01161256
alp31 <- 0.84830386
lam31 <- 0.01978104

# Hazard functions for Black group
H11 <- lam11 * t^alp11
H21 <- lam21 * t^alp21
H31 <- lam31 * t^alp31

# White group parameters
c2 <- 0.25560349
alp12 <- 0.87427706
lam12 <- 0.21790935
alp22 <- 0.63086710
lam22 <- 0.01187222
alp32 <- 0.92533822
lam32 <- 0.01691911

# Hazard functions for White group
H12 <- lam12 * t^alp12
H22 <- lam22 * t^alp22
H32 <- lam32 * t^alp32

# Asian group parameters
c3 <- 0.257575758
alp13 <- 0.890349402
lam13 <- 0.192416958
alp23 <- 1.620762834
lam23 <- 0.003888169
alp33 <- 1.249929009
lam33 <- 0.014994453

# Hazard functions for Asian group
H13 <- lam13 * t^alp13
H23 <- lam23 * t^alp23
H33 <- lam33 * t^alp33

# Survival functions for each group
S1 <- (c1 * exp(-H21 - H31)) + ((1 - c1) * exp(-H11 - H21 - H31))
S2 <- (c2 * exp(-H22 - H32)) + ((1 - c2) * exp(-H12 - H22 - H32))
S3 <- (c3 * exp(-H23 - H33)) + ((1 - c3) * exp(-H13 - H23 - H33))

# Plot the survival curves
plot(t, S1, xlab = "Time (in years)", ylab = "Overall survival probability", type = "l", lty = 1, col = "red", main = "Survival Curves by Race")
lines(t, S2, lty = 2, col = "blue")
lines(t, S3, lty = 3, col = "purple")

# Add legend to the plot
legend("topright", legend = c("Black", "White", "Asian"), lty = c(1, 2, 3), col = c("red", "blue", "purple"))

############################################################################################################



# Time vector
t <- seq(0, 25, by = 0.1)

# LUM-A group parameters
c1 <- 0.254309018
alp11 <- 0.978947912
lam11 <- 0.184105689
alp21 <- 1.390431102
lam21 <- 0.003239407
alp31 <- 1.015996346
lam31 <- 0.017234042

# Hazard functions for LUM-A group
H11 <- lam11 * t^alp11
H21 <- lam21 * t^alp21
H31 <- lam31 * t^alp31

# LUM-B group parameters
c2 <- 0.432988363
alp12 <- 0.901797300
lam12 <- 0.229996582
alp22 <- 1.330212358
lam22 <- 0.005025655
alp32 <- 0.606073636
lam32 <- 0.023643003

# Hazard functions for LUM-B group
H12 <- lam12 * t^alp12
H22 <- lam22 * t^alp22
H32 <- lam32 * t^alp32

# HER2-ER group parameters
c3 <- 0.380995787
alp13 <- 0.939429489
lam13 <- 0.211528249
alp23 <- 1.017305892
lam23 <- 0.009686423
alp33 <- 0.892428937
lam33 <- 0.027429371

# Hazard functions for HER2-ER group
H13 <- lam13 * t^alp13
H23 <- lam23 * t^alp23
H33 <- lam33 * t^alp33

# TRI-NEG group parameters
c4 <- 0.23154762
alp14 <- 0.68888979
lam14 <- 0.37159333
alp24 <- 0.76465161
lam24 <- 0.01513079
alp34 <- 0.64417738
lam34 <- 0.02665233

# Hazard functions for TRI-NEG group
H14 <- lam14 * t^alp14
H24 <- lam24 * t^alp24
H34 <- lam34 * t^alp34

# UN group parameters
c5 <- 0.14814815
alp15 <- 0.74805626
lam15 <- 0.29436863
alp25 <- 0.43048614
lam25 <- 0.04588792
alp35 <- 0.97632903
lam35 <- 0.01327434

# Hazard functions for UN group
H15 <- lam15 * t^alp15
H25 <- lam25 * t^alp25
H35 <- lam35 * t^alp35

# Survival functions for each group
S1 <- (c1 * exp(-H21 - H31)) + ((1 - c1) * exp(-H11 - H21 - H31))
S2 <- (c2 * exp(-H22 - H32)) + ((1 - c2) * exp(-H12 - H22 - H32))
S3 <- (c3 * exp(-H23 - H33)) + ((1 - c3) * exp(-H13 - H23 - H33))
S4 <- (c4 * exp(-H24 - H34)) + ((1 - c4) * exp(-H14 - H24 - H34))
S5 <- (c5 * exp(-H25 - H35)) + ((1 - c5) * exp(-H15 - H25 - H35))

# Plot the survival curves
plot(t, S1, xlab = "Time (in years)", ylab = "Overall survival probability", type = "l", lty = 1, col = "red", main = "Survival Curves by Cancer Type")
lines(t, S2, lty = 2, col = "green")
lines(t, S3, lty = 3, col = "orange")
lines(t, S4, lty = 4, col = "purple")
lines(t, S5, lty = 5, col = "blue")

# Add legend to the plot
legend("topright", legend = c("LUM-A", "LUM-B", "HER2-ER", "TRI-NEG", "UN"), lty = c(1, 2, 3, 4, 5), col = c("red", "green", "orange", "purple", "blue"))


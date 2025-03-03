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
    
    #w0 = ((1-p.old[1,1])*exp(-H1.0-H2.0-H3.0))/((p.old[1,1]*exp(-H2.0-H3.0)) + ((1-p.old[1,1])*exp(-H1.0-H2.0-H3.0))) 
    #w2 = ((1-p.old[1,1])*exp(-H1.2-H2.2-H3.2))/((p.old[1,1]*exp(-H2.2-H3.2)) + ((1-p.old[1,1])*exp(-H1.2-H2.2-H3.2))) 
    #w3 = ((1-p.old[1,1])*exp(-H1.3-H2.3-H3.3))/((p.old[1,1]*exp(-H2.3-H3.3)) + ((1-p.old[1,1])*exp(-H1.3-H2.3-H3.3)))
    
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
  print(p.old)
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

#brst.data.new


c.init = 0.289
alp1.init = 1.29
lam1.init = 0.29
alp2.init = 1.02
lam2.init = 0.46
alp3.init = 1.26
lam3.init = 0.26

EM.MCCR3(brst.data.new,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results

#[1] "iteration"
#[1] 35
#[1] 0.204402264 1.032766447 0.273091477 0.709014669 0.013384210 0.947013010 0.017734819 0.024220248 0.036541631 0.015333101 0.088791114
#[12] 0.002565078 0.083279863 0.002946646


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

data.race = data.frame(Y=brst.data$T/12,J,Race=brst.data$Race,D=brst.data$D)

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


# Race Black

# alp1.init = 1.14; lam1.init = 0.42
# alp2.init = 0.79; lam2.init = 0.76
# alp3.init = 0.97; lam3.init = 0.49

data = data.frame(y=data.black$Y,d=data.black$D)

c.init = 0.135
alp1.init = 1.14
lam1.init = 0.42
alp2.init = 0.79
lam2.init = 0.76
alp3.init = 0.97
lam3.init = 0.49

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] "iteration"
#[1] 28
#[1] 0.088380123 0.994808862 0.377142158 0.665895983 0.020979589 0.865604803 0.034543139 0.036716548 0.070087618 0.039451081 0.181286365 0.007935594
#[13] 0.155538221 0.010225787


# Race White

# alp1.init = 1.33; lam1.init = 0.26
# alp2.init = 0.99; lam2.init = 0.50
# alp3.init = 1.36; lam3.init = 0.23

data = data.frame(y=data.white$Y,d=data.white$D)

c.init = 0.314
alp1.init = 1.33
lam1.init = 0.26
alp2.init = 0.99
lam2.init = 0.50
alp3.init = 1.36
lam3.init = 0.23

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] "iteration"
#[1] 40
#[1] 0.218755569 1.049604421 0.250918856 0.671087595 0.012522370 0.995602144 0.013838098 0.030290756 0.044231354 0.017029247 0.101874213 0.002833254  
#[13] 0.107230235 0.002957499


# Race Asian

# alp1.init = 1.34; lam1.init = 0.26
# alp2.init = 2.08; lam2.init = 0.04 
# alp3.init = 1.31; lam3.init = 0.19


data = data.frame(y=data.asian$Y,d=data.asian$D)

c.init = 0.396
alp1.init = 1.34
lam1.init = 0.26
alp2.init = 2.08
lam2.init = 0.04
alp3.init = 1.31
lam3.init = 0.19

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] "iteration"
#[1] 29
#[1] 0.325692036 1.103772572 0.231449825 1.820595227 0.002691591 1.016861220 0.022025064 0.097447852 0.165438076 0.057514698 0.679856373 0.003546364
#[13] 0.280197525 0.012389639


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

# LUM-A

data = data.frame(y=data.luma$Y,d=data.luma$D)

c.init = 0.281
alp1.init = 1.49
lam1.init = 0.17
alp2.init = 1.59
lam2.init = 0.14
alp3.init = 1.37
lam3.init = 0.22

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] "iteration"
#[1] 115
#[1] 0.104222536 1.121769997 0.165421155 1.254496884 0.003768942 1.029101623 0.017680682 0.068151680 0.065424849 0.017450162 0.259144462 0.001875976
#[13] 0.120805512 0.004196136


# LUM-B

data = data.frame(y=data.lumb$Y,d=data.lumb$D)

c.init = 0.431
alp1.init = 1.43
lam1.init = 0.22
alp2.init = 3.09
lam2.init = 0.03
alp3.init = 1.25
lam3.init = 0.27

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

# results
#[1] "iteration"
#[1] 49
#[1] 0.308830699 1.010266970 0.228783565 1.211568544 0.002610899 0.763143939 0.014036371 0.092377731 0.116783291 0.043776835 0.533208095 0.002776264
#[13] 0.219669969 0.006732345


# HER2-ER

data = data.frame(y=data.her$Y,d=data.her$D)

c.init = 0.380
alp1.init = 1.44
lam1.init = 0.27
alp2.init = 0.78
lam2.init = 1.00
alp3.init = 0.98
lam3.init = 0.35

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] "iteration"
#[1] 8
#[1] 0.366439240 1.300857871 0.286725862 0.454605669 0.024444242 0.844201306 0.015974418 0.044655323 0.135965868 0.050709464 0.172966208 0.011532396
#[13] 0.279199096 0.009301916


# TRI-NEG

data = data.frame(y=data.tri$Y,d=data.tri$D)

c.init = 0.226
alp1.init = 1.01
lam1.init = 0.67
alp2.init = 1.04
lam2.init = 0.42
alp3.init = 0.85
lam3.init = 0.63

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] "iteration"
#[1] 5
#[1] 0.224909359 1.067240138 0.587983373 0.959309751 0.012413281 0.808868394 0.024671415 0.028122307 0.064176394 0.053216665 0.265006532 0.006345383
#[13] 0.181107972 0.008917510



# UN

data = data.frame(y=data.un$Y,d=data.un$D)

c.init = 0.157
alp1.init = 1.12
lam1.init = 0.43
alp2.init = 0.57
lam2.init = 1.22
alp3.init = 1.41
lam3.init = 0.18

EM.MCCR3(data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init,alpha3=alp3.init,lambda3=lam3.init)

#results
#[1] "iteration"
#[1] 29
#[1] 0.10558595 0.96228586 0.38772107 0.50608746 0.05009890 1.19014910 0.01376427 0.04760868 0.08034249 0.04921023 0.11632126 0.01425732 0.29638470
#[14] 0.00749338

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


######################################### 
#Data Analysis                          #
#Version: 04/11/2019                    #
#Author: Changee Chang, Jeong Hoon Jang #
#########################################

options(scipen=9)
#Install Packages




##### Call Functions #####
source("Nuclear_Function_05062018.R")
source("Simulation_10-12_Function.R")
source("Misc_Function.R")




##### Read in the data ######
renaldata <- read.delim("Pseudo_Dataset.txt",header=T,sep="")

# Get the time point for baseline and diuretic renogram curves
Baseline_Time_Interval <- as.numeric(substring(colnames(renaldata)[10:68],2))
Diuretic_Time_Interval <- as.numeric(substring(colnames(renaldata)[69:108],2))


#### Turn it into an analyzable form ####

### Read in each variable ###
n <- 216 # Size of total dataset
n1 <- 140 # Size of training dataset
n2 <- 76 # Size of testing dataset

y1 <- t(as.matrix(renaldata[,10:68])) # Read in Baseline renogram curves
y2 <- t(as.matrix(renaldata[,69:108])) # Read in diuretic renogram curves

X <- t(as.matrix(renaldata[,c(5,6)])) # Read in clinical variables

W <- t(as.matrix(renaldata[,c(7:9)])) # Read in Expert Ratings




### Basic Descriptive Analysis ###
# Age
range(X[1,])
mean(X[1,])

# Gender (0: male; 1: Female)
table(X[2,])


# Unanimous Classification
con <- rep(0,n)
for (i in 1:n) {
  if (W[1,i] == 0 & W[2,i] == 0 & W[3,i] == 0) {
    con[i] = 0
  } else if (W[1,i] == 2 & W[2,i] == 2 & W[3,i] == 2) {
    con[i] = 2 }
  else {
    con[i] = 1 
  }
}

table(con)




### Normalize Renogram curves and clinical variables ###
# Normalize curves (each baseline curve has maximum 25)
for(i in 1:n) {
  y2[,i] <- y2[,i]/max(y1[,i])*25
  y1[,i] <- y1[,i]/max(y1[,i])*25
}

X[1,] <- (X[1,] - mean(X[1,]))/sd(X[1,]) # Normalize age










##Divide into training and testing datasets ##
set.seed(90220)
kidney_index <- 1:n
train <- sort(sample(kidney_index,n1))
#93 %in% train
test <- kidney_index[-train]



y1_training <- y1[,train]
y2_training <- y2[,train]
y1_testing <- y1[,test]
y2_testing <- y2[,test]

X_training <- X[,train]
X_testing <- X[,test]

W_training <- W[,train]
W_testing <- W[,test]





# # Unanimous Classification for testing dataset
# con2 <- rep(0,n2)
# for (i in 1:n2) {
#   if (W_testing[1,i] == 0 & W_testing[2,i] == 0 & W_testing[3,i] == 0) {
#     con2[i] = 0
#   } else if (W_testing[1,i] == 2 & W_testing[2,i] == 2 & W_testing[3,i] == 2) {
#     con2[i] = 2 }
#   else {
#     con2[i] = 1 
#   }
# }
# table(con2)








### Start Data Analysis using our algorithm ###
set.seed(227920) # set seed


nsample <- 10000 # number of mcmc sample
burn <- 2000 # burn-in

# Normalize the time intervals to be [0,1]
t1 <- as.numeric(Baseline_Time_Interval)/max(as.numeric(Baseline_Time_Interval))
t2 <- as.numeric(Diuretic_Time_Interval)/max(as.numeric(Diuretic_Time_Interval))

p <- 15 # number of basis functions
q <- c(6,6) # number of latent factors
nu <- 3 # 1st basis function parameter (controls smoothness)
psi <- 1/p # 2nd basis function parameter2 (equally spaced kernel)


# Other tuning parameters
ups <- 5
a <-  c(1.5,1.5)
a_sigma <- 1
b_sigma <- 1
tau <- 5
kappa <- 5



# Setting theta based on expert ratings
sum <- apply(W_training,2,sum)
mj <- rep(0, n1)
for (i in 1:n1) {
  if (sum[i] < 4) {
    mj[i] = 0
  }
  else {
    mj[i] = 1
  }
}
ep <- sum(mj)/n1
theta <- log(-(ep-1)/ep)




## Run the proposed algorithm on the training dataset ##
Result <- Nuclear(nsample,y1_training,y2_training,X_training,W_training,t1,t2,p,q,nu,psi,ups,a,a_sigma,b_sigma,tau,kappa,theta)














## Fitted Curves (Revision) ##
# Generate basis functions
b1 <- make_basis(p,nu,psi, t1,t2)$b1
b2 <- make_basis(p,nu,psi, t1,t2)$b2


#posterior beta mean
beta1mean <- apply(Result$beta[,29,1,(burn+1):nsample],1,mean)
beta2mean <- apply(Result$beta[,29,2,(burn+1):nsample],1,mean)

# Fitted and raw curve for each subject
par(mfrow=c(1,2))
plot(Baseline_Time_Interval, y1_training[,29], type='l',lwd=2,col=1, ylab="MAG3/1000",
     xlab="Time since injection (mins)", main="Baseline Renogram")
lines(Baseline_Time_Interval, b1%*%beta1mean, type='l',lwd=2,col=2)

plot(Diuretic_Time_Interval, y2_training[,29], type='l',lwd=2,col=1,ylim=c(0,50),ylab="MAG3/1000",
     xlab="Time since injection (mins)", main="Diuretic Renogram")
lines(Diuretic_Time_Interval, b2%*%beta2mean, type='l',lwd=2,col=2)








## Predict the outcome on the testing dataset ##
pr0 <- matrix(rep(0,n2*(nsample-burn)),nrow=n2) #store prediction results based on PR0 scheme
pr1 <- matrix(rep(0,n2*(nsample-burn)),nrow=n2) #store prediction results based on PR1 scheme

# Randomly select one out of the three expert ratings for each kidney in the testing dataset
# Selected expert ratings will be used for prediction using PR1 scheme
set.seed(1234)
g <- sample(c(1,2,3), n2, replace=TRUE) 

for(s in 1:n2) {
  pr0[s,] <- Nuclear_predict(fit=Result,bi=burn,y1_testing[,s],y2_testing[,s],X_testing[,s],W=NULL)
  pr1[s,] <- Nuclear_predict(fit=Result,bi=burn,y1_testing[,s],y2_testing[,s],X_testing[,s],W=W_testing[g[s],s])
  print(s)
}




## Predictive Probablities and upper & lower credible limits for the predictive probabilities
mpr0 <- apply(pr0,1,mean) # Mean Predictive Probabilities of PR0 results
upr0 <- apply(pr0,1,quantile, 0.975) # Upper limit of PR0 results
lpr0 <- apply(pr0,1,quantile, 0.025) # Lower limit of PR0 results

mpr1 <- apply(pr1,1,mean) # Mean Predictive Probabilities of PR1 results
upr1 <- apply(pr1,1,quantile, 0.975) # Upper limit of PR1 results
lpr1 <- apply(pr1,1,quantile, 0.025) # Lower limit of PR1 resilts



### Prediction result by imposing cutoff point at 0.5 ###
## PR0
rpr0 <- rep(0,n2) # Store binary prediction results under the PR0 scheme
for (i in 1:n2) {
  if (mpr0[i] < 0.5){
    rpr0[i] = 0
  } else {
    rpr0[i]=1
  }
}

# PR1
rpr1 <- rep(0,n2) # Store binary prediction results under the PR1 scheme
for (i in 1:n2) {
  if (mpr1[i] < 0.5){
    rpr1[i] = 0
  } else {
    rpr1[i]=1
  }
}



## Draw the mean curves for each obstruction status given by PR1 (Figure 4) ##
# PR1
y1rpc <- cbind(t(y1_testing),rpr1)
y10 <- y1rpc[which(rpr1=="0"),-length(y1rpc[1,])]
y11 <- y1rpc[which(rpr1=="1"),-length(y1rpc[1,])]
my10 <- apply(y10,2,mean)
my11 <- apply(y11,2,mean)

y2rpc <- cbind(t(y2_testing),rpr1)
y20 <- y2rpc[which(rpr1=="0"),-length(y2rpc[1,])]
y21 <- y2rpc[which(rpr1=="1"),-length(y2rpc[1,])]
my20 <- apply(y20,2,mean)
my21 <- apply(y21,2,mean)

par(mfrow=c(1,2))
plot(Baseline_Time_Interval,my10,col=1,ylim=c(0,30),type="l",lwd=2,lty=1,main="Baseline Renogram Curve",
     ylab="MAG3/1000",xlab="Time Since Injection (mins)")
lines(Baseline_Time_Interval,my11,col=2,lwd=2,lty=2)
plot(Diuretic_Time_Interval,my20,col=1,ylim=c(0,60),type="l",lwd=2,lty=1,main="Diuretic Renogram Curve",
     ylab="MAG3/1000",xlab="Time Since Injection (mins)")
lines(Diuretic_Time_Interval,my21,col=2,lwd=2,lty=2)
par(mfrow=c(1,1))









## Based on K-mean clustering on the curves ##
set.seed(227920)
mdata <- t(rbind(y1_training,y2_training))
kmean <-   kmeans(mdata, 2, nstart=100)
cc <- kmean[[1]]
kcenter <- kmean[[2]]

train <- Result$pC[,-(1:burn)]
train <- apply(train,1,mean)
train_result <- rep(0,n1)
for (i in 1:n1) {
  if (train[i] <= 0.5){
    train_result[i] = 0
  } else {
    train_result[i] =1
  }
}


rtable <- table(train_result,cc)
kc <- rep(0,n1)
if (rtable[1,1]+rtable[2,2] >= rtable[1,2]+rtable[2,1]) {
  for (i in 1:n1) {
    if (cc[i]==1){
      kc[i] = 0
    } else {
      kc[i] =1
    }
    
  }
  
} else if (rtable[1,1]+rtable[2,2] < rtable[1,2]+rtable[2,1]) {
  
  for (i in 1:n1) {
    if (cc[i]==1){
      kc[i] = 1
    } else {
      kc[i]=0
    }
  }
  A <- kcenter[1,] 
  B <- kcenter[2,]
  kcenter <- rbind(B,A)
}

kct <- rep(0,n2)
for (i in 1:n2) {
  diff1 <- euc.dist(c(y1_testing[,i], y2_testing[,i]),kcenter[1,])
  diff2 <- euc.dist(c(y1_testing[,i], y2_testing[,i]),kcenter[2,])
  if (diff1 < diff2) {
    kct[i] = 0
  }
  else {
    kct[i] = 1
  }
}










## Based on multivariate functional clustering clustering on the curves ##
basis1 <- create.bspline.basis(c(0,Baseline_Time_Interval[59]), nbasis=15)
basiseval1 <- eval.basis(Baseline_Time_Interval,basis1)

basis2 <- create.bspline.basis(c(0,Diuretic_Time_Interval[40]), nbasis=15)
basiseval2 <- eval.basis(Diuretic_Time_Interval,basis2)



# training 
fy1_training <- smooth.basis(Baseline_Time_Interval,y1_training,fdParobj=basis1)$fd
fy2_training <- smooth.basis(Diuretic_Time_Interval,y2_training,fdParobj=basis2)$fd
fy_training <- list(fy1_training,fy2_training)


set.seed(227920)
fcluster <- funHDDC(fy_training,K=2)
fclass <- fcluster$class-1

fmean01 <- basiseval1%*%fcluster$mu[1,1:15]
fmean11 <- basiseval1%*%fcluster$mu[2,1:15]

fmean02 <- basiseval2%*%fcluster$mu[1,16:30]
fmean12 <- basiseval2%*%fcluster$mu[2,16:30]


ftable_1 <- table(train_result,fclass)
fc <- c()
if (ftable_1[1,1]+ftable_1[2,2] >= ftable_1[1,2]+ftable_1[2,1]) {
  for (i in 1:n1) {
    if (fclass[i]==1){
      fc[i] = 1
    } else {
      fc[i] =0
    }
    
  }
  
} else if (ftable_1[1,1]+ftable_1[2,2] < ftable_1[1,2]+ftable_1[2,1]) {
  
  for (i in 1:n1) {
    if (fclass[i]==1){
      fc[i] = 0
    } else {
      fc[i]=1
    }
  }
  A <- fmean01 
  B <- fmean11
  C <- fmean02
  D <- fmean12
  
  fmean01 <- B
  fmean11 <- A
  fmean02 <- D
  fmean12 <- C
}


# testing
coefs_y1_testing <- smooth.basis(Baseline_Time_Interval,y1_testing,fdParobj=basis1)$fd$coefs
smooth_y1_testing <- basiseval1%*%coefs_y1_testing
coefs_y2_testing <- smooth.basis(Diuretic_Time_Interval,y2_testing,fdParobj=basis2)$fd$coefs
smooth_y2_testing <- basiseval2%*%coefs_y2_testing


fct <- c()
for (i in 1:n2) {
  fdiff0 <- hilbert_dist(smooth_y1_testing[,i],fmean01,t1, smooth_y2_testing[,i],fmean02,t2)
  fdiff1 <- hilbert_dist(smooth_y1_testing[,i],fmean11,t1, smooth_y2_testing[,i],fmean12,t2)
  if (fdiff0 < fdiff1) {
    fct[i] = 0
  }
  else {
    fct[i] = 1
  }
}




table(fct,kct)









### Artificial Gold Standard ###
nkct <- ifelse(kct==1,2,0)
nrpr1 <- ifelse(rpr1==1,2,0)
newdata <- cbind(W_testing[1,],W_testing[2,],W_testing[3,],nkct,nrpr1)
goldd <- rep(0,n2)
for (i in 1:n2) {
  goldd[i] <- majority(newdata[i,]) 
}
gold <- goldd[-(which(goldd==1))]
gold <- ifelse(gold==2,1,0)
length(gold)


diagnostic_performance(table(rpr0[-(which(goldd==1))],gold))
diagnostic_performance(table(rpr1[-(which(goldd==1))],gold))
diagnostic_performance2(table(W_testing[1,-(which(goldd==1))],gold))
diagnostic_performance2(table(W_testing[2,-(which(goldd==1))],gold))
diagnostic_performance2(table(W_testing[3,-(which(goldd==1))],gold))
diagnostic_performance(table(kct[-(which(goldd==1))],gold))
diagnostic_performance(table(fct[-(which(goldd==1))],gold))





#write.csv(newdata, file = "goldstandard.csv")






## Extract consensus expert ratings (0: if all 0; 1: if no consensus; 2: if all 2) fot testing data
con <- rep(0,n2)
for (i in 1:n2) {
  if (W_testing[1,i] == 0 & W_testing[2,i] == 0 & W_testing[3,i] == 0) {
    con[i] = 0
  } else if (W_testing[1,i] == 2 & W_testing[2,i] == 2 & W_testing[3,i] == 2) {
    con[i] = 2 }
  else {
    con[i] = 1 
  }
}
table(con)

## Predictive probabilities and credible intervals colored by consensus ratings (Figure 5) ##
par(mfrow=c(2,1))
ti <- 1:n2

# PR0
ppdata0 <- cbind(mpr0,upr0,lpr0)
plotCI(ti[which(con==2)],ppdata0[which(con==2),1],ui=ppdata0[which(con==2),2],li=ppdata0[which(con==2),3],col="red",main="CAD-P(0)",
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",ylim=c(0,1),xlim=c(1,n2)) # CI plot
plotCI(ti[which(con==0)],ppdata0[which(con==0),1],ui=ppdata0[which(con==0),2],li=ppdata0[which(con==0),3],col="blue",
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",add=TRUE) # CI plot
plotCI(ti[which(con==1)],ppdata0[which(con==1),1],ui=ppdata0[which(con==1),2],li=ppdata0[which(con==1),3],
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",add=TRUE) # CI plot
axis(1,1:n2,1:n2)
abline(h=0.5,lty=2,lwd=1,col="red")
text(8,0.55,"Cutoff = 0.5",col="red")


#PR1
ppdata1 <- cbind(mpr1,upr1,lpr1)
plotCI(ti[which(con==2)],ppdata1[which(con==2),1],ui=ppdata1[which(con==2),2],li=ppdata1[which(con==2),3],col="red",main="CAD-P(1)",
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",ylim=c(0,1),xlim=c(1,n2)) # CI plot
plotCI(ti[which(con==0)],ppdata1[which(con==0),1],ui=ppdata1[which(con==0),2],li=ppdata1[which(con==0),3],col="blue",
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",add=TRUE) # CI plot
plotCI(ti[which(con==1)],ppdata1[which(con==1),1],ui=ppdata1[which(con==1),2],li=ppdata1[which(con==1),3],
       sfrac=0.003,cex=0.7, ylab="Predictive Probability of Obstruction",xlab="Kidney Index",xaxt = "n",add=TRUE) # CI plot
axis(1,1:n2,1:n2)
abline(h=0.5,lty=2,lwd=1,col="red")
text(8,0.55,"Cutoff = 0.5",col="red")
par(mfrow=c(1,1))







## Analyze kidneys that were misclassified with the consensus ##
misclass <- rep(0, n2)
for (i in 1:n2) {
  if (con[i] == 2 & rpr1[i] == 0) {
    misclass[i] <- 1
  } else if (con[i]==0 & rpr1[i]== 1) {
    misclass[i] <- 2
  }  else {
    misclass[i] = 0
  }
}
# Consensus  = 0 but our method=1 (Figure 6)
case2 <- which(misclass==2)

par(mfrow=c(1,2))
plot(Baseline_Time_Interval,y1_testing[,case2[1]],col=1,type="l",ylim=c(0,30),lwd=2,lty=1,main="Baseline Renogram Curve",ylab="Normalized MAG3",xlab="Time (mins)")
for (i in 1:(length(case2)-1)) {
  lines(Baseline_Time_Interval, y1_testing[,case2[1+i]],lwd=2,lty=1+i,col=1+i)
}
plot(Diuretic_Time_Interval,y2_testing[,case2[1]],col=1,type="l",ylim=c(0,50),lwd=2,lty=1,main="Diuretic Renogram Curve",ylab="Normalized MAG3",xlab="Time (mins)")
for (i in 1:(length(case2)-1)) {
  lines(Diuretic_Time_Interval, y2_testing[,case2[1+i]],lwd=2,lty=1+i,col=1+i)
}
par(mfrow=c(1,1))






















#Fitted lambda
dim(Result$lam)
lam_un_1 <- Result$lam[,,1,1,nsample]
lam_un_2 <- Result$lam[,,1,2,nsample]

lam_ob_1 <- Result$lam[,,2,1,nsample]
lam_ob_2 <- Result$lam[,,2,2,nsample]

# Fitted gamma
dim(Result$gamma)
gamma_1 <- Result$gamma[,,1,nsample]
gamma_2 <- Result$gamma[,,2,nsample]


# Fitted Curves for 30 years old male
curve01 <- b1%*%lam_un_1%*%gamma_1%*%c((30 - 57.15741)/16.97504,0)
curve02 <- b2%*%lam_un_2%*%gamma_2%*%c((30 - 57.15741)/16.97504,0)

curve11 <- b1%*%lam_ob_1%*%gamma_1%*%c((30 - 57.15741)/16.97504,0)
curve12 <- b2%*%lam_ob_2%*%gamma_2%*%c((30 - 57.15741)/16.97504,0)

par(mfrow=c(1,2))
plot(Baseline_Time_Interval,curve01,col=1,type="l",lwd=2,lty=1,main="Baseline Renogram Curve",
     ylab="Normalized MAG3",xlab="Normalized Time (mins)",ylim=c(0,50))
lines(Baseline_Time_Interval, curve11, lty=1, lwd=2, type='l',col='red')
plot(Diuretic_Time_Interval,curve02,col=1,type="l",lwd=2,lty=1,main="Diuretic Renogram Curve",
     ylab="Normalized MAG3",xlab="Normalized Time (mins)",ylim=c(0,50))
lines(Diuretic_Time_Interval, curve12, lty=1, lwd=2, type='l',col='red')
par(mfrow=c(1,1))

















## Find kidneys that were predicted differently between PR0 and PR1 ##
d <-rep(0,n2)
for (i in 1:n2) {
  if (rpr0[i] == rpr1[i]) {
    d[i] <- 1
  }
  else {
    d[i] <- 0
  }
}

dd <- which(d == 0)


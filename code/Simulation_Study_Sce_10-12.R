########################## Simulations for Scenarios 10-12 ########################################

# written by Jeong Hoon Jang
# ver 4/12/2019


##### Call Functions #####
source("Nuclear_Function_05062018.R")
source("Simulation_10-12_Function.R")
source("Misc_Function.R")


##### Read in the data ######
renaldata <- read.delim("group1234_curve_expert_data.txt",header=T,sep="")

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



### Normalize Renogram curves and clinical variables ###
# Normalize curves (each baseline curve has maximum 25)
for(i in 1:n) {
  y2[,i] <- y2[,i]/max(y1[,i])*25
  y1[,i] <- y1[,i]/max(y1[,i])*25
}

X[1,] <- (X[1,] - mean(X[1,]))/sd(X[1,]) # Normalize age



## Set hyperparameters for extracting posterior means ##
nsample <- 10000
burn <- 2000
t1 <- as.numeric(Baseline_Time_Interval)/max(as.numeric(Baseline_Time_Interval))
t2 <- as.numeric(Diuretic_Time_Interval)/max(as.numeric(Diuretic_Time_Interval))
p <- 15
q <- c(3,3)
nu <- 3 
psi <- 1/p
ups <- 1
a <-  c(1.5,1.5)
a_sigma <- 1
b_sigma <- 1
tau <- 5
kappa <- 5


sum <- apply(W,2,sum)
mj <- rep(0, n)
for (i in 1:n) {
  if (sum[i] < 4) {
    mj[i] = 0
  }
  else {
    mj[i] = 1
  }
}
ep <- sum(mj)/n
theta <- log(-(ep-1)/ep)


## Extract Posterior Means for renogram-curve related parameters (from the full dataset above) ##
set.seed(227920)
Result <- Nuclear(nsample,y1,y2,X,W,t1,t2,p,q,nu,psi,ups,a,a_sigma,b_sigma,tau,kappa,theta)

# C
result_pC <- Result$pC
rpC <- apply(result_pC[,-(1:burn)],1,mean)
rpCC <- rep(0,n1)
for (i in 1:n1) {
  if (rpC[i] < 0.5){
    rpCC[i] = 0
  } else {
    rpCC[i]=1
  }
}
prob <- length(rpCC[which(rpCC == "1")])/n1


# delta
result_delta <- Result$delta
delta0 <- apply(result_delta[,1,-(1:burn)],1,mean)
delta1 <- apply(result_delta[,2,-(1:burn)],1,mean)


# S
result_S <- Result$S
S0 <- apply(result_S[,,1,-(1:burn)],c(1,2),mean)
S1 <- apply(result_S[,,2,-(1:burn)],c(1,2),mean)


# gamma
result_gamma <- Result$gamma
gamma1 <- apply(result_gamma[,,1,-(1:burn)],c(1,2),mean)
gamma2 <- apply(result_gamma[,,2,-(1:burn)],c(1,2),mean)


#lambda
result_lambda <- Result$lam
lambda01 <- apply(result_lambda[,,1,1,-(1:burn)],c(1,2),mean)
lambda02 <- apply(result_lambda[,,1,2,-(1:burn)],c(1,2),mean)
lambda11 <- apply(result_lambda[,,2,1,-(1:burn)],c(1,2),mean)
lambda12 <- apply(result_lambda[,,2,2,-(1:burn)],c(1,2),mean)

#Sigma
result_Sigma <- Result$sig2
Sigma01 <- apply(result_Sigma[,1,1,-(1:burn)],1,mean)
Sigma02 <- apply(result_Sigma[,1,2,-(1:burn)],1,mean)
Sigma11 <- apply(result_Sigma[,2,1,-(1:burn)],1,mean)
Sigma12 <- apply(result_Sigma[,2,2,-(1:burn)],1,mean)

# Time points
t1 <- as.numeric(Baseline_Time_Interval)/max(as.numeric(Baseline_Time_Interval))
t2 <- as.numeric(Diuretic_Time_Interval)/max(as.numeric(Diuretic_Time_Interval))

#sigma
result_sigma <- Result$sigma2[-(1:burn)]
sigma <- mean(result_sigma[-(1:burn)])



#xi: Cutoff point for expert rating
result_xi <- Result$xi
xi <- apply(result_xi[,-(1:burn)],1,mean)








## Extract posterior means or set expert rating related parameters that differ by scenario 

# Poor Agreement among experts (scenario 10)
result_mu_rho <- Result$mu_rho
mu_rho0 <- apply(result_mu_rho[,1,-(1:burn)],1,mean)
mu_rho1 <- apply(result_mu_rho[,2,-(1:burn)],1,mean)
sig2_rho0 <- c(1.5,1.5)
sig2_rho1 <- c(1.5,1.5)



# Moderate Agreement among experts (scenario 11)
#result_mu_rho <- Result$mu_rho
#mu_rho0 <- apply(result_mu_rho[,1,-(1:burn)],1,mean)
#mu_rho1 <- c(0.5,1.7) 
#sig2_rho0 <- c(0.2,0.2)
#sig2_rho1 <- c(0.2,0.2)


# Strong Agreement among experts (scenario 12)
#mu_rho0 <- c(-0.5,-5)
#mu_rho1 <- c(0.5,5) 
#sig2_rho0 <- c(0.01,0.01)
#sig2_rho1 <- c(0.01,0.01)









## Hyperparamter settings for our MCMC algorithm
n1 <- 80 # of training set
n2 <- 40 # of testing set
m <- 100 # Number of monte carlo sets
nsample <- 10000 # Number of MCMC samples
burn <- 2000 # burn-in
nb <- 15 # of basis functions used as hyperparameter
nl <- c(3,3) # of latent factors used as hyperparameter
nl_d <- c(6,6) # of latent factors used as hyperparameter



## Conduct simulation studies ##
simresult <- simulation10(prob,delta0,delta1,S0,S1,gamma1,gamma2,lambda01,lambda11,lambda02,lambda12,
                      Sigma01,Sigma02,Sigma11,Sigma12,t1,t2,sigma,mu_rho0,mu_rho1,sig2_rho0,sig2_rho1,xi,
                      nsample,burn,nb,nl,nl_d,m,n1,n2)




## Store the results ##

kap <- simresult$kap
er0 <- simresult$er0
er1 <- simresult$er1

tr3 <- simresult$tr3
pr0 <- simresult$pr0
pr1 <- simresult$pr1
tr3_d <- simresult$tr3_d
pr0_d <- simresult$pr0_d
pr1_d <- simresult$pr1_d
mv <- simresult$mv
kc <- simresult$kc
kct <- simresult$kct

C_training <- simresult$C_training
C_testing <- simresult$C_testing



result1 <- matrix(rep(0,n1*m),nrow=n1)
mv_cc <- rep(0,m)
mv_sensitivity <- rep(0,m)
mv_specificity <- rep(0,m)
mv_ppv <- rep(0,m)
mv_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n1) {
    if (C_training[i,c] == 0 & mv[i,c] == 0) {
      result1[i,c] <- 1
    } else if (C_training[i,c] == 0 & mv[i,c] == 1) {
      result1[i,c] <- 2
    } else if (C_training[i,c] == 0 & mv[i,c] == 2) {
      result1[i,c] <- 3
    } else if (C_training[i,c] == 1 & mv[i,c] == 0) {
      result1[i,c] <- 4
    } else if (C_training[i,c] == 1 & mv[i,c] == 1) {
      result1[i,c] <- 5
    } else {
      result1[i,c] <- 6
    }
  }
  mv_cc[c] <- (sum(result1[,c]==1)+sum(result1[,c]==6))/n1
  mv_sensitivity[c] <- sum(result1[,c]==6)/(sum(result1[,c]==4)+sum(result1[,c]==5)+sum(result1[,c]==6))
  mv_specificity[c] <- sum(result1[,c]==1)/(sum(result1[,c]==1)+sum(result1[,c]==2)+sum(result1[,c]==3))
  mv_ppv[c] <- sum(result1[,c]==6)/(sum(result1[,c]==3)+sum(result1[,c]==6))
  mv_npv[c] <- sum(result1[,c]==1)/(sum(result1[,c]==1)+sum(result1[,c]==4))
}



result2 <- matrix(rep(0,n1*m),nrow=n1)
kc_cc <- rep(0,m)
kc_sensitivity <- rep(0,m)
kc_specificity <- rep(0,m)
kc_ppv <- rep(0,m)
kc_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n1) {
    if (simresult$C_training[i,c] == 1 & simresult$kc[i,c] == 1) {
      result2[i,c] <- 1
    } else if (simresult$C_training[i,c] == 1 & simresult$kc[i,c] == 0) {
      result2[i,c] <- 2
    } else if (simresult$C_training[i,c] == 0 & simresult$kc[i,c] == 1) {
      result2[i,c] <- 3
    } else {
      result2[i,c] <- 4
    }
  }
  kc_cc[c] <- (sum(result2[,c]==1)+sum(result2[,c]==4))/n1
  kc_sensitivity[c] <- sum(result2[,c]==1)/(sum(result2[,c]==1)+sum(result2[,c]==2))
  kc_specificity[c] <- sum(result2[,c]==4)/(sum(result2[,c]==4)+sum(result2[,c]==3))
  kc_ppv[c] <- sum(result2[,c]==1)/(sum(result2[,c]==1)+sum(result2[,c]==3))
  kc_npv[c] <- sum(result2[,c]==4)/(sum(result2[,c]==4)+sum(result2[,c]==2))
}



result3 <- matrix(rep(0,n2*m),nrow=n2)
kct_cc <- rep(0,m)
kct_sensitivity <- rep(0,m)
kct_specificity <- rep(0,m)
kct_ppv <- rep(0,m)
kct_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$kct[i,c] == 1) {
      result3[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$kct[i,c] == 0) {
      result3[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$kct[i,c] == 1) {
      result3[i,c] <- 3
    } else {
      result3[i,c] <- 4
    }
  }
  kct_cc[c] <- (sum(result3[,c]==1)+sum(result3[,c]==4))/n2
  kct_sensitivity[c] <- sum(result3[,c]==1)/(sum(result3[,c]==1)+sum(result3[,c]==2))
  kct_specificity[c] <- sum(result3[,c]==4)/(sum(result3[,c]==4)+sum(result3[,c]==3))
  kct_ppv[c] <- sum(result3[,c]==1)/(sum(result3[,c]==1)+sum(result3[,c]==3))
  kct_npv[c] <- sum(result3[,c]==4)/(sum(result3[,c]==4)+sum(result3[,c]==2))
}





result4 <- matrix(rep(0,n1*m),nrow=n1)
tr3_cc <- rep(0,m)
tr3_sensitivity <- rep(0,m)
tr3_specificity <- rep(0,m)
tr3_ppv <- rep(0,m)
tr3_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n1) {
    if (simresult$C_training[i,c] == 1 & simresult$tr3[i,c] == 1) {
      result4[i,c] <- 1
    } else if (simresult$C_training[i,c] == 1 & simresult$tr3[i,c] == 0) {
      result4[i,c] <- 2
    } else if (simresult$C_training[i,c] == 0 & simresult$tr3[i,c] == 1) {
      result4[i,c] <- 3
    } else {
      result4[i,c] <- 4
    }
  }
  tr3_cc[c] <- (sum(result4[,c]==1)+sum(result4[,c]==4))/n1
  tr3_sensitivity[c] <- sum(result4[,c]==1)/(sum(result4[,c]==1)+sum(result4[,c]==2))
  tr3_specificity[c] <- sum(result4[,c]==4)/(sum(result4[,c]==4)+sum(result4[,c]==3))
  tr3_ppv[c] <- sum(result4[,c]==1)/(sum(result4[,c]==1)+sum(result4[,c]==3))
  tr3_npv[c] <- sum(result4[,c]==4)/(sum(result4[,c]==4)+sum(result4[,c]==2))
}





result5 <- matrix(rep(0,n2*m),nrow=n2)
pr0_cc <- rep(0,m)
pr0_sensitivity <- rep(0,m)
pr0_specificity <- rep(0,m)
pr0_ppv <- rep(0,m)
pr0_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$pr0[i,c] == 1) {
      result5[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$pr0[i,c] == 0) {
      result5[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$pr0[i,c] == 1) {
      result5[i,c] <- 3
    } else {
      result5[i,c] <- 4
    }
  }
  pr0_cc[c] <- (sum(result5[,c]==1)+sum(result5[,c]==4))/n2
  pr0_sensitivity[c] <- sum(result5[,c]==1)/(sum(result5[,c]==1)+sum(result5[,c]==2))
  pr0_specificity[c] <- sum(result5[,c]==4)/(sum(result5[,c]==4)+sum(result5[,c]==3))
  pr0_ppv[c] <- sum(result5[,c]==1)/(sum(result5[,c]==1)+sum(result5[,c]==3))
  pr0_npv[c] <- sum(result5[,c]==4)/(sum(result5[,c]==4)+sum(result5[,c]==2))
}




result6 <- matrix(rep(0,n2*m),nrow=n2)
pr1_cc <- rep(0,m)
pr1_sensitivity <- rep(0,m)
pr1_specificity <- rep(0,m)
pr1_ppv <- rep(0,m)
pr1_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$pr1[i,c] == 1) {
      result6[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$pr1[i,c] == 0) {
      result6[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$pr1[i,c] == 1) {
      result6[i,c] <- 3
    } else {
      result6[i,c] <- 4
    }
  }
  pr1_cc[c] <- (sum(result6[,c]==1)+sum(result6[,c]==4))/n2
  pr1_sensitivity[c] <- sum(result6[,c]==1)/(sum(result6[,c]==1)+sum(result6[,c]==2))
  pr1_specificity[c] <- sum(result6[,c]==4)/(sum(result6[,c]==4)+sum(result6[,c]==3))
  pr1_ppv[c] <- sum(result6[,c]==1)/(sum(result6[,c]==1)+sum(result6[,c]==3))
  pr1_npv[c] <- sum(result6[,c]==4)/(sum(result6[,c]==4)+sum(result6[,c]==2))
}





result7 <- matrix(rep(0,n1*m),nrow=n1)
tr3_d_cc <- rep(0,m)
tr3_d_sensitivity <- rep(0,m)
tr3_d_specificity <- rep(0,m)
tr3_d_ppv <- rep(0,m)
tr3_d_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n1) {
    if (simresult$C_training[i,c] == 1 & simresult$tr3_d[i,c] == 1) {
      result7[i,c] <- 1
    } else if (simresult$C_training[i,c] == 1 & simresult$tr3_d[i,c] == 0) {
      result7[i,c] <- 2
    } else if (simresult$C_training[i,c] == 0 & simresult$tr3_d[i,c] == 1) {
      result7[i,c] <- 3
    } else {
      result7[i,c] <- 4
    }
  }
  tr3_d_cc[c] <- (sum(result7[,c]==1)+sum(result7[,c]==4))/n1
  tr3_d_sensitivity[c] <- sum(result7[,c]==1)/(sum(result7[,c]==1)+sum(result7[,c]==2))
  tr3_d_specificity[c] <- sum(result7[,c]==4)/(sum(result7[,c]==4)+sum(result7[,c]==3))
  tr3_d_ppv[c] <- sum(result7[,c]==1)/(sum(result7[,c]==1)+sum(result7[,c]==3))
  tr3_d_npv[c] <- sum(result7[,c]==4)/(sum(result7[,c]==4)+sum(result7[,c]==2))
}






result8 <- matrix(rep(0,n2*m),nrow=n2)
pr0_d_cc <- rep(0,m)
pr0_d_sensitivity <- rep(0,m)
pr0_d_specificity <- rep(0,m)
pr0_d_ppv <- rep(0,m)
pr0_d_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$pr0_d[i,c] == 1) {
      result8[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$pr0_d[i,c] == 0) {
      result8[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$pr0_d[i,c] == 1) {
      result8[i,c] <- 3
    } else {
      result8[i,c] <- 4
    }
  }
  pr0_d_cc[c] <- (sum(result8[,c]==1)+sum(result8[,c]==4))/n2
  pr0_d_sensitivity[c] <- sum(result8[,c]==1)/(sum(result8[,c]==1)+sum(result8[,c]==2))
  pr0_d_specificity[c] <- sum(result8[,c]==4)/(sum(result8[,c]==4)+sum(result8[,c]==3))
  pr0_d_ppv[c] <- sum(result8[,c]==1)/(sum(result8[,c]==1)+sum(result8[,c]==3))
  pr0_d_npv[c] <- sum(result8[,c]==4)/(sum(result8[,c]==4)+sum(result8[,c]==2))
}




result9 <- matrix(rep(0,n2*m),nrow=n2)
pr1_d_cc <- rep(0,m)
pr1_d_sensitivity <- rep(0,m)
pr1_d_specificity <- rep(0,m)
pr1_d_ppv <- rep(0,m)
pr1_d_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$pr1_d[i,c] == 1) {
      result9[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$pr1_d[i,c] == 0) {
      result9[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$pr1_d[i,c] == 1) {
      result9[i,c] <- 3
    } else {
      result9[i,c] <- 4
    }
  }
  pr1_d_cc[c] <- (sum(result9[,c]==1)+sum(result9[,c]==4))/n2
  pr1_d_sensitivity[c] <- sum(result9[,c]==1)/(sum(result9[,c]==1)+sum(result9[,c]==2))
  pr1_d_specificity[c] <- sum(result9[,c]==4)/(sum(result9[,c]==4)+sum(result9[,c]==3))
  pr1_d_ppv[c] <- sum(result9[,c]==1)/(sum(result9[,c]==1)+sum(result9[,c]==3))
  pr1_d_npv[c] <- sum(result9[,c]==4)/(sum(result9[,c]==4)+sum(result9[,c]==2))
}








result10 <- matrix(rep(0,n1*m),nrow=n1)
fc_cc <- rep(0,m)
fc_sensitivity <- rep(0,m)
fc_specificity <- rep(0,m)
fc_ppv <- rep(0,m)
fc_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n1) {
    if (simresult$C_training[i,c] == 1 & simresult$fc[i,c] == 1) {
      result10[i,c] <- 1
    } else if (simresult$C_training[i,c] == 1 & simresult$fc[i,c] == 0) {
      result10[i,c] <- 2
    } else if (simresult$C_training[i,c] == 0 & simresult$fc[i,c] == 1) {
      result10[i,c] <- 3
    } else {
      result10[i,c] <- 4
    }
  }
  fc_cc[c] <- (sum(result10[,c]==1)+sum(result10[,c]==4))/n1
  fc_sensitivity[c] <- sum(result10[,c]==1)/(sum(result10[,c]==1)+sum(result10[,c]==2))
  fc_specificity[c] <- sum(result10[,c]==4)/(sum(result10[,c]==4)+sum(result10[,c]==3))
  fc_ppv[c] <- sum(result10[,c]==1)/(sum(result10[,c]==1)+sum(result10[,c]==3))
  fc_npv[c] <- sum(result10[,c]==4)/(sum(result10[,c]==4)+sum(result10[,c]==2))
}



result11 <- matrix(rep(0,n2*m),nrow=n2)
fct_cc <- rep(0,m)
fct_sensitivity <- rep(0,m)
fct_specificity <- rep(0,m)
fct_ppv <- rep(0,m)
fct_npv <- rep(0,m)
for (c in 1:m) {
  for(i in 1:n2) {
    if (simresult$C_testing[i,c] == 1 & simresult$fct[i,c] == 1) {
      result11[i,c] <- 1
    } else if (simresult$C_testing[i,c] == 1 & simresult$fct[i,c] == 0) {
      result11[i,c] <- 2
    } else if (simresult$C_testing[i,c] == 0 & simresult$fct[i,c] == 1) {
      result11[i,c] <- 3
    } else {
      result11[i,c] <- 4
    }
  }
  fct_cc[c] <- (sum(result11[,c]==1)+sum(result11[,c]==4))/n2
  fct_sensitivity[c] <- sum(result11[,c]==1)/(sum(result11[,c]==1)+sum(result11[,c]==2))
  fct_specificity[c] <- sum(result11[,c]==4)/(sum(result11[,c]==4)+sum(result11[,c]==3))
  fct_ppv[c] <- sum(result11[,c]==1)/(sum(result11[,c]==1)+sum(result11[,c]==3))
  fct_npv[c] <- sum(result11[,c]==4)/(sum(result11[,c]==4)+sum(result11[,c]==2))
}














####### Simulation results ########
er0 <- mean(er0) # Average Expert Rating for non-obstructed kidneys
er1 <- mean(er1) # Average Expert Rating for obstructed kidneys
kap <- mean(simresult$kap) # Kappa agreement measure between the three experts



# TR3 - True number of latent factors
tr3cc <- mean(tr3_cc,na.rm = TRUE)
tr3cc_sd <- sd(tr3_cc,na.rm = TRUE)/sqrt(m)

tr3sensitivity <- mean(tr3_sensitivity,na.rm = TRUE)
tr3sensitivity_sd <- sd(tr3_sensitivity,na.rm = TRUE)/sqrt(m)


tr3specificity <- mean(tr3_specificity,na.rm = TRUE)
tr3specificity_sd <- sd(tr3_specificity,na.rm = TRUE)/sqrt(m)

tr3ppv <- mean(tr3_ppv,na.rm = TRUE)
tr3ppv_sd <- sd(tr3_ppv,na.rm = TRUE)/sqrt(m)

tr3npv <- mean(tr3_npv,na.rm = TRUE)
tr3npv_sd <- sd(tr3_npv,na.rm = TRUE)/sqrt(m)



#PR0 - True number of latent factors
pr0cc <- mean(pr0_cc,na.rm = TRUE)
pr0cc_sd <- sd(pr0_cc,na.rm = TRUE)/sqrt(m)

pr0sensitivity <- mean(pr0_sensitivity,na.rm = TRUE)
pr0sensitivity_sd <- sd(pr0_sensitivity,na.rm = TRUE)/sqrt(m)


pr0specificity <- mean(pr0_specificity,na.rm = TRUE)
pr0specificity_sd <- sd(pr0_specificity,na.rm = TRUE)/sqrt(m)

pr0ppv <- mean(pr0_ppv,na.rm = TRUE)
pr0ppv_sd <- sd(pr0_ppv,na.rm = TRUE)/sqrt(m)

pr0npv <- mean(pr0_npv,na.rm = TRUE)
pr0npv_sd <- sd(pr0_npv,na.rm = TRUE)/sqrt(m)




#PR1 - True number of latent factors
pr1cc <- mean(pr1_cc,na.rm = TRUE)
pr1cc_sd <- sd(pr1_cc,na.rm = TRUE)/sqrt(m)

pr1sensitivity <- mean(pr1_sensitivity,na.rm = TRUE)
pr1sensitivity_sd <- sd(pr1_sensitivity,na.rm = TRUE)/sqrt(m)


pr1specificity <- mean(pr1_specificity,na.rm = TRUE)
pr1specificity_sd <- sd(pr1_specificity,na.rm = TRUE)/sqrt(m)

pr1ppv <- mean(pr1_ppv,na.rm = TRUE)
pr1ppv_sd <- sd(pr1_ppv,na.rm = TRUE)/sqrt(m)

pr1npv <- mean(pr1_npv,na.rm = TRUE)
pr1npv_sd <- sd(pr1_npv,na.rm = TRUE)/sqrt(m)





# TR3 - Wrong number of latent factors
tr3_dcc <- mean(tr3_d_cc,na.rm = TRUE)
tr3_dcc_sd <- sd(tr3_d_cc,na.rm = TRUE)/sqrt(m)

tr3_dsensitivity <- mean(tr3_d_sensitivity,na.rm = TRUE)
tr3_dsensitivity_sd <- sd(tr3_d_sensitivity,na.rm = TRUE)/sqrt(m)


tr3_dspecificity <- mean(tr3_d_specificity,na.rm = TRUE)
tr3_dspecificity_sd <- sd(tr3_d_specificity,na.rm = TRUE)/sqrt(m)

tr3_dppv <- mean(tr3_d_ppv,na.rm = TRUE)
tr3_dppv_sd <- sd(tr3_d_ppv,na.rm = TRUE)/sqrt(m)

tr3_dnpv <- mean(tr3_d_npv,na.rm = TRUE)
tr3_dnpv_sd <- sd(tr3_d_npv,na.rm = TRUE)/sqrt(m)





# PR0 - Wrong number of latent factors
pr0_dcc <- mean(pr0_d_cc,na.rm = TRUE)
pr0_dcc_sd <- sd(pr0_d_cc,na.rm = TRUE)/sqrt(m)

pr0_dsensitivity <- mean(pr0_d_sensitivity,na.rm = TRUE)
pr0_dsensitivity_sd <- sd(pr0_d_sensitivity,na.rm = TRUE)/sqrt(m)


pr0_dspecificity <- mean(pr0_d_specificity,na.rm = TRUE)
pr0_dspecificity_sd <- sd(pr0_d_specificity,na.rm = TRUE)/sqrt(m)

pr0_dppv <- mean(pr0_d_ppv,na.rm = TRUE)
pr0_dppv_sd <- sd(pr0_d_ppv,na.rm = TRUE)/sqrt(m)

pr0_dnpv <- mean(pr0_d_npv,na.rm = TRUE)
pr0_dnpv_sd <- sd(pr0_d_npv,na.rm = TRUE)/sqrt(m)





# PR1 - Wrong number of latent factors
pr1_dcc <- mean(pr1_d_cc,na.rm = TRUE)
pr1_dcc_sd <- sd(pr1_d_cc,na.rm = TRUE)/sqrt(m)

pr1_dsensitivity <- mean(pr1_d_sensitivity,na.rm = TRUE)
pr1_dsensitivity_sd <- sd(pr1_d_sensitivity,na.rm = TRUE)/sqrt(m)

pr1_dspecificity <- mean(pr1_d_specificity,na.rm = TRUE)
pr1_dspecificity_sd <- sd(pr1_d_specificity,na.rm = TRUE)/sqrt(m)

pr1_dppv <- mean(pr1_d_ppv,na.rm = TRUE)
pr1_dppv_sd <- sd(pr1_d_ppv,na.rm = TRUE)/sqrt(m)

pr1_dnpv <- mean(pr1_d_npv,na.rm = TRUE)
pr1_dnpv_sd <- sd(pr1_d_npv,na.rm = TRUE)/sqrt(m)






# Majority Voting
mvcc <- mean(mv_cc,na.rm = TRUE)
mvcc_sd <- sd(mv_cc,na.rm = TRUE)/sqrt(m)

mvsensitivity <- mean(mv_sensitivity,na.rm = TRUE)
mvsensitivity_sd <- sd(mv_sensitivity,na.rm = TRUE)/sqrt(m)

mvspecificity <- mean(mv_specificity,na.rm = TRUE)
mvspecificity_sd <- sd(mv_specificity,na.rm = TRUE)/sqrt(m)

mvppv <- mean(mv_ppv,na.rm = TRUE)
mvppv_sd <- sd(mv_ppv,na.rm = TRUE)/sqrt(m)

mvnpv <- mean(mv_npv,na.rm = TRUE)
mvnpv_sd <- sd(mv_npv,na.rm = TRUE)/sqrt(m)






# Result of K-mean clustering of curves and clinical variables on training data
kccc <- mean(kc_cc,na.rm = TRUE)
kccc_sd <- sd(kc_cc,na.rm = TRUE)/sqrt(m)

kcsensitivity <- mean(kc_sensitivity,na.rm = TRUE)
kcsensitivity_sd <- sd(kc_sensitivity,na.rm = TRUE)/sqrt(m)

kcspecificity <- mean(kc_specificity,na.rm = TRUE)
kcspecificity_sd <- sd(kc_specificity,na.rm = TRUE)/sqrt(m)

kcppv <- mean(kc_ppv,na.rm = TRUE)
kcppv_sd <- sd(kc_ppv,na.rm = TRUE)/sqrt(m)

kcnpv <- mean(kc_npv,na.rm = TRUE)
kcnpv_sd <- sd(kc_npv,na.rm = TRUE)/sqrt(m)





#  Prediction of K-mean clustering of curves and clinical variables on testing data
kctcc <- mean(kct_cc,na.rm = TRUE)
kctcc_sd <- sd(kct_cc,na.rm = TRUE)/sqrt(m)

kctsensitivity <- mean(kct_sensitivity,na.rm = TRUE)
kctsensitivity_sd <- sd(kct_sensitivity,na.rm = TRUE)/sqrt(m)

kctspecificity <- mean(kct_specificity,na.rm = TRUE)
kctspecificity_sd <- sd(kct_specificity,na.rm = TRUE)/sqrt(m)

kctppv <- mean(kct_ppv,na.rm = TRUE)
kctppv_sd <- sd(kct_ppv,na.rm = TRUE)/sqrt(m)

kctnpv <- mean(kct_npv,na.rm = TRUE)
kctnpv_sd <- sd(kct_npv,na.rm = TRUE)/sqrt(m)




# Result of multivariate functional clustering of curves on training data
fccc <- mean(fc_cc,na.rm = TRUE)
fccc_sd <- sd(fc_cc,na.rm = TRUE)/sqrt(m)

fcsensitivity <- mean(fc_sensitivity,na.rm = TRUE)
fcsensitivity_sd <- sd(fc_sensitivity,na.rm = TRUE)/sqrt(m)

fcspecificity <- mean(fc_specificity,na.rm = TRUE)
fcspecificity_sd <- sd(fc_specificity,na.rm = TRUE)/sqrt(m)

fcppv <- mean(fc_ppv,na.rm = TRUE)
fcppv_sd <- sd(fc_ppv,na.rm = TRUE)/sqrt(m)

fcnpv <- mean(fc_npv,na.rm = TRUE)
fcnpv_sd <- sd(fc_npv,na.rm = TRUE)/sqrt(m)





#  Prediction of multivariate functional clustering of curves on testing data
fctcc <- mean(fct_cc,na.rm = TRUE)
fctcc_sd <- sd(fct_cc,na.rm = TRUE)/sqrt(m)

fctsensitivity <- mean(fct_sensitivity,na.rm = TRUE)
fctsensitivity_sd <- sd(fct_sensitivity,na.rm = TRUE)/sqrt(m)

fctspecificity <- mean(fct_specificity,na.rm = TRUE)
fctspecificity_sd <- sd(fct_specificity,na.rm = TRUE)/sqrt(m)

fctppv <- mean(fct_ppv,na.rm = TRUE)
fctppv_sd <- sd(fct_ppv,na.rm = TRUE)/sqrt(m)

fctnpv <- mean(fct_npv,na.rm = TRUE)
fctnpv_sd <- sd(fct_npv,na.rm = TRUE)/sqrt(m)
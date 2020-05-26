########################## Simulations for Scenarios 1-9 ########################################
# written by Jeong Hoon Jang
# ver 20190304


##### Call Functions #####
source("Nuclear_Function_05062018.R")
source("Simulation_1-9_Function.R")
source("Misc_Function.R")










# For small discrepancy between renogram curves
lambda01 <- matrix(c(3,1.5,0.5,-1.5,-3,-2,-1,0.5,1,2,1,0.7,0.5,-0.7,-1),nrow=5)
lambda11 <- matrix(c(3,1.5,0.5,-1.5,-3,-2,-1,0.5,1,2,1,0.7,0.5,-0.7,-1),nrow=5)
lambda02 <- matrix(c(2,1,0.5,-1,-2,-1,-0.7,0.5,0.7,1,3,1.5,0.5,-1.5,-3),nrow=5)
lambda12 <- matrix(c(2,1,0.5,-1,-2,-1,-0.7,0.5,0.7,1,3,1.5,0.5,-1.5,-3),nrow=5)


# For moderate discrepancy between renogram curves
#lambda01 <- matrix(c(3,1.5,0.5,-1.5,-3,-2,-1,0.5,1,2,1,0.7,0.5,-0.7,-1),nrow=5)
#lambda11 <- matrix(c(3.5,2,1,-1.5,-3,-2,-1,1,1.5,2.5,1.5,1.2,1,-0.7,-1),nrow=5)
#lambda02 <- matrix(c(2,1,0.5,-1,-2,-1,-0.7,0.5,0.7,1,3,1.5,0.5,-1.5,-3),nrow=5)
#lambda12 <- matrix(c(2.5,1.5,1,-1,-2,-1,-0.7,1,1.2,1.5,3.5,2,1,-1.5,-3),nrow=5)


# For large discrepancy between renogram curves
#lambda01 <- matrix(c(3,1.5,0.5,-1.5,-3,-2,-1,0.5,1,2,1,0.7,0.5,-0.7,-1),nrow=5)
#lambda11 <- matrix(c(3,3,1,-1.5,-3,-2,-1,1,2,4,2,1.4,1,-0.7,-1),nrow=5)
#lambda02 <- matrix(c(2,1,0.5,-1,-2,-1,-0.7,0.5,0.7,1,3,1.5,0.5,-1.5,-3),nrow=5)
#lambda12 <- matrix(c(2,2,1,-1,-2,-1,-0.7,1,1.4,2,6,3,1,-1.5,-3),nrow=5)



# For poor agreement among expert ratings
sigma_rho0g <- 1
sigma_rho1g <- 1

# For moderate agreement among expert ratings
#sigma_rho0g <- 0.4
#sigma_rho1g <- 0.4

# For strong agreement among expert ratings
#sigma_rho0g <- 0.01
#sigma_rho1g <- 0.01






# Specify scenario conditions and hyperparmeters
n1 <- 80 # of training set
n2 <- 40 # of testing set
m <- 100 # Number of monte carlo sets
nsample <- 10000 # Number of MCMC samples
burn <- 2000 # burn-in
nb <- 15 # of basis functions used as hyperparameter
nl <- c(3,3) # of latent factors used as hyperparameter
nl_d <- c(6,6) # of latent factors used as hyperparameter





# Conduct simulation
simresult <- simulation(lambda01,lambda11,lambda02,lambda12,sigma_rho0g,sigma_rho1g,
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
fc <- simresult$fc
fct <- simresult$fct

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
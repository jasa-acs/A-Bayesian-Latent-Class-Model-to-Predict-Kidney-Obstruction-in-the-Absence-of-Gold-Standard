##### Function for Conducting Simulation Study Under Scenarios 10-12 #########
# written by Jeong Hoon Jang
# ver 04/11/2019
# prob: Extracted Prevalence of kidney obstruction used to generate underlying obstruction status for simulation
# delta0, delta 1: Extracted Posterior Mean of mean vector of clinical variables used to generate clinical variables for simulation
# S0, S1: Extracted posterior mean of covariance matrix clinical variables used to generate data clinical variables for simulation
# gamma1, gamma2: Extracted posterior mean of gammas
# lambda01,lambda11,lambda02,lambda12: Extracted Posterior mean of latent factor loadings used to generate renogram curves for simulation
# Sigma01,Sigma02,Sigma11,Sigma12: Extracted Posterior mean of Diagnoal Covariance of beta coefficients used to generate renogram curves for simulation
# t1,t2: Observed Normalized Time points in [0,1] where renogram curves will be generated on
# sigma: Extracted Posterior mean of measurement error variance in renogram curves
# mu_rho0,mu_rho1: Mean Parameters used to generate expert ratings (can be posterior means or artifical values depending on the scenario)
# sig2_rho0,sig2_rho1: Variance Parameters used to generate expert ratings (can be posterior means or artifical values depending on the scenario)
# xi: Extracted Posterior mean of Cutoff points used to generate expert ratings for simulation
# nsample: Number of MCMC iterations of our algorithm applied to generated simulation data
# burn: Number of burn-ins used for our MCMC algorithm
# nb: number of basis functions that will be used for implementing our MCMC algorithm
# nl: True number of latent factors that will be used for implementing our MCMC algorithm
# nl_d:  Arbitual number of latent factors that will be used for implementing our MCMC algorithm
# m: Number of monte carlo datasets for evaluating simulation results
# n1, n2: number of training and testing dataset for each monte carlo dataset

source("Nuclear_Function_05062018.R")
source("Misc_Function.R")



if ( !require(LaplacesDemon) )
{
  install.packages("LaplacesDemon")
  library(LaplacesDemon)
}

if ( !require(MASS) )
{
  install.packages("MASS")
  library(MASS)
}

if ( !require(irr) )
{
  install.packages("irr")
  library(irr)
}




if ( !require(funHDDC) )
{
  install.packages("funHDDC")
  library(funHDDC)
}






simulation10 <- function(prob,delta0,delta1,S0,S1,gamma1,gamma2,lambda01,lambda11,lambda02,lambda12,
                         Sigma01,Sigma02,Sigma11,Sigma12,t1,t2,sigma,mu_rho0,mu_rho1,sig2_rho0,sig2_rho1,xi,
                         nsample,burn,nb,nl,nl_d,m,n1,n2)
                           {
  kap <- rep(0,m)
  er0 <- rep(0,m)
  er1 <- rep(0,m)
  
  
  C_training <- matrix(rep(0,n1*m),nrow=n1)
  C_testing <- matrix(rep(0,n2*m),nrow=n2)
  
  
  tr3 <- matrix(rep(0,n1*m),nrow=n1)
  pr0 <- matrix(rep(0,n2*m),nrow=n2)
  pr1 <- matrix(rep(0,n2*m),nrow=n2)
  
  
  tr3_d <- matrix(rep(0,n1*m),nrow=n1)
  pr0_d <- matrix(rep(0,n2*m),nrow=n2)
  pr1_d <- matrix(rep(0,n2*m),nrow=n2)
  
  mv <- matrix(rep(0,n1*m),nrow=n1)
  
  kc <- matrix(rep(0,n1*m),nrow=n1)
  kct <- matrix(rep(0,n2*m),nrow=n2)
  
  fc <- matrix(rep(0,n1*m),nrow=n1)
  fct <- matrix(rep(0,n2*m),nrow=n2)
  
  n <- n1+n2
  
  for (c in 1:m){
    
    set.seed(c*200)
    ### Step 1 ###
    
    C <- rbern(n, prob) # Generation of true obstruction status for the sample (C=1 Obstructed; C=0 Non-Obstructed)
    
    
    
    
    ### Step 2 ###
    r <- 2   # Number of Covariates
    
    # Generate Covariates X and store it to r X n matrix
    X <- matrix(rep(0, n*r), nrow=r)
    for (i in 1:n) {
      if (C[i] < 1) {
        X[,i] <- mvrnorm(1, delta0, S0)
      }
      else if (C[i] > 0) {
        X[,i] <-  mvrnorm(1, delta1, S1)
      }
    }
    
    
    
    ### Step 3 ###
    q1 <- 3 # Number of latent factors for Baseline Scan (s=1)
    q2 <- 3 # Number of latent factors for Post-furosemide Scan (s=2)
    
    
    # Generate q(s) X n latent factors for both scans s=1 and s=2
    eta1 <- matrix(rep(0,n*q1),nrow=q1)
    for (i in 1:n){
      eta1[,i] <- mvrnorm(1, gamma1%*%X[,i], diag(1,nrow=q1,ncol=q1))
    }
    eta2 <- matrix(rep(0,n*q2),nrow=q2)
    for (i in 1:n){
      eta2[,i] <- mvrnorm(1, gamma2%*%X[,i], diag(1,nrow=q2,ncol=q2))
    }
    
    
    
    
    ### Step 4 ###
    p <- 15  # Number of Basis Functions
    
    
    # Generate n x p matrix of beta coefficients for both scans (s= 1,2)
    beta1 <- matrix(rep(0,n*p),nrow=n)
    beta2 <- matrix(rep(0,n*p),nrow=n)
    for (i in 1:n) {
      if (C[i] > 0){
        beta1[i,] <- mvrnorm(1,lambda11%*%eta1[,i], diag(Sigma11))
        beta2[i,] <- mvrnorm(1,lambda12%*%eta2[,i], diag(Sigma12))
      }
      if (C[i] < 1){
        beta1[i,] <- mvrnorm(1,lambda01%*%eta1[,i], diag(Sigma01))
        beta2[i,] <- mvrnorm(1,lambda02%*%eta2[,i], diag(Sigma02))
      }
    }
    
    
    
    ### Step 5 ###
    m1 <- length(t1) # Number of time points for baseline scan s=1
    m2 <- length(t2) # Number of time points for post-furosemide scan s=2
    nu <- 3 # bandwidth 
    psi <- 1/p
    
    # Generate basis functions
    b1 <- make_basis(p,nu,psi, t1,t2)$b1
    b2 <- make_basis(p,nu,psi, t1,t2)$b2
    
    
    
    # Generate n X m1 observarions for baseline scan s=1 and n X m2 observations for post-furosemide scan
    y1 <- matrix(rep(0,n*m1),nrow=m1)
    y2 <- matrix(rep(0,n*m2),nrow=m2)
    for (i in 1:n) {
      y1[,i] <- mvrnorm(1,b1%*%beta1[i,],diag(sigma,nrow=m1, ncol=m1))
      y2[,i] <- mvrnorm(1,b2%*%beta2[i,],diag(sigma,nrow=m2, ncol=m2))
    }
    
    ### Step 6 ###
    L <- 3 # Number of Experts
    
    
    Sigma_rho0 <- diag(sig2_rho0, r)
    Sigma_rho1 <- diag(sig2_rho1, r)
    
    rho0 <- mvrnorm(L,mu_rho0, Sigma_rho0) #L X r Coefficient for expert mean given c=0
    rho1 <- mvrnorm(L,mu_rho1, Sigma_rho1) #L X r Coefficient for expert mean given c=1
    
    
    # Generate n X L matrix of underlying continuous expert ratings
    Z <- matrix(rep(0,n*L),nrow=L)
    for (i in 1:n) {
      if (C[i] < 1) {
        Z[,i] <- mvrnorm(1,rho0%*%X[,i],diag(1,L))
      }
      if (C[i] > 0) {
        Z[,i] <- mvrnorm(1,rho1%*%X[,i],diag(1,L))
      }
    }
    
    
    
    ### Step 7 ###
    # Generate a n X L matrix of actual ordinal ratings from the experts
    W <- matrix(rep(0,n*L),nrow=L)
    for (i in 1:L) {
      for (j in 1:n) {
        if (Z[i,j] <= 0) {
          W[i,j] = 0
        }
        if (0 < Z[i,j] & Z[i,j] <= xi[i]) {
          W[i,j] = 1
        } 
        if (xi[i] < Z[i,j]) {
          W[i,j] = 2
        } 
      }
    }
    

    
    
    ## Dividing into testing and training dataset
    index <- 1:n
    train <- sort(sample(index,n1))
    test <- index[-train]
    
    
    y1_training <- y1[,train]
    y2_training <- y2[,train]
    y1_testing <- y1[,test]
    y2_testing <- y2[,test]
    
    X_training <- X[,train]
    X_testing <- X[,test]
    
    W_training <- W[,train]
    W_testing <- W[,test]
    
    C_training[,c] <- C[train]
    C_testing[,c] <- C[test]
    
    
    
    # Average and agreement for generated expert ratings
    av  <- apply(W_training,2,mean)
    avv <- cbind(C_training[,c],av)
    avv0 <- avv[which(C_training[,c]=="0"),]
    avv1 <- avv[which(C_training[,c]=="1"),]
    er0[c] <- mean(avv0[,2]) # Mean of three expert ratings in the training data given c=0
    er1[c] <- mean(avv1[,2]) # Mean of three expert ratings in the training data given c=1
    kap[c] <- kappam.fleiss(t(W_training))[[5]] # agreement of three expert ratings in the training data
    
    
    
    ### Our proposed algorithm
    
    #selecting theta based on expert rating
    es <- rep(0,n1)
    WW <- apply(W_training,2,sum)
    for (i in 1:n1) {
      if (WW[i] < 4) {
        es[i] = 0
      } else {
        es[i] = 1
      }
    }
    
    ep <- sum(es)/n1
    emptheta <- log(-(ep-1)/ep)
    
    
    
    ## True Number of Latent Factors
    # Training 
    Result <- Nuclear(nsample=nsample,y1=y1_training,y2=y2_training,X=X_training,W=W_training,t1=t1,t2=t2,
                      p=nb,q=nl,nu=3,psi=1/nb,ups=1,a=c(1.5,1.5),a_sigma=1,b_sigma=1,tau=5,kappa=5,theta=emptheta)
    
    train_pc <- Result$pC
    train_pcc <- apply(train_pc[,(burn+1):nsample],1,mean)
    for (i in 1:n1) {
      if (train_pcc[i] < 0.5){
        tr3[i,c] = 0
      } else {
        tr3[i,c] = 1
      }
    }
    
    
    # Testing
    set.seed(c*10)
    g <- sample(c(1,2,3), n2, replace=TRUE)
    test1_pc <- matrix(rep(0,n2*(nsample-burn)),nrow=n2)
    test2_pc <- matrix(rep(0,n2*(nsample-burn)),nrow=n2)
    for(s in 1:n2) {
      test1_pc[s,] <- Nuclear_predict(fit=Result,bi=burn,y1=y1_testing[,s],y2=y2_testing[,s],x=X_testing[,s],W=NULL)
      test2_pc[s,] <- Nuclear_predict(fit=Result,bi=burn,y1=y1_testing[,s],y2=y2_testing[,s],x=X_testing[,s],W=W_testing[g[s],s])
      print(s)
    }
    
    test1_pcc <- apply(test1_pc,1,mean)
    test2_pcc <- apply(test2_pc,1,mean)
    
    for (i in 1:n2) {
      if (test1_pcc[i] < 0.5){
        pr0[i,c] = 0
      } else {
        pr0[i,c]=1
      }
    }
    
    
    for (i in 1:n2) {
      if (test2_pcc[i] < 0.5){
        pr1[i,c] = 0
      } else {
        pr1[i,c]=1
      }
    }
    


    
    
    
    ## Different number of Latent Factors
    Result_d <- Nuclear(nsample=nsample,y1=y1_training,y2=y2_training,X=X_training,W=W_training,t1=t1,t2=t2,
                        p=nb,q=nl_d,nu=3,psi=1/nb,ups=1,a=c(1.5,1.5),a_sigma=1,b_sigma=1,tau=5,kappa=5,theta=emptheta)
    
    train_pc_d <- Result_d$pC
    train_pcc_d <- apply(train_pc_d[,(burn+1):nsample],1,mean)
    for (i in 1:n1) {
      if (train_pcc_d[i] < 0.5){
        tr3_d[i,c] = 0
      } else {
        tr3_d[i,c]=1
      }
    }
    


    
    # Testing
    test1_pc_d <- matrix(rep(0,n2*(nsample-burn)),nrow=n2)
    test2_pc_d <- matrix(rep(0,n2*(nsample-burn)),nrow=n2)
    for(s in 1:n2) {
      test1_pc_d[s,] <- Nuclear_predict(fit=Result_d,bi=burn,y1=y1_testing[,s],y2=y2_testing[,s],x=X_testing[,s],W=NULL)
      test2_pc_d[s,] <- Nuclear_predict(fit=Result_d,bi=burn,y1=y1_testing[,s],y2=y2_testing[,s],x=X_testing[,s],W=W_testing[g[s],s])
      print(s)
    }
    
    test1_pcc_d <- apply(test1_pc_d,1,mean)
    test2_pcc_d <- apply(test2_pc_d,1,mean)
    

    for (i in 1:n2) {
      if (test1_pcc_d[i] < 0.5){
        pr0_d[i,c] = 0
      } else {
        pr0_d[i,c]=1
      }
    }
    
 
    
  for (i in 1:n2) {
      if (test2_pcc_d[i] < 0.5){
        pr1_d[i,c] = 0
      } else {
        pr1_d[i,c]=1
      }
    }

    
    
    ### Naive Methods
    
    # Based on Majority Voting
    for (i in 1:n1) {
        mv[i,c] <- majority(W_training[,i]) 
    }
    
    
    
    
    
    # Based on K-mean clustering of C(Y,X)
    mdata <- cbind(t(y1_training),t(y2_training),t(X_training))
    
    set.seed(c*10)
    kmean <-   kmeans(mdata, 2, nstart=5)
    cc <- kmean[[1]]
    kcenter <- kmean[[2]]
    
    rtable_1 <- table(C_training[,c],cc)
    
    if (rtable_1[1,1]+rtable_1[2,2] >= rtable_1[1,2]+rtable_1[2,1]) {
      for (i in 1:n1) {
        if (cc[i]==1){
          kc[i,c] = 0
        } else {
          kc[i,c] =1
        }
        
      }
      
    } else if (rtable_1[1,1]+rtable_1[2,2] < rtable_1[1,2]+rtable_1[2,1]) {
      
      for (i in 1:n1) {
        if (cc[i]==1){
          kc[i,c] = 1
        } else {
          kc[i,c]=0
        }
      }
      A <- kcenter[1,] 
      B <- kcenter[2,]
      kcenter <- rbind(B,A)
    }
    
    
    for (i in 1:n2) {
      diff1 <- euc.dist(c(y1_testing[,i], y2_testing[,i], X_testing[,i]),kcenter[1,])
      diff2 <- euc.dist(c(y1_testing[,i], y2_testing[,i], X_testing[,i]),kcenter[2,])
      if (diff1 < diff2) {
        kct[i,c] = 0
      }
      else {
        kct[i,c] = 1
      }
    }
    
    
    
    
    
    
    
    
    
    ## Based on multivariate functional clustering
    basis1 <- create.bspline.basis(c(0,t1[m1]), nbasis=10)
    basiseval1 <- eval.basis(t1,basis1)
    
    basis2 <- create.bspline.basis(c(0,t2[m2]), nbasis=10)
    basiseval2 <- eval.basis(t2,basis2)
    
    
    
    
    #training 
    fy1_training <- smooth.basis(t1,y1_training,fdParobj=basis1)$fd
    fy2_training <- smooth.basis(t2,y2_training,fdParobj=basis2)$fd
    fy_training <- list(fy1_training,fy2_training)
    
    
    set.seed(c*10)
    fcluster <- funHDDC(fy_training,K=2)
    fclass <- fcluster$class-1
    
    fmean01 <- basiseval1%*%fcluster$mu[1,1:10]
    fmean11 <- basiseval1%*%fcluster$mu[2,1:10]
    
    fmean02 <- basiseval2%*%fcluster$mu[1,11:20]
    fmean12 <- basiseval2%*%fcluster$mu[2,11:20]
    
    
    ftable_1 <- table(C_training[,c],fclass)
    
    if (ftable_1[1,1]+ftable_1[2,2] >= ftable_1[1,2]+ftable_1[2,1]) {
      for (i in 1:n1) {
        if (fclass[i]==1){
          fc[i,c] = 1
        } else {
          fc[i,c] =0
        }
        
      }
      
    } else if (ftable_1[1,1]+ftable_1[2,2] < ftable_1[1,2]+ftable_1[2,1]) {
      
      for (i in 1:n1) {
        if (fclass[i]==1){
          fc[i,c] = 0
        } else {
          fc[i,c]=1
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
    coefs_y1_testing <- smooth.basis(t1,y1_testing,fdParobj=basis1)$fd$coefs
    smooth_y1_testing <- basiseval1%*%coefs_y1_testing
    coefs_y2_testing <- smooth.basis(t2,y2_testing,fdParobj=basis2)$fd$coefs
    smooth_y2_testing <- basiseval2%*%coefs_y2_testing
    
    for (i in 1:n2) {
      fdiff0 <- hilbert_dist(smooth_y1_testing[,i],fmean01,t1, smooth_y2_testing[,i],fmean02,t2)
      fdiff1 <- hilbert_dist(smooth_y1_testing[,i],fmean11,t1, smooth_y2_testing[,i],fmean12,t2)
      if (fdiff0 < fdiff1) {
        fct[i,c] = 0
      }
      else {
        fct[i,c] = 1
      }
    }  
    
    
   print(c)
    
  }
  list(kap=kap,er0=er0, er1=er1, C_training=C_training, C_testing=C_testing,
       tr3=tr3,pr0=pr0,pr1=pr1,tr3_d=tr3_d,pr0_d=pr0_d,pr1_d=pr1_d,mv=mv, 
       kc=kc, kct=kct,fc=fc,fct=fct)
}

























































temp <- function(prob,delta0,delta1,S0,S1,gamma1,gamma2,lambda01,lambda11,lambda02,lambda12,
                         Sigma01,Sigma02,Sigma11,Sigma12,t1,t2,sigma,mu_rho0,mu_rho1,sig2_rho0,sig2_rho1,xi,
                         m,n1,n2)
{
  kap <- rep(0,m)
  er0 <- rep(0,m)
  er1 <- rep(0,m)
  
  
  C_training <- matrix(rep(0,n1*m),nrow=n1)
  C_testing <- matrix(rep(0,n2*m),nrow=n2)
  
  

  
  mv <- matrix(rep(0,n1*m),nrow=n1)
  
  kc <- matrix(rep(0,n1*m),nrow=n1)
  kct <- matrix(rep(0,n2*m),nrow=n2)
  
  fc <- matrix(rep(0,n1*m),nrow=n1)
  fct <- matrix(rep(0,n2*m),nrow=n2)
  
  n <- n1+n2
  
  for (c in 1:m){
    
    set.seed(c*200)
    ### Step 1 ###
    
    C <- rbern(n, prob) # Generation of true obstruction status for the sample (C=1 Obstructed; C=0 Non-Obstructed)
    
    
    
    
    ### Step 2 ###
    r <- 2   # Number of Covariates
    
    # Generate Covariates X and store it to r X n matrix
    X <- matrix(rep(0, n*r), nrow=r)
    for (i in 1:n) {
      if (C[i] < 1) {
        X[,i] <- mvrnorm(1, delta0, S0)
      }
      else if (C[i] > 0) {
        X[,i] <-  mvrnorm(1, delta1, S1)
      }
    }
    
    
    
    ### Step 3 ###
    q1 <- 3 # Number of latent factors for Baseline Scan (s=1)
    q2 <- 3 # Number of latent factors for Post-furosemide Scan (s=2)
    
    
    # Generate q(s) X n latent factors for both scans s=1 and s=2
    eta1 <- matrix(rep(0,n*q1),nrow=q1)
    for (i in 1:n){
      eta1[,i] <- mvrnorm(1, gamma1%*%X[,i], diag(1,nrow=q1,ncol=q1))
    }
    eta2 <- matrix(rep(0,n*q2),nrow=q2)
    for (i in 1:n){
      eta2[,i] <- mvrnorm(1, gamma2%*%X[,i], diag(1,nrow=q2,ncol=q2))
    }
    
    
    
    
    ### Step 4 ###
    p <- 15  # Number of Basis Functions
    
    
    # Generate n x p matrix of beta coefficients for both scans (s= 1,2)
    beta1 <- matrix(rep(0,n*p),nrow=n)
    beta2 <- matrix(rep(0,n*p),nrow=n)
    for (i in 1:n) {
      if (C[i] > 0){
        beta1[i,] <- mvrnorm(1,lambda11%*%eta1[,i], diag(Sigma11))
        beta2[i,] <- mvrnorm(1,lambda12%*%eta2[,i], diag(Sigma12))
      }
      if (C[i] < 1){
        beta1[i,] <- mvrnorm(1,lambda01%*%eta1[,i], diag(Sigma01))
        beta2[i,] <- mvrnorm(1,lambda02%*%eta2[,i], diag(Sigma02))
      }
    }
    
    
    
    ### Step 5 ###
    m1 <- length(t1) # Number of time points for baseline scan s=1
    m2 <- length(t2) # Number of time points for post-furosemide scan s=2
    nu <- 3 # bandwidth 
    psi <- 1/p
    
    # Generate basis functions
    b1 <- make_basis(p,nu,psi, t1,t2)$b1
    b2 <- make_basis(p,nu,psi, t1,t2)$b2
    
    
    
    # Generate n X m1 observarions for baseline scan s=1 and n X m2 observations for post-furosemide scan
    y1 <- matrix(rep(0,n*m1),nrow=m1)
    y2 <- matrix(rep(0,n*m2),nrow=m2)
    for (i in 1:n) {
      y1[,i] <- mvrnorm(1,b1%*%beta1[i,],diag(sigma,nrow=m1, ncol=m1))
      y2[,i] <- mvrnorm(1,b2%*%beta2[i,],diag(sigma,nrow=m2, ncol=m2))
    }
    
    ### Step 6 ###
    L <- 3 # Number of Experts
    
    
    Sigma_rho0 <- diag(sig2_rho0, r)
    Sigma_rho1 <- diag(sig2_rho1, r)
    
    rho0 <- mvrnorm(L,mu_rho0, Sigma_rho0) #L X r Coefficient for expert mean given c=0
    rho1 <- mvrnorm(L,mu_rho1, Sigma_rho1) #L X r Coefficient for expert mean given c=1
    
    
    # Generate n X L matrix of underlying continuous expert ratings
    Z <- matrix(rep(0,n*L),nrow=L)
    for (i in 1:n) {
      if (C[i] < 1) {
        Z[,i] <- mvrnorm(1,rho0%*%X[,i],diag(1,L))
      }
      if (C[i] > 0) {
        Z[,i] <- mvrnorm(1,rho1%*%X[,i],diag(1,L))
      }
    }
    
    
    
    ### Step 7 ###
    # Generate a n X L matrix of actual ordinal ratings from the experts
    W <- matrix(rep(0,n*L),nrow=L)
    for (i in 1:L) {
      for (j in 1:n) {
        if (Z[i,j] <= 0) {
          W[i,j] = 0
        }
        if (0 < Z[i,j] & Z[i,j] <= xi[i]) {
          W[i,j] = 1
        } 
        if (xi[i] < Z[i,j]) {
          W[i,j] = 2
        } 
      }
    }
    
    
    
    
    ## Dividing into testing and training dataset
    index <- 1:n
    train <- sort(sample(index,n1))
    test <- index[-train]
    
    
    y1_training <- y1[,train]
    y2_training <- y2[,train]
    y1_testing <- y1[,test]
    y2_testing <- y2[,test]
    
    X_training <- X[,train]
    X_testing <- X[,test]
    
    W_training <- W[,train]
    W_testing <- W[,test]
    
    C_training[,c] <- C[train]
    C_testing[,c] <- C[test]
    
    
    
    # Average and agreement for generated expert ratings
    av  <- apply(W_training,2,mean)
    avv <- cbind(C_training[,c],av)
    avv0 <- avv[which(C_training[,c]=="0"),]
    avv1 <- avv[which(C_training[,c]=="1"),]
    er0[c] <- mean(avv0[,2]) # Mean of three expert ratings in the training data given c=0
    er1[c] <- mean(avv1[,2]) # Mean of three expert ratings in the training data given c=1
    kap[c] <- kappam.fleiss(t(W_training))[[5]] # agreement of three expert ratings in the training data
    
    
    
    
    
    
    
    ### Naive Methods
    
    # Based on Majority Voting
    for (i in 1:n1) {
      mv[i,c] <- majority(W_training[,i]) 
    }
    
    
    
    
    
    # Based on K-mean clustering of C(Y,X)
    mdata <- cbind(t(y1_training),t(y2_training),t(X_training))
    
    set.seed(c*10)
    kmean <-   kmeans(mdata, 2, nstart=5)
    cc <- kmean[[1]]
    kcenter <- kmean[[2]]
    
    rtable_1 <- table(C_training[,c],cc)
    
    if (rtable_1[1,1]+rtable_1[2,2] >= rtable_1[1,2]+rtable_1[2,1]) {
      for (i in 1:n1) {
        if (cc[i]==1){
          kc[i,c] = 0
        } else {
          kc[i,c] =1
        }
        
      }
      
    } else if (rtable_1[1,1]+rtable_1[2,2] < rtable_1[1,2]+rtable_1[2,1]) {
      
      for (i in 1:n1) {
        if (cc[i]==1){
          kc[i,c] = 1
        } else {
          kc[i,c]=0
        }
      }
      A <- kcenter[1,] 
      B <- kcenter[2,]
      kcenter <- rbind(B,A)
    }
    
    
    for (i in 1:n2) {
      diff1 <- euc.dist(c(y1_testing[,i], y2_testing[,i], X_testing[,i]),kcenter[1,])
      diff2 <- euc.dist(c(y1_testing[,i], y2_testing[,i], X_testing[,i]),kcenter[2,])
      if (diff1 < diff2) {
        kct[i,c] = 0
      }
      else {
        kct[i,c] = 1
      }
    }
    
    
    
    
    
    
    
    
    
    ## Based on multivariate functional clustering
    basis1 <- create.bspline.basis(c(0,t1[m1]), nbasis=10)
    basiseval1 <- eval.basis(t1,basis1)
    
    basis2 <- create.bspline.basis(c(0,t2[m2]), nbasis=10)
    basiseval2 <- eval.basis(t2,basis2)
    
    
    
    
    #training 
    fy1_training <- smooth.basis(t1,y1_training,fdParobj=basis1)$fd
    fy2_training <- smooth.basis(t2,y2_training,fdParobj=basis2)$fd
    fy_training <- list(fy1_training,fy2_training)
    
    
    set.seed(c*10)
    fcluster <- funHDDC(fy_training,K=2)
    fclass <- fcluster$class-1
    
    fmean01 <- basiseval1%*%fcluster$mu[1,1:10]
    fmean11 <- basiseval1%*%fcluster$mu[2,1:10]
    
    fmean02 <- basiseval2%*%fcluster$mu[1,11:20]
    fmean12 <- basiseval2%*%fcluster$mu[2,11:20]
    
    
    ftable_1 <- table(C_training[,c],fclass)
    
    if (ftable_1[1,1]+ftable_1[2,2] >= ftable_1[1,2]+ftable_1[2,1]) {
      for (i in 1:n1) {
        if (fclass[i]==1){
          fc[i,c] = 1
        } else {
          fc[i,c] =0
        }
        
      }
      
    } else if (ftable_1[1,1]+ftable_1[2,2] < ftable_1[1,2]+ftable_1[2,1]) {
      
      for (i in 1:n1) {
        if (fclass[i]==1){
          fc[i,c] = 0
        } else {
          fc[i,c]=1
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
    coefs_y1_testing <- smooth.basis(t1,y1_testing,fdParobj=basis1)$fd$coefs
    smooth_y1_testing <- basiseval1%*%coefs_y1_testing
    coefs_y2_testing <- smooth.basis(t2,y2_testing,fdParobj=basis2)$fd$coefs
    smooth_y2_testing <- basiseval2%*%coefs_y2_testing
    
    for (i in 1:n2) {
      fdiff0 <- hilbert_dist(smooth_y1_testing[,i],fmean01,t1, smooth_y2_testing[,i],fmean02,t2)
      fdiff1 <- hilbert_dist(smooth_y1_testing[,i],fmean11,t1, smooth_y2_testing[,i],fmean12,t2)
      if (fdiff0 < fdiff1) {
        fct[i,c] = 0
      }
      else {
        fct[i,c] = 1
      }
    }  
    
    
    print(c)
    
  }
  list(kap=kap,er0=er0, er1=er1, C_training=C_training, C_testing=C_testing,
       mv=mv, kc=kc, kct=kct,fc=fc,fct=fct)
}

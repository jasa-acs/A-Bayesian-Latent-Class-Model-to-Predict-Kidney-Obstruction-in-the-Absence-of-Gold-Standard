# Nuclear Medicine Project
# MCMC software
# written by Changgee Chang
# ver 20180506
#
# Function Nuclear
# Arguments:
# y1, y2: m1 (and m2) by n observation matrix of pre-furosemide (and post-furosemide)
# X: r by n covariates matrix
# W: L by n matrix; Expert Ratings; Values 0, 1, or 2
# t1, t2: m1 (and m2) time points at which the data are collected for pre- (and post-) furosemide
# p, nu, psi: spline basis function parameters; p is the number of basis functions, and see the technical documents for nu and psi
# q: c(q1,q2)
# ups : tuning parameter upsilon
# a : tuning parameters a = c(a1,a2). If NULL, adaptive with initial c(2,2)
# tau : expert variability parameter
# theta: sparsity tuning parameter
#

if ( !require(msm) )
{
  install.packages("msm")
  library(msm)
}

if ( !require(tmvtnorm) )
{
  install.packages("tmvtnorm")
  library(tmvtnorm)
}

Nuclear <- function(nsample,y1,y2,X,W,t1,t2,p,q,nu,psi,ups,a,a_sigma,b_sigma,tau,kappa,theta)
{
  n = ncol(X)
  r = nrow(X)
  L = nrow(W)

  m1 = length(t1)
  m2 = length(t2)
  
  A1 = matrix(1,m1,p)
  A2 = matrix(1,m2,p)
  for ( k in 2:p )
  {
    A1[,k] = exp(-nu*(t1-(k-1)*psi)^2)
    A2[,k] = exp(-nu*(t2-(k-1)*psi)^2)
  }

  AA = array(0,c(p,p,2))
  AA[,,1] = t(A1)%*%A1
  AA[,,2] = t(A2)%*%A2
  
  Ay = array(0,c(p,n,2))
  Ay[,,1] = t(A1)%*%y1
  Ay[,,2] = t(A2)%*%y2
  
  # adaptive a
  if ( is.null(a) )
  {
    ada = TRUE
    a = c(2,2)
  }
  else
    ada = FALSE

  # initialize C
  sumW = apply(W,2,sum)
  C = as.integer(sumW>L)
  U = (sumW==0 | sumW==2*L)
  pC = rep(0,n)
  H = matrix(c(1-C,C),n)
  Cwh = H==1
  N = apply(H,2,sum)

  # initialize sigma2
  beta = array(0,c(p,n,2))
  for ( s in 1:2 )
    beta[,,s] = chol2inv(chol(AA[,,s]+diag(0.01,p)))%*%Ay[,,s]
  eps1 = y1 - A1%*%beta[,,1]
  eps2 = y2 - A2%*%beta[,,2]
  sigma2 = (1+sum(eps1^2)+sum(eps2^2))/(1+n*(m1+m2))

  # initialize Lambda and gamma
  lam = array(0,c(p,max(q),2,2))
  gamma = array(0,c(max(q),r,2))
  iXX0 = chol2inv(chol(X%*%(t(X)*H[,1])+diag(0.01,r)))
  iXX1 = chol2inv(chol(X%*%(t(X)*H[,2])+diag(0.01,r)))
  for ( s in 1:2 )
  {
    tmp0 = beta[,,s]%*%(t(X)*H[,1])%*%iXX0
    tmp1 = beta[,,s]%*%(t(X)*H[,2])%*%iXX1
    tmp = svd(rbind(tmp0,tmp1),q[s],q[s])
    lam[,1:q[s],1,s] = tmp$u[1:p,] %*% diag(tmp$d/sqrt(r),q[s])
    lam[,1:q[s],2,s] = tmp$u[1:p+p,] %*% diag(tmp$d/sqrt(r),q[s])
    if ( q[s] <= r )
      gamma[1:q[s],,s] = t(tmp$v)*sqrt(r)
    else
    {
      gamma[1:r,,s] = t(tmp$v)*sqrt(r)
      gamma[(r+1):q[s],,s] = rnorm((q[s]-r)*r)
    }    
  }    
  
  # initialize eta
  eta = array(0,c(max(q),n,2))
  for ( s in 1:2 )
    eta[1:q[s],,s] = gamma[1:q[s],,s]%*%X
  
  # initialize beta
  for ( s in 1:2 )
  {
    V = chol2inv(chol(AA[,,s]/sigma2+diag(a_sigma/b_sigma,p)))
    for ( i in 1:n )
      beta[,i,s] =  V %*% (Ay[,i,s]/sigma2+a_sigma*lam[,1:q[s],C[i]+1,s]%*%eta[1:q[s],i,s]/b_sigma)
  }
  
  # initialize mu  
  mu = array(0,c(max(q),2,2))
  mu[1,,] = a[1]
  pi = array(a[1],c(max(q),2,2))
  for ( s in 1:2 )
    if ( q[s]>1 )
    {
      mu[2:q[s],,s] = a[2]
      for ( h in 2:q[s] )
        pi[h,,s] = pi[h-1,,s]*mu[h,,s]
    }
  
  # initialize phi    
  phi = (ups+1) / (ups+lam^2*rep(pi,each=p))
  
  # initialize omega
  omega = 2/(gamma^2+1)

  # initialize sig2
  sig2 = array(b_sigma/a_sigma,c(p,2,2))
  for ( s in 1:2 )
    for ( c in 0:1 )
    {
      zeta = beta[,,s] - lam[,1:q[s],c+1,s] %*% eta[1:q[s],,s]
      sig2[,c+1,s] = (1+apply(t(zeta)^2*H[,c+1],2,sum))/(1+N[c+1])
    }
  
  # initialize Z
  Z = W-1/2

  # initialize xi
  xi = matrix(0,L,4)
  xi[,1] = -Inf
  xi[,3] = 1
  xi[,4] = Inf
  
  mu_xi = 1
  sig2_xi = 1/tau

  # initialize rho
  rho = array(0,c(L,r,2))
  rho[,,1] = (Z %*% (t(X)*H[,1])) %*% iXX0
  rho[,,2] = (Z %*% (t(X)*H[,2])) %*% iXX1
  rhoX = array(0,c(L,n,2))
  rhoX[,,1] = rho[,,1] %*% X
  rhoX[,,2] = rho[,,2] %*% X
  
  mu_rho = apply(rho,c(2,3),mean)
  sig2_rho = matrix(1/tau,r,2)
  
  # initialize delta
  delta = X%*%H/rep(N,each=r)

  # initialize S
  S = array(0,c(r,r,2))
  iS = S  
  for ( c in 0:1 )
  {
    S[,,c+1] = (diag(kappa,r)+(X-delta[,c+1])%*%(t(X-delta[,c+1])*H[,c+1]))/(kappa+N[c+1])
    iS[,,c+1] = chol2inv(chol(S[,,c+1]))
  }
  
  # repositories
  rep_beta = array(0,c(p,n,2,nsample))
  rep_lam = array(0,c(p,max(q),2,2,nsample))
  rep_eta = array(0,c(max(q),n,2,nsample))
  rep_gamma = array(0,c(max(q),r,2,nsample))
  rep_C = matrix(1,n,nsample)
  rep_pC = matrix(1,n,nsample)
  rep_phi = array(1,c(p,max(q),2,2,nsample))
  rep_pi = array(1,c(max(q),2,2,nsample))
  rep_omega = array(1,c(max(q),r,2,nsample))
  rep_sig2 = array(0,c(p,2,2,nsample))
  rep_sigma2 = rep(0,nsample)
  rep_mu = array(1,c(max(q),2,2,nsample))
  
  rep_Z = array(0,c(L,n,nsample))
  rep_xi = matrix(0,L,nsample)
  rep_mu_xi = rep(0,nsample)
  rep_sig2_xi = rep(0,nsample)
  rep_rho = array(0,c(L,r,2,nsample))
  rep_mu_rho = array(0,c(r,2,nsample))
  rep_sig2_rho = array(0,c(r,2,nsample))
  rep_delta = array(0,c(r,2,nsample))
  rep_S = array(0,c(r,r,2,nsample))
  rep_iS = rep_S
  
  if ( ada )
    rep_a = matrix(0,2,nsample)
  else
    rep_a = a
  

  for ( it in 1:nsample )
  {
    # beta
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        iV = AA[,,s]/sigma2 + diag(1/sig2[,c+1,s],p)
        ciV = chol(iV)
        iVM = Ay[,Cwh[,c+1],s]/sigma2 + (lam[,,c+1,s]%*%eta[,Cwh[,c+1],s])/sig2[,c+1,s]
        M = backsolve(ciV,forwardsolve(t(ciV),iVM))
        beta[,Cwh[,c+1],s] = M + backsolve(ciV,matrix(rnorm(p*N[c+1]),p))
      }

    # eta
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        iV = t(lam[,1:q[s],c+1,s]) %*% (lam[,1:q[s],c+1,s]/sig2[,c+1,s]) + diag(1,q[s])
        ciV = chol(iV)
        iVM = t(lam[,1:q[s],c+1,s]) %*% (beta[,Cwh[,c+1],s]/sig2[,c+1,s]) + gamma[1:q[s],,s]%*%X[,Cwh[,c+1]]
        M = backsolve(ciV,forwardsolve(t(ciV),iVM))
        eta[1:q[s],Cwh[,c+1],s] = M + backsolve(ciV,matrix(rnorm(q[s]*N[c+1]),q[s]))
      }

    # lambda
    for ( s in 1:2 )
      for ( c in 0:1 )
        for ( k in 1:p )
        {
          iV = eta[1:q[s],,s]%*%(t(eta[1:q[s],,s])*H[,c+1])/sig2[k,c+1,s] + diag(phi[k,1:q[s],c+1,s]*pi[1:q[s],c+1,s],q[s])
          iVM = eta[1:q[s],,s]%*%(beta[k,,s]*H[,c+1])/sig2[k,c+1,s]
          ciV = chol(iV)
          M = backsolve(ciV,forwardsolve(t(ciV),iVM))
          lam[k,1:q[s],c+1,s] = M + backsolve(ciV,rnorm(q[s]))
        }

    # gamma
    for ( s in 1:2 )
      for ( h in 1:q[s] )
      {
        iV = X%*%t(X) + diag(omega[h,,s],r)
        iVM = X%*%eta[h,,s]
        ciV = chol(iV)
        M = backsolve(ciV,forwardsolve(t(ciV),iVM))
        gamma[h,,s] = M + backsolve(ciV,rnorm(r))
      }

    # phi
    for ( s in 1:2 )
      phi[,1:q[s],,s] = rgamma(p*q[s]*2,(ups+1)/2,(ups+lam[,1:q[s],,s]^2*rep(pi[1:q[s],,s],each=p))/2)
    
    # mu
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        if ( q[s] > 1 )
          philam2 = apply(phi[,1:q[s],c+1,s]*lam[,1:q[s],c+1,s]^2,2,sum)
        else
          philam2 = sum(phi[,1,c+1,s]*lam[,1,c+1,s]^2)
        
        for ( d in 1:q[s] )
        {
          pibar = pi[d:q[s],c+1,s]/mu[d,c+1,s]
          if ( ada == FALSE )
          {
            if ( d==1 )
              mu[d,c+1,s] = rgamma(1,a[1]+p*q[s]/2,1+sum(pibar*philam2)/2)
            else
              mu[d,c+1,s] = rgamma(1,a[2]+p*(q[s]-d+1)/2,1+sum(pibar*philam2[d:q[s]])/2)
          }
          else
          {
            if ( d==1 )
              mu[d,c+1,s] = rgig(1,4*a[1],4/a[1]+sum(pibar*philam2),p*q[s]/2-1/2)
            else
              mu[d,c+1,s] = rgig(1,4*a[2],4/a[2]+sum(pibar*philam2[d:q[s]]),p*(q[s]-d+1)/2-1/2)
          }
          pi[d:q[s],c+1,s] = pibar*mu[d,c+1,s]
        }
      }

        
    # a
    if ( ada )
    {
      mu1 = mu[1,,]
      a[1] = rgig(1,8+4*sum(mu1),2+4*sum(1/mu1),3/2)
      mu2 = mu[-1,,]
      a[2] = rgig(1,8+4*sum(mu2),2+4*sum(1/mu2[mu2!=0]),sum(q)-5/2)
    }
        
        
    # omega
    for ( s in 1:2 )
      omega[1:q[s],,s] = rexp(q[s]*r,(gamma[1:q[s],,s]^2+1)/2)

    # sig2
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        zeta = beta[,,s] - lam[,1:q[s],c+1,s] %*% eta[1:q[s],,s]
        sig2[,c+1,s] = 1/rgamma(p,(a_sigma+N[c+1])/2,(b_sigma+apply(t(zeta)^2*H[,c+1],2,sum))/2)
      }
    
    # sigma2
    eps1 = y1 - A1%*%beta[,,1]
    eps2 = y2 - A2%*%beta[,,2]
    sigma2 = 1/rgamma(1,(1+n*(m1+m2))/2,(1+sum(eps1^2)+sum(eps2^2))/2)
    
    # Z
    for ( c in 0:1 )
      for ( l in 1:L )
        Z[l,Cwh[,c+1]] = rtnorm(N[c+1],rhoX[l,Cwh[,c+1],c+1],1,xi[l,W[l,Cwh[,c+1]]+1],xi[l,W[l,Cwh[,c+1]]+2])

    # xi
    for ( l in 1:L )
    {
      lb = max(0,Z[l,W[l,]==1])
      ub = min(Inf,Z[l,W[l,]==2])
      xi[l,3] = rtnorm(1,mu_xi,sqrt(sig2_xi),lb,ub)
    }
    
    # mu_xi
    mu_xi = rtnorm(1,sum(xi[,3])/(sig2_xi+L),sqrt(sig2_xi/(sig2_xi+L)),0)
    while ( TRUE )
    {
      if ( runif(1) < (2*pnorm(mu_xi/sqrt(sig2_xi)))^(-L) )
        break
      mu_xi = rtnorm(1,sum(xi[,3])/(sig2_xi+L),sqrt(sig2_xi/(sig2_xi+L)),0)
    }
    
    # sig2_xi
    sig2_xi = 1/rgamma(1,(tau+L)/2,(1+sum((xi[,3]-mu_xi)^2))/2)
    while ( TRUE )
    {
      if ( runif(1) < (2*pnorm(mu_xi/sqrt(sig2_xi)))^(-L) )
        break
      sig2_xi = 1/rgamma(1,(tau+L)/2,(1+sum((xi[,3]-mu_xi)^2))/2)
    }
    
    # rho
    for ( c in 0:1 )
    {
      iV = X%*%(t(X)*H[,c+1]) + diag(1/sig2_rho[,c+1],r)
      V = chol2inv(chol(iV))
      MiV =  t(X%*%(t(Z)*H[,c+1]) + mu_rho[,c+1]/sig2_rho[,c+1])
      M = MiV %*% V
      cV = kronecker(t(chol(V)),diag(1,L))
      rho[,,c+1] = M + as.vector(cV %*% rnorm(L*r))
      rhoX[,,c+1] = rho[,,c+1] %*% X
    }
    
    # mu_rho
    for ( c in 0:1 )
    {
      M = apply(rho[,,c+1],2,sum)/(sig2_rho[,c+1]+L)
      V = sig2_rho[,c+1]/(sig2_rho[,c+1]+L)
      mu_rho[,c+1] = rnorm(r,M,sqrt(V))
    }

    # sig2_rho
    for ( c in 0:1 )
      sig2_rho[,c+1] = 1/rgamma(r,(tau+L)/2,(1+apply((t(rho[,,c+1])-mu_rho[,c+1])^2,1,sum))/2)

    # delta
    for ( c in 0:1 )
    {
      iV = diag(1,r) + N[c+1]*iS[,,c+1]
      iVM = iS[,,c+1] %*% (X %*% H[,c+1])
      ciV = chol(iV)
      cV = backsolve(ciV,diag(1,r))
      M = backsolve(ciV,forwardsolve(t(ciV),iVM))
      delta[,c+1] = M + cV %*% rnorm(r)
    }
    
    # S
    for ( c in 0:1 )
    {
      V = chol2inv(chol(diag(kappa,r) + (X-delta[,c+1]) %*% (t(X-delta[,c+1])*H[,c+1])))
      iS[,,c+1] = rWishart(1,kappa+N[c+1],V)
      S[,,c+1] = chol2inv(chol(iS[,,c+1]))
    }
    
    # C
    V = array(0,c(p,p,2,2))
    for ( c in 0:1 )
      for ( s in 1:2 )
      {
        iV = AA[,,s]/sigma2 + diag(1/sig2[,c+1,s],p)
        V[,,c+1,s] = chol2inv(chol(iV))
      }
    
    for ( i in 1:n )
    {
      if ( it<500 & U[i] )
        next
      llk = rep(0,2)
      for ( c in 0:1 )
      {
        llk[c+1] = -log(det(S[,,c+1]))/2 - sum(log(sig2[,c+1,]))/2 - c*theta
        llk[c+1] = llk[c+1] - sum((Z[,i]-rhoX[,i,c+1])^2)/2
        llk[c+1] = llk[c+1] - sum((X[,i]-delta[,c+1])*(iS[,,c+1]%*%(X[,i]-delta[,c+1])))/2
        for ( s in 1:2 )
        {
          lameta = lam[,,c+1,s]%*%eta[,i,s]
          iVM = Ay[,i,s]/sigma2 + lameta/sig2[,c+1,s]
          M = V[,,c+1,s]%*%iVM
          llk[c+1] = llk[c+1] + log(det(V[,,c+1,s]))/2 + sum(iVM*M)/2 - sum(lameta^2/sig2[,c+1,s])/2
        }
      }
      pC[i] = 1/(1+exp(llk[1]-llk[2]))
      C[i] = rbinom(1,1,pC[i])
    }
    H = matrix(c(1-C,C),n)
    Cwh = H==1
    N = apply(H,2,sum)
    
    # store
    rep_beta[,,,it] = beta
    rep_lam[,,,,it] = lam
    rep_eta[,,,it] = eta
    rep_gamma[,,,it] = gamma
    rep_C[,it] = C
    rep_pC[,it] = pC
    rep_phi[,,,,it] = phi
    rep_pi[,,,it] = pi
    rep_omega[,,,it] = omega
    rep_sig2[,,,it] = sig2
    rep_sigma2[it] = sigma2
    rep_mu[,,,it] = mu
    if ( ada )
      rep_a[,it] = a
    
    rep_Z[,,it] = Z
    rep_xi[,it] = xi[,3]
    rep_mu_xi[it] = mu_xi
    rep_sig2_xi[it] = sig2_xi
    rep_rho[,,,it] = rho
    rep_mu_rho[,,it] = mu_rho
    rep_sig2_rho[,,it] = sig2_rho

    rep_delta[,,it] = delta
    rep_S[,,,it] = S
    rep_iS[,,,it] = iS
    
    print(it)
  }

  ret = list(nsample=nsample,y1=y1,y2=y2,m=c(m1,m2),n=n,r=r,L=L,X=X,W=W,t1=t1,t2=t2,
        p=p,q=q,nu=nu,psi=psi,ups=ups,tau=tau,kappa=kappa,theta=theta,
        A1=A1,A2=A2,AA=AA,
        beta=rep_beta,lam=rep_lam,eta=rep_eta,gamma=rep_gamma,C=rep_C,pC=rep_pC,
        phi=rep_phi,pi=rep_pi,omega=rep_omega,sig2=rep_sig2,sigma2=rep_sigma2,mu=rep_mu,a=rep_a,
        Z=rep_Z,xi=rep_xi,mu_xi=rep_mu_xi,sig2_xi=rep_sig2_xi,rho=rep_rho,mu_rho=rep_mu_rho,sig2_rho=rep_sig2_rho,
        delta=rep_delta,S=rep_S,iS=rep_iS)
}


# Function Nuclear_predict
# Arguments:
# bi: burn-ins; the first bi samples are ignored
# y1, y2: m1 (and m2) by 1 new observation vector of pre-furosemide (and post-furosemide)
# x: r by 1 covariates vector
# W: Expert Rating; Values 0, 1, 2, or NULL(unavailable,default)
#

Nuclear_predict <- function(fit,bi,y1,y2,x,W=NULL)
{
  p = fit$p
  q = fit$q
  N = fit$nsample
  Ay = cbind(t(fit$A1)%*%y1,t(fit$A2)%*%y2)

  llk = -apply(log(fit$sig2[,,,(bi+1):N]),c(4,2),sum)/2
  llk[,2] = llk[,2] - fit$theta
  
  for ( it in (bi+1):N )
  {
    for ( c in 0:1 )
    {
      for ( s in 1:2 )
      {
        V = matrix(0,p+q[s],p+q[s])
        V[1:p,1:p] = fit$AA[,,s]/fit$sigma2[it] + diag(1/fit$sig2[,c+1,s,it],p)
        V[1:p,p+1:q[s]] = -fit$lam[,1:q[s],c+1,s,it] / fit$sig2[,c+1,s,it]
        V[p+1:q[s],1:p] = t(V[1:p,p+1:q[s]])
        V[p+1:q[s],p+1:q[s]] = -t(fit$lam[,1:q[s],c+1,s,it]) %*% V[1:p,p+1:q[s]] + diag(1,q[s])
        cV = chol(V)
        m = c(Ay[,s]/fit$sigma2[it],fit$gamma[1:q[s],,s,it]%*%x)
        llk[it-bi,c+1] = llk[it-bi,c+1] - sum(log(diag(cV))) + sum(forwardsolve(t(cV),m)^2)/2
      }
      
      if ( is.null(W) )
      {
        pW = 1
      }
      else if ( W==0 )
        pW = pnorm(0,sum(fit$mu_rho[,c+1,it]*x),sqrt(1+sum(x*x*fit$sig2_rho[,c+1,it]))) * pnorm(0,fit$mu_xi[it],sqrt(fit$sig2_xi[it]))
      else if ( W==1 )
      {
        V = matrix(c(1+sum(x*x*fit$sig2_rho[,c+1,it])+fit$sig2_xi[it],-fit$sig2_xi[it],-fit$sig2_xi[it],fit$sig2_xi[it]),2,2)
        m = c(sum(fit$mu_rho[,c+1,it]*x)-fit$mu_xi[it],fit$mu_xi[it])
        pW = ptmvnorm(c(-Inf,0),c(0,Inf),m,V)
        pW = pW - pnorm(0,sum(fit$mu_rho[,c+1,it]*x),sqrt(1+sum(x*x*fit$sig2_rho[,c+1,it]))) * pnorm(0,fit$mu_xi[it],sqrt(fit$sig2_xi[it]))
      }
      else
      {
        V = matrix(c(1+sum(x*x*fit$sig2_rho[,c+1,it])+fit$sig2_xi[it],-fit$sig2_xi[it],-fit$sig2_xi[it],fit$sig2_xi[it]),2,2)
        m = c(sum(fit$mu_rho[,c+1,it]*x)-fit$mu_xi[it],fit$mu_xi[it])
        pW = ptmvnorm(c(0,0),c(Inf,Inf),m,V)
      }
      llk[it-bi,c+1] = llk[it-bi,c+1] + log(as.numeric(pW))
      
      llk[it-bi,c+1] = llk[it-bi,c+1] - log(det(fit$S[,,c+1,it]))/2 - sum((x-fit$delta[,c+1,it])*(fit$iS[,,c+1,it]%*%(x-fit$delta[,c+1,it])))/2
    }
  }
  
  1/(1+exp(llk[,1]-llk[,2]))
}




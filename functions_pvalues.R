#|**********************************************************************;
#  FUNCTIONS used for Shrinkage based Gaussian Graphical Models (GGMs) 
#
#  * Project: EXACT HYPOTHESIS TESTING FOR SHRINKAGE BASED GGM (Bernal et al).
#  * Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Marco Grzegorczyk, Peter Horvatovich
#  * Date created      : 2018-11-12

#***************************************************
# (functions) Shrinkage based Gaussian Graphical Models

# this script has 3 main functions for GGM inference (i.e. partial correlations).
# 1. p values from the shrunk density
# 2. p values from Monte Carlo
# 3. p values from standard density

# r = real partial corr coeficients
# p = number of variables (e.g. genes)
# n = sample size
# lambda= shrinkage value 

#************************************************************************
#  * Revision History  : 30 JAN 2018 replace sd by se=sd/sqrt(sample), and (minor edits e.g. font size)
#  **********************************************************************
#  * Details
#  The shrunk probability density is presented here  
#  For empirical p value with Monte Carlo see  [2] 
#  The standard probability density is presented in [3] 
#  The simulation GGMs with the optimal shrinkage (lambda)
#  is performed with GeneNet [Sch?fer,J. and Strimmer,K. (2005a),  GeneNet 1.2.13. CRAN.]
#************************************************************************;
# Refererences
# [1] Sch?fer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
# [2] Martinez, Wendy L., and Angel R. Martinez. Computational statistics handbook with MATLAB. Chapman and Hall/CRC (2007)
# [3] Fisher, Ronald Aylmer. "The distribution of the partial correlation coefficient." Metron 3 (1924): 329-332.
#**************************************************
library(stats4)

# 1. p values from the SHRUNK density
pcor.function = function(d, lambda, corrected = T){
  
  R = cor(d)
  
  if(corrected == T){
    
    R = sm2vec(R)*(1 + (sm2vec(R)^2)/2/(nrow(d)-1))
    R = vec2sm(R)
    diag(R) = 1
    
  }
  
  C = (1-lambda)*(R)+(lambda)*diag(diag(R)) # shrinkage
  D = chol(x = C) # Cholebsky decomposition
  omega.chol = solve(D, diag(ncol(d))) %*% solve(t(D), diag(ncol(d)))# precision matrix
  
  if(lambda!= 0){
    
    chol.PCOR = -cov2cor( omega.chol )#/(1-z)
    
  } else{ 
    
    chol.PCOR =  -cov2cor( omega.chol ) 
  }
  
  
  # if(corrected == T){
  # 
  #   chol.PCOR = sm2vec(chol.PCOR)*(1 + (sm2vec(chol.PCOR)^2)/2/k.shrunk(p=ncol(d) , n=nrow(d), lambda=z))
  #   chol.PCOR = vec2sm(chol.PCOR)
  # 
  # }
  
  
  diag(chol.PCOR) = 1
  
  return(chol.PCOR)
  
}
# Simulate a GGM under the null hypothesis.
# Use this to estimate the degrees of freedom k
# for the shrunk probability density 
# and compute p values
k.shrunk <- function( p , n,lambda){
  ## Start by simulating null hypothetic data
  sim.pcor.S <- ggm.simulate.pcor(p, 0) # simulate identity GGM
  sim.data.S <- ggm.simulate.data( n, sim.pcor.S) # data from identity GGM (sim.pcor.S)
  GGM.S <- pcor.function(d = sim.data.S, lambda = lambda)#pcor.shrink(sim.data.S,lambda , verbose =FALSE)  ## infer GGM structure from sim.data.S
  r.S <- sm2vec(GGM.S) # vectorize the GGM (symmetric matrix)
  
  # compute the negative log likelihood
  nlogL.shrunk <- function(k) {
    
    #density.shrunk <- function(r.S) {(  ((1)^2-(r.S/(1-lambda))^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))*(1-lambda) )}
    log.density.shrunk <- function(x) {  ((k-3)*0.5)*log(1-(x/(1-lambda))^2) - log( beta(0.5, 0.5*(k-1)))-log(1-lambda) }
    log.f<-log.density.shrunk(r.S)
    
    return(-sum(log.f))
  }
  
  # minimize the negative log likelihood to find the degrees of freedom (k.fit.shrunk)
  k.fit.shrunk <-  mle(nlogL.shrunk, start = list(k = 100), method = "L-BFGS-B", lower = c(5),
                      upper = c(Inf))
  #k.fit.shrunk <-  mle(minuslogl = nlogL.shrunk, start = list(k = 100),method = "BFGS")
  #k.fit.shrunk <-  mle(nlogL.shrunk, start = list(k = 100), method = "Nelder-Mead")
  return(k.fit.shrunk@coef[1])}


p.shrunk <- function( r, p, n ,lambda){
        
        
        k<-mean(unlist(lapply(c(1:25), function(w){ k.shrunk( p , n, lambda) })))
        print(k)
      #   ## Start by simulating null hypothetic data
      #   sim.pcor.S <- ggm.simulate.pcor(p, 0) # simulate identity GGM
      #   sim.data.S <- ggm.simulate.data( n, sim.pcor.S) # data from identity GGM (sim.pcor.S)
      #   GGM.S <- pcor.shrink(sim.data.S,lambda , verbose =FALSE)  ## infer GGM structure from sim.data.S
      #   r.S <- sm2vec(GGM.S) # vectorize the GGM (symmetric matrix)
      #  
      #   # compute the negative log likelihood
      #   nlogL.shrunk <- function(k) {
      #         density.shrunk <- function(r.S) {(  ((1)^2-(r.S/(1-lambda))^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))*(1-lambda) )}
      #         f<-density.shrunk(r.S)
      #         return(-sum(log(f)))
      #   }
      #   
      #   # minimize the negative log likelihood to find the degrees of freedom (k.fit.shrunk)
      # k.fit.shrunk <-  mle(nlogL.shrunk, start = list(k = 100), method = "L-BFGS-B", lower = c(20),
      #                      upper = c(Inf))
      
      #summary(k.fit.shrunk)
     
      ## Compute p values with k.fit.shrunk
        
      if(lambda!=0){  
            density.shrunk <- function(r) {(  ((1)^2-(r/(1-lambda))^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))*(1-lambda)) }
      
            pval.shrunk<-matrix(Inf, length(r),1)
            
            for (i in 1:length(r)){
              int<-  integrate(density.shrunk , lower = -(1-lambda), upper = -abs( r[i] )   )
              pval.shrunk[i]<- 2*(int$value) 
            }
      }else {
            density <- function(r) {(  ((1)^2-(r)^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))) }
            
            pval.shrunk<-matrix(Inf, length(r),1)
            
            for (i in 1:length(r)){
              int<-  integrate(density , lower = -(1), upper = -abs( r[i] )   )
              pval.shrunk[i]<- 2*(int$value) 
            }
            
     }
        
    
    return(pval.shrunk)
}


      
################################
## MONTECARLO p-values

      
  p.montecarlo <- function( r, number, p, n ,lambda){

        # Initialize variables
    r.monte<-matrix(NA, length(r),1) 
    p.values<-matrix(NA,length(r),ncol(r.monte))
    cum.pv<-matrix(0,length(r),1)
    
    # Simulate null hypothetic GGM coefficients for "number" times
      for (i in 1:number){
        
        r.data <- ggm.simulate.data(n,diag(p))
        r.monte.GGM <- pcor.function(d = r.data, lambda = lambda) #ggm.estimate.pcor(r.data, lambda=lambda , verbose = FALSE )
        r.monte <- sm2vec(r.monte.GGM)
        
        # compare the real coefficients against r.monte  
        pv<-sapply(r, function(x) sum( abs(r.monte) >=abs(x) )/length(r.monte))  
        cum.pv<-cum.pv+pv
        
      }
  
     # p values  
     p.monte<-cum.pv/number
     return (p.monte)
}

#********************************************************
# STANDARD p-values 
 
   p.standard <- function( r, p,n){

    # Same procedure as p values from the SHRUNK density
    # using the standard density (i.e. lambda =0)
     
          nlogL.std <- function(k) {
            
            density.std <- function(r) {(  (1-r^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1)))}
            
            f<-density.std(r)
            
            return(-sum(log(f)))
          }
    
    # fit neg log Likelihood to the null hypothetic coefficients
    k.fit.std <-  mle(nlogL.std, start = list(k = 100), method = "L-BFGS-B", lower = c(20),
                         upper = c(Inf))
    
    #summary(k.fit.std)
    
    ## the standard density is
    density.std <- function(r) {(  (1-(r)^2) ^((k.fit.std@coef[1]-3)*0.5)  )/( beta(0.5, 0.5*(k.fit.std@coef[1]-1))) }
    
    ### shrunk p values
    pval.std<-matrix(Inf, length(r),1)
    
          for (i in 1:length(r)){
            int<-  integrate(density.std, lower = -(1), upper = -abs( r[i] )   )
            pval.std[i]<- 2*(int$value) 
          }
    
    
    return(pval.std)
  }
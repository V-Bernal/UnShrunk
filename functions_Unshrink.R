#|**********************************************************************
#  FUNCTIONS to "un-shrink" the partial correlation in Gaussian Graphical Models 
#
#  * Project: "un-shrink" the partial correlation in Gaussian Graphical Models (Bernal et al).
#  * Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Peter Horvatovich, Marco Grzegorczyk.
#  * Date created      : 27 FEB 2020

#***************************************************
# Variables:
# d = data matrix with p columns and n rows
# p = number of variables (e.g. genes)
# n = sample size
# lambda = shrinkage value 

# Note: d a is a n x p matrix of data 
# Note: the p columns of d are the random variables 
# Note: the n rows of d are samples

#************************************************************************
#  * Revision History  : 27 FEB 2020
#************************************************************************;
# Refererences
# [1] Schäfer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
#**************************************************



## The tidyverse style guide
# If files should be run in a particular order, prefix them with numbers.
# Variable names should be nouns
# function names should be verbs.
# Variable and function names should use only lowercase letters, numbers, and _. 
# reserve dots exclusively for the S3 object system
# Always put a space after a comma,
# loops have space between parethesis, function not.
# operators (==, +, -, <-, etc.) should always be surrounded by spaces:
# indent the second line to where the definition starts.
#Comments should be in sentence case, and only end with a full stop if they contain at least two sentences

library(stats4)
library(GeneNet)

unshrink_GGM <- function ( d , PLOT = FALSE ) {
  
  # Note: the shrinkage = lambda
  # Note: d a is a n x p matrix of data 
  # Note: the p columns of d are the random variables 
  # Note: the n rows of d are samples
  
  # Error handling
  
  
  
  # 1. Define a range of 5 equally spaced points between the 0.1 and 1 (the maximum shrinkage)
  cat("static", "\n")
  
  opt_lambda= attr(pcor.shrink(d, verbose = FALSE), "lambda")
  
  cat("optimal shrinkage =",opt_lambda,"\n")
  
  lambda = seq(from = 0.1, to = 1, length.out = 5)
  
  # 2. Compute the shrunk partial correlations [1] for each shrinkage from Step 1
  
  GGM.list<-lapply(lambda , function(z) {
    
    partial.corr <- pcor.shrink(d, z, verbose=FALSE )
    
    return(sm2vec(partial.corr))
    
  }
  )
  
  # 3. Fit a Polynomial for the shrunk partial correlations as function of the shrinkage 
  # Each column is one partial correlation (for different lambdas).
  # Note: poly retrieves error for degrees >20.
  # To (over)fit we want degrees as close to length(lambda)
  
  # Reshape
  pcors = matrix(c(unlist(GGM.list)), nrow = 0.5*(ncol(d)*(ncol(d)-1)),ncol = length(lambda) )
  
  # 4. Polynomial extrapolation to zero shrinkage
  dff = data.frame('lambda' = c(lambda),'pcor' = t(pcors))
  
  degrees = length(lambda) - 1      
  
  unshrunk_pcor = apply(dff, 2 , function(y) {
    
    predict( lm(y ~ poly(lambda, degree = degrees)) ,
             newdata =  data.frame(lambda = 0)) 
    
  }
  
  )
  
  # Remove the label "lambda"
  unshrunk_pcor = unshrunk_pcor[-c(1)]
  
  # Optional: Plot the partial correlation (y axis) vs the shrinkage (x axis)
  if( PLOT ){
    
    plot.ts(t(pcors), plot.type = "single",
            type = "b",
            lw=2, lty = 1:2, col=rainbow(0.5*ncol(d)*(ncol(d)-1)),
            ylab="partial correlation",xlab="shrinkage",
            axes=F, cex=1 , ylim=c(-1, 1), xlim=c(0, length(lambda)))
    
    points(rep(0, length(unshrunk_pcor) ), unshrunk_pcor , pch=20, cex=2, 
           col=rainbow(0.5*ncol(d)*(ncol(d)-1)))
    axis(2)
    axis(1, labels=c(0, signif(lambda, 2)), at= seq(from=0, by=1, length.out=c(length(lambda)+1)) 
    )
    
    box()    
    
  }
  
  # un_pcor is the extrapolated value and is an approximation the "un-shrunk" partial correlations.
  
  return(unshrunk_pcor)
  
  
}


#..................................................


k.shrunk <- function(p , n, lambda) {
  ## Start by simulating null hypothetic data
  sim.pcor.S <- ggm.simulate.pcor(p, 0) # simulate identity GGM
  sim.data.S <- ggm.simulate.data( n, sim.pcor.S) # data from identity GGM (sim.pcor.S)
  GGM.S <- pcor.shrink(sim.data.S, lambda , verbose =FALSE)  ## infer GGM structure from sim.data.S
  r.S <- sm2vec(GGM.S) # vectorize the GGM (symmetric matrix)

  # compute the negative log likelihood
  nlogL.shrunk <- function(k) {
    density.shrunk <- function(r.S) {(  ((1)^2-(r.S / (1-lambda)) ^ 2) ^((k-3) * 0.5)  )/( beta(0.5, 0.5 * (k-1)) * (1-lambda) )}
    f<-density.shrunk(r.S)
    return(-sum(log(f)))
  }

  # minimize the negative log likelihood to find the degrees of freedom (k.fit.shrunk)
  k.fit.shrunk <-  mle(nlogL.shrunk, start = list(k = 100), method = "L-BFGS-B", lower = c(5),
                       upper = c(Inf))

  return(k.fit.shrunk@coef)}



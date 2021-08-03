#------------------
# - Title/Project: The "un-shrunk" partial correlation in Gaussian Graphical Models 
# - Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Peter Horvatovich, Marco Grzegorczyk.
# - Date created      : 27 FEB 2020
# - Revision History  : 29 MARCH 2021, 24 MAY 2021
#------------------
# Description:
# This scripts contains the FUNCTIONS developed for the work "The "un-shrunk" partial correlation in Gaussian Graphical Models" 
#------------------
# Variables:
# - data = data matrix with p columns and n rows
# - p = number of variables (e.g. genes)
# - n = sample size
# - lambda = shrinkage value 
#-----------------
# Notes: 
# - d a is a n x p matrix of data 
# - the p columns of d are the random variables 
# - the n rows of d are samples
#-----------------
# References
# [1] Sch?fer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
#------------------

# Libraries
library(stats4)
library(GeneNet)
library(pracma)
#---------------------------------

unshrink_GGM <- function ( data , l_0 = 0.01,  PLOT = FALSE, corrected = FALSE ) {
  
  #------------------
  # Parameters:
  #------------------
  # lambda = the shrinkage vaue
  # data = the n x p matrix of data 
  # p = the number of random variables, equal to the number of columns of data 
  # n = the number of samples, equal to the number of row of data 
  #------------------
  # Output: 
  #------------------
  # unshrunk_pcor = the approximation of the "un-shrunk" partial correlations
  #------------------
  
  #------------------
  # Error handling
  #------------------
  #
  #
  #------------------
  
  #------------------
  # 1. Define a range of 50 equally spaced points between l_0 and 1 (the maximum shrinkage)
  #------------------
  opt_lambda = attr(pcor.shrink(data, verbose = FALSE), "lambda")
  
  cat("optimal shrinkage =",opt_lambda,"\n")
  
  lambda = seq(from = l_0, to = 1, length.out = 50)
  
  #------------------
  # Compute the condition number with varying shrinkage
  kappa.list<-lapply(lambda , function(z) {

    R = cor(data)
    C = (1-z)*(R)+(z)*diag(diag(R))

    return(rcond(C))

    }

    )
  #------------------
  # function to estimate shrunk partial correlations
  pcor.function = function(z){
    
    # Correlation matrix
    R = cor(data)
    
    # if( corrected == T ){
    #   
    #   R = sm2vec(R)*( 1 + ( sm2vec(R)^2 )/2/(nrow(data)-1) )
    #   R = vec2sm(R)
    #   diag(R) = 1
    #   
    # }
    
    # Ledoit-WOlf shrinkage on the correlation matrix
    C = ( 1 - z )*(R) + (z)*diag( diag(R) ) # shrinkage
    D = chol(x = C) # Cholesky decomposition
    
    # Precision matrix
    omega.chol = solve( D, diag(ncol(data)) ) %*% solve( t(D), diag(ncol(data)) )# precision matrix
    
    # Partial correlations
    chol.PCOR = -cov2cor( omega.chol )
    diag( chol.PCOR ) = 1
    
    return( chol.PCOR )
    
  }

  #------------------
  # 2. Compute the shrunk partial correlations in the shrinkage range from Step 1
  #------------------
  GGM.list<-lapply( lambda , function(z) {
    
      if((z == 1)){ 
        
        partial.corr = sm2vec(matrix(0, nrow = ncol(data), ncol = ncol(data)) ) 
      
        } else {
          
        partial.corr = sm2vec(pcor.function(z)) 
          
               }

        return(partial.corr)
      
    }
  )
  
  #------------------
  # 3. Apply Fisher transformation
  #------------------
  pcors = matrix( c(unlist(GGM.list) ), 
                  nrow = 0.5*( ncol(data)*(ncol(data)-1) ), 
                  ncol = length( lambda ) ) #Reshape. Each pcor (row) with varying lambda (cols)
  
  pcors = atanh(pcors) # fisher trasnform 
  
  #------------------------------------
  # 4. Fit a smoothing splines function to the pcors and extrapolate to zero shrinkage
  #------------------------------------
  dff = data.frame('lambda' = c(lambda),'pcor' = t(pcors))  # Note: Each column is one partial correlation (for different lambdas).
  
  wght = c(unlist(kappa.list))^2.5 # weights
  
  unshrunk_pcor = apply(dff, 2 , function(y) {

                         fit = smooth.spline( x = lambda, y = y , 
                                              w = wght, all.knots = T,
                                              df = 2, keep.data = T)
                         
                         predicted = predict(fit, x =  0 , deriv = 0 )$y

                         return( predicted )
                         
                         }
  
  
                        )
  
  
  unshrunk_pcor = unshrunk_pcor[-c(1)] # Remove the label "lambda"
  
  
  #------------------------------------    
  # 5 . Back transform fisher
  #------------------------------------
  unshrunk_pcor = tanh(unshrunk_pcor)
  
  pcors = tanh(pcors) 
  
  #--------------------------    
  # Optional: Plot the partial correlation (y axis) vs the shrinkage (x axis)
  #--------------------------
  if( PLOT ){
    
    cat( 'Plotting is only recommended for small networks, e.g. 10 nodes' )
    
    plot.ts( t(pcors) , plot.type = "single",
            type = "b",
            lw=2, lty = 1:2, col=rainbow(0.5*ncol(data)*(ncol(data)-1)),
            ylab="partial correlation",xlab="shrinkage",
            axes=F, cex=1 , ylim=c(-1, 1), xlim=c(0, length(lambda)))
    
        points(rep(0, length(unshrunk_pcor) ), unshrunk_pcor , pch=20, cex=2, 
               col = rainbow(0.5*ncol(data)*(ncol(data)-1)))
        axis(2)
        axis(1, labels=c(0, signif(lambda, 2)), at= seq(from=0, by=1, length.out=c(length(lambda)+1)) )
    
        box()    
    
  }
  
  #--------------------------
  # Return the final splines approximation to lambda = 0
  return(unshrunk_pcor)
  
  #--------------------------
  
}

#-------------------------
# ROC and PR curves
#-------------------------
# adapted functions from https://blog.revolutionanalytics.com/2016/08/roc-curves-in-two-lines-of-code.html
#-------------------------
simple_roc<-function(labels,scores){
  labels<-labels[order(scores,decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels),FPR=cumsum(!labels)/sum(!labels),labels)
}

simple_PR<-function(labels,scores){
  labels<-labels[order(scores,decreasing=TRUE)]
  TPR=cumsum(labels)/sum(labels)
  FPR=cumsum(!labels)/sum(!labels)
  PPV=cumsum(labels)/sum(labels)/((cumsum(labels)/sum(labels))+(cumsum(!labels)/sum(!labels)))
  data.frame(TPR, PPV, labels)
}

simple_F1<-function(labels,scores){
  labels<-labels[order(scores,decreasing=TRUE)]
  TPR=cumsum(labels)/sum(labels)
  FPR=cumsum(!labels)/sum(!labels)
  PPV=cumsum(labels)/sum(labels)/((cumsum(labels)/sum(labels))+(cumsum(!labels)/sum(!labels)))
  data.frame(2*(TPR*PPV)/(TPR+PPV), labels)
}

simple_auc<-function(TPR,FPR){
  #inputs already sorted, best scores first
  dFPR<-c(diff(FPR),0)
  dTPR<-c(diff(TPR),0)
  sum(TPR*dFPR)+sum(dTPR*dFPR)/2
}

truncated_simple_auc<-function(TPR,FPR){
  #inputs already sorted, best scores first
  idu = unique(FPR)
  FPR = FPR[idu]
  TPR = TPR[idu]
  #dFPR<-c(diff(FPR),0)
  #dTPR<-c(diff(TPR),0)
  id = approx(x = FPR, y = TPR, xout= seq(from = 0,to = 0.5,length.out = 10) )
  return(trapz(id$x, id$y))
  #sum(TPR*dFPR)+sum(dTPR*dFPR)/2
}

#-------------------------------
roc_fun = function(n){
  cat(n)
  num = 25
  
  # Simulate data
  sim.data1<-
    lapply(X = c(1:num) ,function(X){
      ggm.simulate.data(n,TrueNET)
    })
  
  # Shrunk partial correlations
  GGM1.shrunk<-
    lapply(X = sim.data1,function(X){
      sm2vec(pcor.shrink(X,verbose=FALSE))
    })
  
  # optimal shrinkage value
  shrunk.lambdas<-
    lapply(X = sim.data1,function(X){
      return(attr(pcor.shrink(X,verbose=FALSE),"lambda"))
    })
  
  # Un-shrunk partial correlation
  unshrunk1<-
    lapply(X = sim.data1,function(X){
      unshrink_GGM(d = X, PLOT = F , l_0 = 0.01, corrected = F)
    })
  
  # Glasso partial correlations
  glasso<-
    lapply(sim.data1,function(X){
      -sm2vec(cov2cor(
        huge.select(huge(X,method="glasso"), criterion="stars")$opt.icov
      ))
    })
  
  
  #--------------------------
  # Compute ROC
  simp_roc1<-
    lapply(GGM1.shrunk,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  Usimp_roc1<-
    lapply(unshrunk1,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  glasso_roc1<-
    lapply(glasso,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  # -------------
  # ROC AUC
  shrunk_auc1<-
    lapply(simp_roc1,function(x){
      signif(simple_auc(x$TPR,x$FPR), digits = 4)
    })
  
  Unshrunk_auc1<-
    lapply(Usimp_roc1,function(x){
      signif(simple_auc(x$TPR,x$FPR), digits = 4)
    })
  
  glasso_auc1<-
    lapply(glasso_roc1,function(x){
      signif(simple_auc(x$TPR,x$FPR), digits = 4)
    })
  
  # -------------
  # average auc
  return( list('Unshrunk' = mean(unlist(Unshrunk_auc1)),
               'shrunk' = mean(unlist(shrunk_auc1) ),
               'glasso' = mean(unlist(glasso_auc1)) ,
               'U-S' = t.test(x = unlist(Unshrunk_auc1) , y = unlist(shrunk_auc1) )$p.value ,
               'U-g' = t.test(x = unlist(Unshrunk_auc1) , y = unlist(glasso_auc1) )$p.value  
  )
  )
  
}

#-------------------------------
truncated_roc_fun = function(n){
  cat(n)
  num = 5
  
  # Simulate data
  sim.data1<-
    lapply(X = c(1:num) ,function(X){
      ggm.simulate.data(n,TrueNET)
    })
  
  # Shrunk partial correlations
  GGM1.shrunk<-
    lapply(X = sim.data1,function(X){
      sm2vec(pcor.shrink(X,verbose=FALSE))
    })
  
  # optimal shrinkage value
  shrunk.lambdas<-
    lapply(X = sim.data1,function(X){
      return(attr(pcor.shrink(X,verbose=FALSE),"lambda"))
    })
  
  # Un-shrunk partial correlation
  unshrunk1<-
    lapply(X = sim.data1,function(X){
      unshrink_GGM(d = X, PLOT = F , l_0 = 0.01, corrected = F)
    })
  
  # Glasso partial correlations
  glasso<-
    lapply(sim.data1,function(X){
      -sm2vec(cov2cor(
        huge.select(huge(X,method="glasso"), criterion="stars")$opt.icov
      ))
    })
  
  
  #--------------------------
  # Compute ROC
  simp_roc1<-
    lapply(GGM1.shrunk,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  Usimp_roc1<-
    lapply(unshrunk1,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  glasso_roc1<-
    lapply(glasso,function(x){
      simple_roc(c(sm2vec(TrueNET)!=0), abs(x))
    })
  
  # -------------
  # Truncated ROC AUC
  shrunk_auc1<-
    lapply(simp_roc1,function(x){
      signif( truncated_simple_auc(x$TPR,x$FPR), digits = 4)
    })
  
  Unshrunk_auc1<-
    lapply(Usimp_roc1,function(x){
      signif(truncated_simple_auc(x$TPR,x$FPR), digits = 4)
    })
  
  glasso_auc1<-
    lapply(glasso_roc1,function(x){
      signif(truncated_simple_auc(x$TPR,x$FPR), digits = 4)
    })  
  
  # -------------
  # average auc
  return( list('Unshrunk' = mean(unlist(Unshrunk_auc1)),
               'shrunk' = mean(unlist(shrunk_auc1) ),
               'glasso' = mean(unlist(glasso_auc1)) #,
               #'U-S' = t.test(x = unlist(Unshrunk_auc1) , y = unlist(shrunk_auc1) )$p.value ,
               #'U-g' = t.test(x = unlist(Unshrunk_auc1) , y = unlist(glasso_auc1) )$p.value  
  )
  )
  
}
#-------------------------------
f1_fun = function(n){
  cat(n)
  num = 25
  
  # Simulate data
  sim.data1<-
    lapply(X = c(1:num) ,function(X){
      ggm.simulate.data(n,TrueNET)
    })
  
  # Shrunk partial correlations
  GGM1.shrunk<-
    lapply(X = sim.data1,function(X){
      sm2vec(pcor.shrink(X,verbose=FALSE))
    })
  
  # optimal shrinkage value
  shrunk.lambdas<-
    lapply(X = sim.data1,function(X){
      return(attr(pcor.shrink(X,verbose=FALSE),"lambda"))
    })
  
  # Un-shrunk partial correlation
  unshrunk1<-
    lapply(X = sim.data1,function(X){
      unshrink_GGM(d = X, PLOT = F , l_0 = 0.01, corrected = F)
    })
  
  # Glasso partial correlations
  glasso<-
    lapply(sim.data1,function(X){
      -sm2vec(cov2cor(
        huge.select(huge(X,method="glasso"), criterion="stars")$opt.icov
      ))
    })
  
  
  #--------------------------
  # Compute ROC
  simp_f1<-
    lapply(GGM1.shrunk,function(x){
      simple_F1(c(sm2vec(TrueNET)!=0),abs(x))
    })
  
  Usimp_f1<-
    lapply(unshrunk1,function(x){
      simple_F1(c(sm2vec(TrueNET)!=0),abs(x))
    })
  
  glasso_f1<-
    lapply(glasso,function(x){
      simple_F1(c(sm2vec(TrueNET)!=0),abs(x))
    })
  
 # -------------
  # average auc
  return( list('Unshrunk' = mean(unlist(Usimp_f1)),
               'shrunk' = mean(unlist(simp_f1) ),
               'glasso' = mean(unlist(glasso_f1)) ,
               'U-S' = t.test(x = unlist(Usimp_f1) , y = unlist(simp_f1) )$p.value ,
               'U-g' = t.test(x = unlist(Usimp_f1) , y = unlist(glasso_f1) )$p.value  
  )
  )
  
}

#------------------
# - Title/Project: Rebuttal to the reviewers / The "un-shrunk" partial correlation in Gaussian Graphical Models 
# - Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Peter Horvatovich, Marco Grzegorczyk.
# - Date created      : 24 JUN 2019
# - Revision History  : 24 MAY 2021
#------------------
# Description:
# This scripts contains the analysis for the rebuttal note to the reviewers of BMC bioinformatics
#------------------
# Variables:
# - data = data matrix with p columns and n rows
# - p = number of variables (e.g. genes)
# - n = sample size
# - lambda = shrinkage value
# - number= Montecarlo iterations
# - rep = Times to repeat the simulation (for fixed p, n)
#
# - etaA = proportion of TP
# - alpha = significance level
#-----------------
# Notes: 
# - d a is a n x p matrix of data 
# - the p columns of d are the random variables 
# - the n rows of d are samples
#-----------------
# References
# [1] Sch?fer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
#------------------

sessionInfo()

#------------------
# 1) Analytical results:
#   The shrinkage shifts the Eigenvalues, but the Eigenvectors stay the same.
#------------------

rm(list = ls())

library(GeneNet)
library(MASS)

set.seed(123)

p = 15
n = 10
lambda = 0.5

# Simulate a GGM and data
true_pcor <- ggm.simulate.pcor(p)
m_sim <- ggm.simulate.data(n, true_pcor)
shrunk_cor <- cor.shrink(m_sim, lambda = lambda, verbose = TRUE)

# EigenValues are scaled and shifted (see Equation X)
a <- eigen(cor(m_sim), symmetric = TRUE)$values
aL <- eigen(shrunk_cor, symmetric = TRUE)$values
evals <- data.frame(aL, ((1 - lambda) * a + lambda))
colnames(evals) <- c("Numerical", "Analytical.")
evals

cumsum(a) / sum(a)
cumsum(aL) / sum(aL)
plot(a, pch = 20, ylab = "eigenvalues")
points(aL, pch = 20, ylab = "eigenvalues", col = "red")


# EigenVectors stay same (see Equation X)
eigen(cor(m_sim), symmetric = TRUE)$vectors -
                eigen(shrunk_cor, symmetric = TRUE)$vectors
Null(cor(m_sim))
Null(shrunk_cor)

# the last p-n-1 columns are the null space of cor
cor(m_sim) %*% eigen(cor(m_sim), symmetric = TRUE)$vectors[, 10]

# shrinkage makes the nullspace of corr to have eigenvalue lambda (instead of 0)
shrunk_cor %*% eigen(cor(m_sim), symmetric = TRUE)$vectors[, 10]
lambda * eigen(cor(m_sim), symmetric = TRUE)$vectors[, 10]

# The EigenValues of the INVERSE's are the reciprocal 1/aL (see Equation X )
ainv <-
  eigen(solve(cor(m_sim)), symmetric = TRUE)$values # doesn not exist for n<p

ainvL <- eigen(solve(shrunk_cor), symmetric = TRUE)$values
ievals <- data.frame(sort(ainvL), 1 / ((1 - lambda) * a + lambda))
colnames(ievals) <- c("Numerical", "Analytical Eq.")
ievals

cumsum(ainv) / sum(ainv)
cumsum(ainvL) / sum(ainvL)

plot(ainv, pch = 20, ylab = "eigenvalues")
points(ainvL,
       pch = 20,
       ylab = "eigenvalues",
       col = "green")

# The INVERSE's eigenvectors are the same
eigen(cor(m_sim), symmetric = TRUE)$vectors # sort columns
eigen(solve(cor(m_sim)), symmetric = TRUE)$vectors


eigen(shrunk_cor, symmetric = TRUE)$vectors # sort columns
eigen(solve(shrunk_cor), symmetric = TRUE)$vectors

# The shrunk Inverse (Equation X)
SymInv <-
  eigen(cor(m_sim), symmetric = TRUE)$vectors %*% solve(diag(eigen(shrunk_cor, symmetric =
                                                                     TRUE)$values)) %*% t(eigen(cor(m.sim), symmetric = TRUE)$vectors)

SymInv -
  solve(shrunk_cor) # and numerically are the same

sum(c(SymInv - solve(shrunk_cor)))



#------------------
# 2) Figure 1. The order of the partial correlations changes with lambda: Toy example
#------------------

library(GeneNet)
library(stats4)
library("devEMF")

# Correlation matrix
m <- matrix(c(
  1,
  1/2, -1/4, -1/8,
  1/2,
  1, -3/4, -3/4, -1/4, -3/4,
  1,
  3/4, -1/8, -3/4,
  3/4,
  1 ), nrow = 4, ncol = 4
)
# The toy matrix has eigenvalues, 2-norm condition number, and the real partial correlation 
m
solve(m)*97

eigen(m)$values
kappa(m)
pcor0 <- sm2vec(cor2pcor(m))

# Reconstruct the shrunk partial correlations
lambda <- c(1:40) / 40

pcors <- lapply(lambda , function(z) {
  return(sm2vec(cor2pcor((1 - z) * m + z * diag(ncol(
    m
  )))))
  
})

pcors <-
  matrix(unlist(pcors), 0.5 * ncol(m) * (ncol(m) - 1), length(lambda))


emf(file = paste0("Figure1_a.emf"), width = 5, height = 5 )

  plot.ts(
    t(pcors),
    plot.type = "single",
    type = "l",
    lw = 3,
    lty = 1,
    col = c("grey40","grey30"),# rainbow(0.5 * ncol(m) * (ncol(m) - 1)),
    ylab = "partial correlation",
    xlab = "LW - shrinkage",
    ylim = c(-1, 1),
    axes = F,
    cex = 1
  )
  points(
    rep(0.35, length(pcor0)),
    pcor0,
    pch = 20,
    cex = 2,#,
    bg="grey",
    col = c("grey40","grey30")# rainbow(0.5 * ncol(m) * (ncol(m) - 1))
  )
  axis(2)
  axis(1,
       labels = c(lambda),
       at = seq(
         from = 1,
         by = 1,
         length.out = length(lambda)
       ))
  box()
dev.off()
  
# Compare with Glasso
library(huge)

pcors <- huge(m, method = "glasso", lambda =sort(lambda))
pcors <-sapply(pcors$icov, function(x){ - sm2vec(cov2cor(x))})  
pcors <-
  matrix(unlist(pcors), 0.5 * ncol(m) * (ncol(m) - 1), length(lambda))

emf(file = paste0("Figure_S1_a.emf"), width = 5, height = 5 )

plot.ts(
  t(pcors),
  plot.type = "single",
  type = "l",
  lw = 4,
  lty = 1,
  col = c("grey40","grey30"),
  ylab = "partial correlation",
  xlab = "glasso shrinkage",
  ylim = c(-1, 1),
  axes = F,
  cex = 1
)
points(
  rep(0.35, length(pcor0)),
  pcor0,
  pch = 20,
  cex = 2,#,
  bg="grey",
  col = c("grey40","grey30")
)
axis(2)
axis(1,
     labels = c(lambda),
     at = seq(
       from = 1,
       by = 1,
       length.out = length(lambda)
     ))
box()
dev.off()

#------------------
# Figure 2.: shrinking the correlation is equivalent to shrinking the data
#------------------

rm(list=ls())

library("GeneNet")
library("MASS")
library("devEMF")

set.seed(123)

# Check pcor between node 2 and node 6
p = 8
n = 10
etaA = 0.1
true.pcor <- ggm.simulate.pcor(p, etaA = etaA)
data.sim <- ggm.simulate.data(n, true.pcor)

#data.sim <- scale(data.sim)
# shrink the correlation
#a <- eigen(  cor(data.sim) , symmetric=TRUE)$values # t(data.sim)%*%data.sim/n-1
#s <- svd(data.sim)
#anim.shrunk.data<- function(L){
#  s$u %*% diag(  sqrt(  (n-1)*( (1-L)*(s$d^2) + L*1  ))  ) %*% t(s$v)
#}

# shrink the covariance
anim.shrunk.data<- function(L){
  
  if( L == 0 ){
    return( data.sim )
  }
  else {
    c = (t(data.sim)%*%data.sim)/(n-1) #cov(data.sim ,method = "pearson")
    cL = (1-L) * c + L * diag(diag(c))
    aL <- eigen(  cL , symmetric=TRUE)$values 
    s <- svd(data.sim)
    
    #sqrt((n-1)*eigen(  c, symmetric=TRUE)$values)
    #s$d
    return( s$u %*% diag(  sqrt((n-1)*aL )) %*% t(s$v)  )
  }
}

par(mar = c(4, 4, .1, .5))

    # if(min(eigen(cov(data.sim))$values)<0){
    # 
    # min_L = abs(min(eigen(cov(data.sim))$values)/(min(eigen(cov(data.sim))$values)-1))
    # 
    # }else{
min_L = 0
    #   
    # }
# un-comment to print/save the plot

emf(file = paste0("Figure1_b.emf"), width = 5, height = 5 )

  plot(x = anim.shrunk.data(0)[,2], y = anim.shrunk.data(0)[,3], 
        xlim = max(abs(anim.shrunk.data(0)[,2]))*c(-1.2, 1.5), 
        ylim = max(abs(anim.shrunk.data(0)[,3]))*c(-1.5, 0.8),
        xlab="node 1", ylab="node 2",  
        cex= 1.75 , pch= 20)
  
  lambda = seq(from = min_L , to = 1, length.out = 200) 
  
  for ( i in c(1:length(lambda) ) ){
     
    points(x = anim.shrunk.data(lambda[i])[, 2] , y = anim.shrunk.data(lambda[i])[, 3], 
            pch=20, cex= i/200 , col=rgb(0.31, 0.31, 0.31, 0.3)) #cex=i/10
    
            }
  
  # Optimal shrinkage
  opt.lambda<-attr(cov.shrink(data.sim, verbose=TRUE), "lambda")
  points(anim.shrunk.data(opt.lambda)[,2] ,
          anim.shrunk.data(opt.lambda)[,3],
          col = rgb(0, 0, 0, 1), bg = "white",
          pch = 21, cex = 1.35 )
  
  
  #text(x = 0.6, y = -0.75, labels = paste0("lambda =",round(opt.lambda, 2)), cex = 1)

dev.off()

#------------------
# Extra: Toy example. The approximation to the limit works well
#------------------
rm(list=ls())

library("GeneNet")

source("functions_Unshrink.R")

p = 10
n =  10

set.seed(1)

TrueNet = ggm.simulate.pcor(num.nodes = p, etaA = 0.1)
dat = ggm.simulate.data(sample.size = n, pcor = TrueNet)

opt2 = attr(pcor.shrink(dat, verbose = FALSE), "lambda")

unshrink_GGM (dat ,PLOT = TRUE, l_0 = 0.01 ,corrected = F)
points(y=sm2vec(pcor.shrink(dat))/(1-opt2), x=rep(0.5, length(sm2vec(TrueNet))), 
       cex=1.2, pch="*", col= rainbow(n = length(sm2vec(TrueNet)) ))
text(x=0.5, y=1,"scaled")
text(x=0.2, y=-0.8,"ordered")

#------------------    
# Figure 3. Trend plot: Comparison for a range of partial correlations
#------------------

#(i) p=100 with eta 0.003
#(ii) p=70 , n=50 with eta  
#(ii) p=50, n=40  with eta 0.01 (change the seed)

rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "reshape", 
              "stats4", 
              "devEMF")

lapply(Packages, library, character.only = TRUE)

               
source("functions_Unshrink.R")
         
               
set.seed(123)
               
p <- 100
n1 <- 70 
etaA <- 0.004

times = 25

TrueNET <- ggm.simulate.pcor(p, etaA)
positives.idx <- which(sm2vec(TrueNET)!=0)
non.positives.idx <- which(sm2vec(TrueNET)==0)   


shrunk.pos<-  unshrunk.pos<-
              shrunk.neg<-
              unshrunk.neg<-
              shrunk.pos.se<-
              unshrunk.pos.se<-
              shrunk.neg.se<-
              unshrunk.neg.se<-
              pval.shrunk<-
              pval.unshrunk<-
              opt.lambdas.mean<-c()


for (samples in 1:length(n1)) {
  
   cat("p=", p , "n=", n1[samples], "\n")
  
    # simulate data             
   sim.data1 <-
     lapply(c(1:times), function(x) {
           ggm.simulate.data(n1[samples] , TrueNET)
            })
                 # shrunk
                 GGM1.shrunk <-
                   lapply(sim.data1, function(x) {
                     sm2vec(pcor.shrink(x, verbose = FALSE))
                   })
                 
                 # optimal lambdas
                 opt.lambdas <-
                   lapply(sim.data1, function(x) {
                     return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
                   })
                 
                 # if any lambda was 1, substitute the simulation
                 while (sum(unlist(opt.lambdas) == 1)) {
                       sim.data1 <-
                         lapply(c(1:times), function(x) {
                           ggm.simulate.data(n1[samples] , TrueNET)
                         })
                       
                       GGM1.shrunk <-
                         lapply(sim.data1, function(x) {
                           sm2vec(pcor.shrink(x, verbose = FALSE))
                         })
                       opt.lambdas <-
                         lapply(sim.data1, function(x) {
                           return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
                         })
                       
                 }
                 
                 # unshrunk
                 unshrunk1 <-
                   lapply(sim.data1, function(x) {
                     unshrink_GGM (x, PLOT = FALSE, l_0 =0.01 , corrected = F)
                   })
                 max(abs(unlist(unshrunk1)))
                 
                 # Means
                 shrunk.pos <-
                   cbind(shrunk.pos, rowMeans(sapply(GGM1.shrunk, function(x) {
                     x[positives.idx]
                   })))
                 unshrunk.pos <-
                   cbind(unshrunk.pos, rowMeans(sapply(unshrunk1, function(x) {
                     x[positives.idx]
                   })))
                 shrunk.neg <-
                   cbind(shrunk.neg, rowMeans(sapply(GGM1.shrunk, function(x) {
                     x[non.positives.idx]
                   })))
                 unshrunk.neg <-
                   cbind(unshrunk.neg, rowMeans(sapply(unshrunk1, function(x) {
                     x[non.positives.idx]
                   })))
                 
                 opt.lambdas.mean <-
                   cbind(opt.lambdas.mean, mean(unlist(opt.lambdas)))
                 
                 # new se
                 shrunk.pos.se <- cbind(shrunk.pos.se, apply(sapply(GGM1.shrunk, function(x) {
                   x[positives.idx]
                 }) , 1 ,
                 function(x) {
                   sd(x) / sqrt(times)
                 }))
                 unshrunk.pos.se <- cbind(unshrunk.pos.se, apply(sapply(unshrunk1, function(x) {
                   x[positives.idx]
                 }) , 1 ,
                 function(x) {
                   sd(x) / sqrt(times)
                 }))
                 shrunk.neg.se <- cbind(shrunk.neg.se, apply(sapply(GGM1.shrunk, function(x) {
                   x[non.positives.idx]
                 }) , 1 ,
                 function(x) {
                   sd(x) / sqrt(times)
                 }))
                 unshrunk.neg.se <- cbind(unshrunk.neg.se, apply(sapply(GGM1.shrunk, function(x) {
                   x[non.positives.idx]
                 }) , 1 ,
                 function(x) {
                   sd(x) / sqrt(times)
                 }))

}

#dev.off()

sample.idx = 1
n1[sample.idx]

spider.top<- c(rep(c(as.character(signif(sm2vec(TrueNET)[positives.idx],digits = 3) ), "0")
                    ,1))
df = data.frame(Group = spider.top,
                "Shrunk" = c(shrunk.pos[,sample.idx], mean(shrunk.neg[,sample.idx])) ,
                "Un-Shrunk" = c(unshrunk.pos[,sample.idx],mean(unshrunk.neg[,sample.idx])),
                #"Linear Shrunk" =c(shrunk.pos[,sample.idx], mean(shrunk.neg[,sample.idx]))/(1-opt.lambdas[[1]]),
                "Actual" =c(sm2vec(TrueNET)[positives.idx],0) )

df.se = data.frame(Group = spider.top,
                   "Shrunk" = c(shrunk.pos.se[,sample.idx], mean(shrunk.neg.se[,sample.idx])) ,
                   "Un-Shrunk" = c(unshrunk.pos.se[,sample.idx],mean(unshrunk.neg.se[,sample.idx])) ,
                   #"Linear Shrunk" = c(unshrunk.pos.se[,sample.idx],mean(unshrunk.neg.se[,sample.idx]))/(1-opt.lambdas[[1]])^2,
                   "Actual" =rep(0, length(c(sm2vec(TrueNET)[positives.idx],0)) ))

df.m <- melt(df, 
             id.vars= c("Group"), 
             measure.vars= c("Shrunk", "Un.Shrunk","Actual"),
             variable.name= "Color",
             value.name=    "val"
)
df.mse <- melt(df.se, 
               id.vars= c("Group"), 
               measure.vars= c("Shrunk", "Un.Shrunk","Actual"),
               variable.name= "Color",
               value.name=    "val"
)


df.m<-df.m[order(as.numeric(as.character(df.m$Group))), ]
df.m$Group<-as.factor(df.m$Group)
df.mse<-df.mse[order(as.numeric(as.character(df.mse$Group))), ]
df.mse$Group<-as.factor(df.mse$Group)

df_fin<-cbind(df.m, df.mse$val)
colnames(df_fin)<-c("Label",  "Method","ave","SE")
df_fin$ave<- df_fin$ave
df_fin$SE<-2*df_fin$SE
df_fin$Label <- factor(df_fin$Label, levels = unique(df_fin$Label))

fig<-ggplot(data = df_fin,  aes(x = Label, y = ave, group = Method, colour = Method, fill=NA)) + 
  geom_line(size= 1, aes(linetype=Method)) +
  geom_point(size= 2) +
  geom_errorbar(aes(ymin= ave - SE, ymax= ave + SE , color=Method), width=0.1, size=1.5)+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.line = element_line(size = 1, colour = "grey80"),
        #axis.title =     element_blank(),
        axis.title.x = element_text(color = "grey20", size = 14),
        axis.title.y = element_text(color = "grey20", size = 14),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "black"))  +
        scale_color_manual(values= c( rgb(0.5, 0.5, 0.5, 0.75), rgb(0.25, 0.25, 0.25, 0.75) , rgb(0, 0, 0, 1) ))+
        scale_fill_manual(values= c( rgb(0.5, 0.5, 0.5, 0.75), rgb(0.25, 0.25, 0.25, 0.1) , rgb(0, 0, 0, 0.1) )) +
        scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+
        scale_linetype_manual(values=c("dashed", "dashed", "solid"))


fig + annotate(geom="text", x= 2, y= 0.85, label = paste0("lambda = ", signif(opt.lambdas.mean[sample.idx],2) ),
                 color="black", size=5 ) +
      annotate(geom="text", x= 2, y= 0.7, label = paste0("n = ",n1[sample.idx]),
                 color="black", size=5)+ 
      annotate(geom="text", x= 2, y= 0.55, label = paste0("p = ", p ),
                                                  color="black", size=5 ) +
      geom_hline(yintercept = 0, linetype="solid", color = "grey")+
      labs(x ="Actual pcor", y = "Inferred pcor")


#emf(file = paste0("Figure3_p",p,"n",n1,".emf"), width = 7, height = 5 )
  fig+ annotate(geom="text", x= 4, y= 0.85, label = paste0("lambda = ", signif(opt.lambdas.mean[sample.idx],2) ),
                color="black", size=5 ) +
    annotate(geom="text", x= 2, y= 0.7, label = paste0("n = ",n1[sample.idx]),
             color="black", size=5)+ 
    annotate(geom="text", x= 2, y= 0.55, label = paste0("p = ", p ),
             color="black", size=5 ) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black")+
    labs(x ="Actual pcor", y = "Inferred pcor")
dev.off()


#------------------   
# Figure 4: Distance to the actual value
#------------------

rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "huge", 
              "stats4", 
              "devEMF",
              "reshape",
              "Hmisc",
              "reshape2")

lapply(Packages, library, character.only = TRUE)

source("functions_Unshrink.R")

#-----------------------------------------
# This is a modified version of function ggm.simulate.pcor (GeneNet)
# It simulates a network with fixed number (instead of of percentage) of edges
ggm.simulate.pcor.modified = function (num.nodes, num.edges) {
    eps = 1e-04
    element.idx = data.frame('row'= seq(1,  num.edges, 2), 
                             'col'= 1 + seq(1, num.edges, 2))
    precision = matrix(0, nrow = num.nodes, ncol = num.nodes)
    precision[ cbind(element.idx$row , element.idx$col ) ] = c(-1+2*rbinom(n = 0.5*num.edges, size = 1,  prob = 0.5)) * seq(from = 0.3, to = 5, length.out = 0.5*num.edges ) 
    precision[lower.tri(precision)] = 0
    diag(precision)=1
    for (i in 1:(num.nodes - 1)) {
      for (j in (i + 1):num.nodes) {
        precision[j, i] = precision[i, j]
      }
    }
    for (i in 1:num.nodes) {
      diag(precision)[i] = sum(abs(precision[, i])) + eps
    }
    pcor = cov2cor(precision)
    return(pcor)
  }
#-----------------------------------------

set.seed(12)

p<-c( 5 )*10   
n<-c( 1:9)*10   
num.edges <- 10

all.lambda<-matrix(NA, length(n),length(p)) 
true.pcor<- ggm.simulate.pcor.modified(p, num.edges)

positives<-which(sm2vec(true.pcor)!= 0)
negatives<-which(sm2vec(true.pcor) == 0)
temp = data.frame(matrix(Inf, nrow = length(positives),
                         ncol= length(n)))



temp= cbind(temp, sm2vec(true.pcor)[positives])
colnames(temp)<-c(n,"Actual")
diff.shrunk.pos<-diff.shrunk.SE<-temp
diff.un.shrunk.pos<-diff.un.shrunk.SE<-temp

for (s.size in 1:length(n)){
  cat('p=', p , '  n=', n[s.size] , '/n')
    # Average the partial corr over several simulations
    repeated = 10
    
    for (i in c(1:repeated)){
      
      null.data<-ggm.simulate.data(n[s.size],true.pcor)
      
          while(attr(pcor.shrink(  null.data, verbose = FALSE  ), 'lambda') > 0.95){
            null.data<-ggm.simulate.data(n[s.size],true.pcor)
          }

      assign( paste0('GGM',i)  , pcor.shrink(  null.data, verbose=FALSE  ))
      assign( paste0('r',i)  , sm2vec(pcor.shrink(  null.data, verbose=FALSE  )))
      
      assign( paste0('UNGGM',i)  , unshrink_GGM(  null.data, PLOT=FALSE, l_0 = 0.01, corrected = F  ))
      
    }
    
    all.lambda[s.size]<- mean(attr(GGM1, 'lambda'),
                                    attr(GGM2, 'lambda'),
                                    attr(GGM3, 'lambda'),
                                    attr(GGM4, 'lambda'),
                                    attr(GGM5, 'lambda'),
                                    attr(GGM6, 'lambda'),
                                    attr(GGM7, 'lambda'),
                                    attr(GGM8, 'lambda'),
                                    attr(GGM9, 'lambda'),
                                    attr(GGM10, 'lambda'))
    

    # L1 distance  Unshrunk  diff.un.shrunk.pos[ s.size, times]  
    diff.un.shrunk.pos[ ,s.size] <- rowMeans( cbind(  c(UNGGM1[positives]  ),
                                               c(UNGGM2[positives] ) ,
                                               c(UNGGM3[positives] ),
                                               c(UNGGM4[positives]  ),
                                               c(UNGGM5[positives] ),
                                               c(UNGGM6[positives] ),
                                               c(UNGGM7[positives] ),
                                               c(UNGGM8[positives] ),
                                               c(UNGGM9[positives] ),
                                               c(UNGGM10[positives])), na.rm = TRUE)  
    
    diff.un.shrunk.SE[ ,s.size]<-(1/sqrt(repeated))* sd(  c(  c(UNGGM1[positives]  ),
                                                                   c(UNGGM2[positives] ) ,
                                                                   c(UNGGM3[positives] ),
                                                                   c(UNGGM4[positives]  ),
                                                                   c(UNGGM5[positives] ),
                                                                   c(UNGGM6[positives] ),
                                                                   c(UNGGM7[positives] ),
                                                                   c(UNGGM8[positives] ),
                                                                   c(UNGGM9[positives] ),
                                                                   c(UNGGM10[positives])))

    # L1 distance  shrunk   
    diff.shrunk.pos[ ,s.size]<-rowMeans( cbind(  c(r1[positives]  ),
                                              c(r2[positives] ) ,
                                              c(r3[positives] ),
                                              c(r4[positives]  ),
                                              c(r5[positives] ),
                                              c(r6[positives] ),
                                              c(r7[positives] ),
                                              c(r8[positives] ),
                                              c(r9[positives] ),
                                              c(r10[positives])), na.rm = TRUE)   
    
  diff.shrunk.SE[ ,s.size]<-(1/sqrt(repeated))* sd(  c(  c(r1[positives]  ),
                                                              c(r2[positives] ) ,
                                                              c(r3[positives] ),
                                                              c(r4[positives]  ),
                                                              c(r5[positives] ),
                                                              c(r6[positives] ),
                                                              c(r7[positives] ),
                                                              c(r8[positives] ),
                                                              c(r9[positives] ),
                                                              c(r10[positives])))
  

}


# Shrunk
df <- melt(diff.shrunk.pos)
df$rowid <- c(1:length(positives))
colnames(df) = c('samples', 'partial_correlation', 'rowid')
df$samples = factor( df$samples, levels = unique(df$samples))

diff.shrunk.SE = melt(diff.shrunk.SE)
diff.shrunk.SE$rowid <- c(1:length(positives))
colnames(diff.shrunk.SE) = c('samples', 'se', 'rowid')
diff.shrunk.SE$samples = factor( diff.shrunk.SE$samples, levels = unique(diff.shrunk.SE$samples))

df = merge(df ,diff.shrunk.SE)


emf(file = paste0("Figure3_pvsn_S",p,".emf"), width = 7, height = 5 )

fig = ggplot(df, aes(samples, partial_correlation, group = factor(rowid))) +
  
  geom_point(data = df[!c(df$samples=='Actual'),],
             aes(size = 2*se , shape = factor(rowid)) , alpha = 0.75)+
  
  geom_line(data = df[!c(df$samples=='Actual'),],
            size = 1, alpha = 0.75 )+
  
  geom_point(data = df[df$samples=='Actual',],
             aes(samples, partial_correlation, shape = factor(rowid)
             ), size= 4 )+
  
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.line = element_line(size = 1, colour = "grey80"),
        axis.title.x = element_text(color = "grey20", size = 14),
        axis.title.y = element_text(color = "grey20", size = 14),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        legend.position = "none"
  )+
  scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+ 
  annotate(geom="text", x= 1, y= 0.55, label = paste0("p = ", p ),
           color="black", size=5 )

fig 
dev.off()

# Un-shrunk 
df <- melt(diff.un.shrunk.pos)
df$rowid <- 1:length(positives)
colnames(df) = c('samples', 'partial_correlation', 'rowid')
df$samples=as.factor( df$samples )

diff.un.shrunk.SE = melt(diff.un.shrunk.SE)
diff.un.shrunk.SE$rowid <- c(1:length(positives))
colnames(diff.un.shrunk.SE) = c('samples', 'se', 'rowid')
diff.un.shrunk.SE$samples = factor( diff.un.shrunk.SE$samples, levels = unique(diff.un.shrunk.SE$samples))

df = merge(df ,diff.un.shrunk.SE)


emf(file = paste0("Figure3_pvsn_US",p,".emf"), width = 7, height = 5 )

fig = ggplot(df, aes(samples, partial_correlation, group = factor(rowid))) +
  
  geom_point(data = df[!c(df$samples=='Actual'),],
             aes(size = 2*se , shape = factor(rowid)) , alpha = 0.75)+
  
  geom_line(data = df[!c(df$samples=='Actual'),],
            size = 1, alpha = 0.75 )+
  
  geom_point(data = df[df$samples=='Actual',],
             aes(samples, partial_correlation, shape = factor(rowid)
             ), size= 4 )+
  
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.line = element_line(size = 1, colour = "grey80"),
        #axis.title =     element_blank(),
        axis.title.x = element_text(color = "grey20", size = 14),
        axis.title.y = element_text(color = "grey20", size = 14),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        legend.position = "none") +
  scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+ 
  annotate(geom="text", x= 1, y= 0.55, label = paste0("p = ", p ),
           color="black", size=5 )

fig

dev.off()

#------------------
# 5) Real data: E. coli example
#------------------
               
rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "igraph", 
              "stats4", 
              "devEMF",
              "reshape",
              "STRINGdb",
              "limma")
               
source("functions_Unshrink.R")
source("functions_pvalues.R")


#Load data
data(ecoli)
               
p <- ncol(ecoli)
n <- nrow(ecoli)
names <- colnames(ecoli)
ecoli <- scale(ecoli)


# Estimate GGM (i.e. the partial correlation coefficients)
unshrunk.pcor <- unshrink_GGM(d = ecoli,PLOT = F, l_0 =0.01, corrected = F )
opt.pcor <- sm2vec(pcor.shrink(ecoli))
opt.lambda <- attr(pcor.shrink(ecoli), "lambda")

# Compare the partial correlation coefficients. 
p.UNSHRUNK <- abs(unshrunk.pcor)
p.opt <- abs(opt.pcor) 

id <- sm.index(matrix(0, p , p), diag = FALSE)
df <- data.frame(id[, 1], id[, 2])

edge <-
  data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))

scat.pval <- 
  data.frame(edge, p.opt, p.UNSHRUNK)
colnames(scat.pval) <- c("Edges", "shrunk", "unshrunk")

scat.pval$Edges <- as.character(scat.pval$Edges)
scat.pval$Edges[which(scat.pval$shrunk < 1 &
                        scat.pval$unshrunk < 1)] <- NA

# plot the partial correlation coefficients. 
fig<-ggplot(scat.pval, aes(x = shrunk, y = unshrunk)) + 
  geom_point(aes(x = scat.pval$shrunk, y = scat.pval$unshrunk), size = 1) +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0.1, linetype="dashed", color = "black")+
  geom_vline(xintercept= 0.1, linetype="dashed", color = "black")+
  geom_text(aes(label = scat.pval$Edges), 
            size = 5, nudge_x = 0, nudge_y = 0)+
  labs(x = "|shrunk pcor|", y = "|un-shrunk pcor|") +
  geom_abline(
    intercept = 0,
    slope = 1,
    size = 1 ,
    linetype = "dashed",
    color = "grey") 


# Uncomment to save as .emf
#emf(file = paste0("Ecoli_scatter.emf"), width = 4, height = 4 )
fig
#dev.off()
fig+scale_x_log10()+ scale_y_log10()

# Compare the p values [2]
rm(p.UNSHRUNK ,p.opt)
p.UNSHRUNK <- p.standard(unshrunk.pcor, ncol(ecoli), nrow(ecoli))
p.opt <- p.shrunk(opt.pcor, ncol(ecoli), nrow(ecoli), opt.lambda)
            
               id <- sm.index(matrix(0, p , p), diag = FALSE)
               df <- data.frame(id[, 1], id[, 2])
               edge <-
                 data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))
               scat.pval <- 
                 data.frame(edge,-log10(p.opt), -log10(p.UNSHRUNK))
               colnames(scat.pval) <- c("Edges", "shrunk", "unshrunk")
               
               scat.pval$Edges <- as.character(scat.pval$Edges)
               scat.pval$Edges[which(scat.pval$shrunk < 100 &
                                       scat.pval$unshrunk < 100)] <- NA
               
               # plot the p values.
               fig<-ggplot(scat.pval, aes(x = shrunk, y = unshrunk) ) + 
                 geom_point(aes(x = scat.pval$shrunk, y = scat.pval$unshrunk), size = 1) +
                 theme_classic(base_size = 14, base_family = "") +
                 geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "black")+
                 geom_vline(xintercept= -log10(0.05), linetype="dashed", color = "black")+
                 geom_text(aes(label = scat.pval$Edges), 
                   size = 5, nudge_x = 0, nudge_y = 0)+
                 labs(x = "-log10 pvals shrunk", y = "-log 10 pvals un-shrunk") +
                  geom_abline(
                    intercept = 0,
                    slope = 1,
                    size = 1 ,
                    linetype = "dashed",
                    color = "grey")
 
               # uncomment to save as .emf
               #emf(file = paste0("Ecoli_pval.emf"), width = 4, height = 4 )
               fig
               #dev.off()
               
               fig+scale_x_log10()+ scale_y_log10()


               
# Plot the network with Benjamini Hochberg adjusted p values
sum(p.opt < 0.05)
sum(p.UNSHRUNK < 0.05)

hist(p.opt )
hist(p.UNSHRUNK, add=T, col=rgb(1, 0, 0, 0.5))

hist(p.adjust(p = p.opt, method = 'BH' ) )
hist(p.adjust(p = p.UNSHRUNK, method = 'BH' ) , add=T, col=rgb(1, 0, 0, 0.5))

sum(p.adjust(p = p.opt, method = 'BH' ) < 0.05)
sum(p.adjust(p = p.UNSHRUNK, method = 'BH' ) < 0.05)

res = data.frame('magnitude' = unshrunk.pcor,
                 'pvalue'= p.standard(r = unshrunk.pcor,
                                      p = p, n = n)
)
res2 = data.frame('magnitude' = opt.pcor,
                  'pvalue'=  p.shrunk(opt.pcor, ncol(ecoli), nrow(ecoli), opt.lambda))

with(res, plot(magnitude,-log10(pvalue), pch=20, main="Volcano plot",
               xlim=c(-1,1), ylim=c(0,8),  col=rgb(0,0,0,0.25) ))
with(res2, points(magnitude,-log10(pvalue), pch=20, main="Volcanoplot",
                  xlim=c(-1,1), col=rgb(1,0,0,0.15)))

abline(h = -log10(.05), lty=2)
abline(v = c(-0.1,0.1), lty=2)

GGM <- vec2sm( p.adjust(p = p.opt, method = 'BH' ) < 0.05) +
   2 * vec2sm( p.adjust(p = p.UNSHRUNK, method = 'BH' ) < 0.05)
 diag(GGM) <- 0
 
 
 new_names = names
 unconnected = which(rowMeans(GGM)==0)
 
 if(length(unconnected)>0){
   GGM = GGM[-unconnected, -unconnected]
   rowMeans(GGM)
   new_names = names[-unconnected]
 }
 
 g1 <-
   graph_from_adjacency_matrix(
     abs(GGM),
     mode = c("undirected"),
     diag = FALSE,
     weighted = TRUE
   )
 V(g1)$label <- as.matrix(names[-unconnected])
 colEdge <- c("blue", "red", "magenta")
 line.thickness <- c(1, 1, 3)
 
plot(
   g1,
  layout = layout_nicely(g1),#layout.sphere(g1),
   edge.width = line.thickness[E(g1)$weight], #0.5 * E(g1)$weight,
   edge.color = colEdge[E(g1)$weight],
   vertex.label.color = "black",
   vertex.label.dist = - 0.55,
   vertex.label.cex = 1 ,
   vertex.color = c("black"),
   vertex.size = 2 ,
   edge.label.family = "Times"
 )
               # legend(
               #   x = -1,
               #   y = -1,
               #   c("Shrunk", "Unshrunk", "Overlap"),
               #   pch = 21,
               #   col = "#777777",
               #   pt.bg = colEdge,
               #   pt.cex = 0.8,
               #   cex = 1,
               #   bty = "n",
               #   ncol = 1
               # )  
 

 #emf(file = paste0("FigureS5_ecoli_net.emf"), width = 6, height = 6 )
 plot(
   g1,
   layout = layout_nicely(g1),#layout.sphere(g1),
   edge.width = line.thickness[E(g1)$weight], #0.5 * E(g1)$weight,
   edge.color = colEdge[E(g1)$weight],
   vertex.label.color = "black",
   vertex.label.dist = - 0.55,
   vertex.label.cex = 1 ,
   vertex.color = c("black"),
   vertex.size = 2 ,
   edge.label.family = "Times"
 )
 dev.off()
               
       # STRINGdb
       dat<-data.frame(colnames(ecoli))
       # Number of edges
       # sum(p.adjust(p = p.opt, method = 'BH', n = length(p.opt) ) < 0.05)
       # sum(p.adjust(p =p.UNSHRUNK, method = 'BH', n = length(p.UNSHRUNK)  )  < 0.05)
       
       # continue with un adjusted 0.05
       
       # Number of nodes
       sum(p.opt<=0.05)
       sum(p.UNSHRUNK<=0.05)
       sig.shrunk = names[ unique( sm.index(diag(p))[ p.opt<=0.05   ]  ) ]
       sig.unshrunk =  names[ unique( sm.index(diag(p))[ p.UNSHRUNK<=0.05 ]  ) ]
       length(sig.shrunk)
       length(sig.unshrunk)
       # present in set 1 but not in set 2
       setdiff(sig.unshrunk, sig.shrunk)
       setdiff(sig.shrunk, sig.unshrunk)
       
       
       dat.unshrunk<-data.frame(sig.unshrunk)
       dat.shrunk<-data.frame(sig.shrunk)
       
       colnames(dat)<-"gene"
       colnames(dat.unshrunk)<-"gene"
       colnames(dat.shrunk)<-"gene"
       head(dat.shrunk)
       
       string_db <- STRINGdb$new( version="10", species=362663, 
                                  score_threshold=0, input_directory="" )
       
       full_ecoli <- string_db$map(dat , "gene", removeUnmappedRows = TRUE )
       backgroundV <- full_ecoli$STRING_id 
       string_db$set_background(backgroundV)
       string_db <- STRINGdb$new( version="10",input_directory="", score_threshold=0, species=362663, backgroundV = backgroundV )
       
       ecoli_mapped.unshrunk <- string_db$map(dat.unshrunk , "gene", removeUnmappedRows = TRUE )
       ecoli_mapped.shrunk <- string_db$map(dat.shrunk , "gene", removeUnmappedRows = TRUE )
       
       hits.unshrunk <- ecoli_mapped.unshrunk$STRING_id
       hits.shrunk <- ecoli_mapped.shrunk$STRING_id
       
       getOption("SweaveHooks")[["fig"]]()
       string_db$plot_network( full_ecoli$STRING_id )
       string_db$plot_network( hits.unshrunk )
       string_db$plot_network( hits.shrunk )
       
       string_db$plot_ppi_enrichment( hits.unshrunk, quiet=TRUE )
       string_db$plot_ppi_enrichment( hits.shrunk, quiet=TRUE )
       
       eh <- string_db$enrichment_heatmap( list( hits.unshrunk, hits.shrunk),
                                            list("list1","list2"), title="My Lists" )
       
       enrichmentGO.unshrunk <- string_db$get_enrichment( hits.unshrunk, category = "Process", methodMT = "fdr", iea = TRUE )
       enrichmentGO.shrunk <- string_db$get_enrichment( hits.shrunk, category = "Process", methodMT = "fdr", iea = TRUE )
       
       enrichmentKEGG.unshrunk <- string_db$get_enrichment( hits.unshrunk, category = "KEGG", methodMT = "fdr", iea = TRUE )
       enrichmentKEGG.shrunk <- string_db$get_enrichment( hits.shrunk, category = "KEGG", methodMT = "fdr", iea = TRUE )
       
       # present in set 1 but not in set 2
       setdiff(enrichmentGO.unshrunk$term_description[enrichmentGO.unshrunk$pvalue_fdr<0.05],
               enrichmentGO.shrunk$term_description[enrichmentGO.shrunk$pvalue_fdr<0.05])
       setdiff(enrichmentGO.shrunk$term_description[enrichmentGO.shrunk$pvalue_fdr<0.05],
               enrichmentGO.unshrunk$term_description[enrichmentGO.unshrunk$pvalue_fdr<0.05])
       
       setdiff(enrichmentKEGG.unshrunk$term_description[enrichmentKEGG.unshrunk$pvalue_fdr<0.05], 
               enrichmentKEGG.shrunk$term_description[enrichmentKEGG.shrunk$pvalue_fdr<0.05])
       setdiff(enrichmentKEGG.shrunk$term_description[enrichmentKEGG.shrunk$pvalue_fdr<0.05], 
               enrichmentKEGG.unshrunk$term_description[enrichmentKEGG.unshrunk$pvalue_fdr<0.05])
       
       
#       write.table(x = enrichmentGO.unshrunk , file = 'GOecoli_UNSHRUNK.txt', quote = F, sep = '\t', col.names = T, row.names = F)
#       write.table(x = enrichmentGO.shrunk, file = 'GOecoli_SHRUNK.txt', quote = F, sep = '\t', col.names = T, row.names = F)
#       write.table(x = enrichmentKEGG.unshrunk , file = 'KEGGecoli_UNSHRUNK.txt', quote = F, sep = '\t', col.names = T, row.names = F)
#       write.table(x = enrichmentKEGG.shrunk , file = 'KEGGecoli_SHRUNK.txt', quote = F, sep = '\t', col.names = T, row.names = F)
       

#------------------
# Real data 2: Botomly
#------------------

rm(list=ls())
#BiocManager::install('biomaRt' )
#BiocManager::install('Biobase' )
#BiocManager::install('limma' )

Packages <- c("GeneNet", 
              "ggplot2",
              "Biobase",
              "igraph", 
              "stats4", 
              "devEMF",
              "reshape",
              "STRINGdb",
              "limma")
       


source("functions_Unshrink.R")
source("functions_pvalues.R")
cohen_criteria = 0.1              

pval.cutoff = 0.05
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)

bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)

# filter low mean, low median per strain
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)

# first 10 vs last 11
pdata_bot$strain

## Differential expression at p- value threshold 
#fit limma model
# perform a differential expression analysis using limma with 
# only the strain variable as an outcome. 
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)

#adjust p value
limma_output = topTable(ebayes_limma, number = dim(edata)[1], adjust.method="BH", sort="none")
names(limma_output)
limma_pvals_adj = limma_output$adj.P.Val
limma_pvals_adj[1:10]

hist(limma_pvals_adj, col = 2)
hist(p.adjust(ebayes_limma$p.value[, 2], 'BH'), add=T, col='blue')
plot(limma_pvals_adj, p.adjust(ebayes_limma$p.value[, 2], 'BH'))

sum(limma_pvals_adj < pval.cutoff)

genes = as.integer(limma_pvals_adj < pval.cutoff)
names(genes) = rownames(edata)
sum(is.na(genes))
not_na = !is.na(genes)
genes = genes[not_na]
head(genes)
sum(genes)


p = sum(genes)
n = ncol(edata)
names = names(genes)[which(limma_pvals_adj < pval.cutoff)]

##  Upper Quantile normalization
# Normalisation is dividing the values by the column-based upper quartile statistic. 
# The matrix needs to be transposed first (with columns as genes) to get the arithmetic

# Un-comment for quantile normalization
#norm_edata = normalize.quantiles(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ]))
#norm_edata = t(normalize.quantiles(t(as.matrix(norm_edata))))
boxplot(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ]))
data.quantileExpressed <- apply(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ])
                    , 2, function(x){quantile(x, 0.75)});


data.norm <- t(t(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ])) / data.quantileExpressed);
boxplot(as.matrix(data.norm))


## Strain effects
plot(data.norm[1, ],col=as.numeric(pdata_bot$strain)) # level shift?
svd1 = svd(data.norm - rowMeans(data.norm))
plot(svd1$v[, 1],svd1$v[, 2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata_bot$strain))

# Correct Strain effects
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(data.norm, mod)
fit_limma$coefficients

plot(data.norm[1, ],col=as.numeric(pdata_bot$strain), xlim=c(0, 22))
abline(c(fit_limma$coefficients[1, 1],0))

yy = data.norm[,as.numeric(pdata_bot$strain) == 2] - fit_limma$coefficients[, 2]
data.norm[,as.numeric(pdata_bot$strain) == 2] = yy

plot(data.norm[2, ],col=as.numeric(pdata_bot$strain),  xlim=c(0, 22))# level shift corrected

svd1 = svd(data.norm - rowMeans(data.norm))
plot(svd1$v[, 1],svd1$v[, 2],xlab = "PC1",ylab = "PC2",
     col=as.numeric(pdata_bot$strain))

## Network analysis
data.norm = t(scale(t(data.norm)))
opt.pcor <-sm2vec(pcor.shrink( t(data.norm ), verbose = FALSE ))
opt.lambda <-attr(pcor.shrink( t(data.norm ), verbose = FALSE ),"lambda")
unshrunk.pcor <- unshrink_GGM( d = t(data.norm), PLOT = FALSE, l_0 = 0.01,corrected = F )

cor(opt.pcor, unshrunk.pcor)
sum(opt.pcor - unshrunk.pcor)^2

idx.shrunk = order( abs(opt.pcor) , decreasing = TRUE)
idx.unshrunk = order( abs(unshrunk.pcor), decreasing = TRUE)

cor(idx.shrunk , idx.unshrunk, method = "kendall")
head(data.frame('shrunk'= idx.shrunk , 'unshrunk'=idx.unshrunk) , n = 10) 

# Compare the coefficients
p.UNSHRUNK <- abs(unshrunk.pcor)
p.opt <- abs(opt.pcor) 

          id <- sm.index(matrix(0, p , p), diag = FALSE)
          df <- data.frame(id[, 1], id[, 2])
          
          short_names = substr(names, start=12, stop=length(names))#names
          
          edge <-
           data.frame(paste (short_names[df[, 1]], short_names[df[, 2]], sep = "-", collapse = NULL))
          scat.pval <- 
           data.frame(edge, p.opt, p.UNSHRUNK)
          colnames(scat.pval) <- c("Edges", "shrunk", "unshrunk")
          
          scat.pval$Edges <- as.character(scat.pval$Edges)
          scat.pval$Edges[which(scat.pval$shrunk < 100 &
                                 scat.pval$unshrunk < 100)] <- NA
          
          fig<-ggplot(scat.pval, aes(x = shrunk, y = unshrunk)) + 
           geom_point(aes(x = scat.pval$shrunk, y = scat.pval$unshrunk), size = 1) +
           geom_hline(yintercept = 0.1 , linetype="dashed", color = "black")+
           geom_vline(xintercept =  0.1 , linetype="dashed", color = "black")+
           theme_classic(base_size = 14, base_family = "") +
           geom_text(aes(label = scat.pval$Edges), 
                     size = 5, nudge_x = 0, nudge_y = 0.1)+
           labs(x = "|shrunk pcor|", y = "|un-shrunk pcor|") +
           geom_abline(
             intercept = 0,
             slope = 1,
             size = 1 ,
             linetype = "dashed",
             color = "grey") #+
            #scale_y_continuous(limits =c(0, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))
          
          
          #emf(file = paste0("MM_scatter.emf"), width = 4, height = 4 )
          fig
          dev.off()
          # and
          fig+scale_x_log10()+ scale_y_log10()

          
# Compare the p values

rm(p.UNSHRUNK ,p.opt)
p.UNSHRUNK <- p.standard(unshrunk.pcor, ncol(t(data.norm)), nrow(t(data.norm)))
p.opt <- p.shrunk(opt.pcor, ncol(t(data.norm)), nrow(t(data.norm)), opt.lambda)
          
          id <- sm.index(matrix(0, p , p), diag = FALSE)
          df <- data.frame(id[, 1], id[, 2])
          edge <-
            data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))
          scat.pval <- 
            data.frame(edge,-log10(p.opt), -log10(p.UNSHRUNK))
          colnames(scat.pval) <- c("Edges", "shrunk", "unshrunk")
          
          scat.pval$Edges <- as.character(scat.pval$Edges)
          scat.pval$Edges[which(scat.pval$shrunk < 100 &
                                  scat.pval$unshrunk < 100)] <- NA
          
          fig<-ggplot(scat.pval, aes(x = shrunk, y = unshrunk),
                      ylim = c(0, max(scat.pval$shrunk , scat.pval$unshrunk)),
                      xlim = c(0, max(scat.pval$shrunk ,  scat.pval$unshrunk))) + 
            geom_point(aes(x = scat.pval$shrunk, y = scat.pval$unshrunk), size = 1) +
            theme_classic(base_size = 14, base_family = "") +
            geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "black")+
            geom_vline(xintercept= -log10(0.05), linetype="dashed", color = "black")+
            geom_text(aes(label = scat.pval$Edges), 
                      size = 5, nudge_x = 0, nudge_y = 0)+
            labs(x = "-log10 pvals shrunk", y = "-log10 pvals un-shrunk") +
            geom_abline(
              intercept = 0,
              slope = 1,
              size = 1 ,
              linetype = "dashed",
              color = "grey")   
         
          
          
          #emf(file = paste0("MM_pval.emf"), width = 4, height = 4 )
          fig
          dev.off()
          # and
          fig+scale_x_log10()+ scale_y_log10()      
                    
##  Plot the network 
sum(p.opt<=0.001)
sum(p.UNSHRUNK<=0.001)
          
hist(p.adjust(p = p.opt, method = 'BH', n = length(p.opt) ) )
hist(p.adjust(p = p.UNSHRUNK, method = 'BH', n = length(p.UNSHRUNK) ) , add=T, col=rgb(1, 0, 0, 0.5))

hist(p.opt, 20 )
hist(p.UNSHRUNK , add=T, col=rgb(1, 0, 0, 0.5))          

res = data.frame('magnitude' = unshrunk.pcor,
                 'pvalue'= p.standard(r = unshrunk.pcor,
                                      p = p, n = n)
)
res2 = data.frame('magnitude' = opt.pcor,
                  'pvalue'=  p.shrunk(r = opt.pcor, n = n, p = p, lambda = opt.lambda))

with(res, plot(magnitude,-log10(pvalue), pch=20, main="Volcano plot",
               xlim=c(-1,1), ylim=c(0,8),  col=rgb(0,0,0,0.25) ))
with(res2, points(magnitude,-log10(pvalue), pch=20, main="Volcanoplot",
                  xlim=c(-1,1), col=rgb(1,0,0,0.15)))
abline(h = -log10(0.001), lty=2)
abline(v = c(-0.1,0.1), lty=2)

GGM <- vec2sm( p.opt < 0.001) +
 2 * vec2sm( p.UNSHRUNK  < 0.001)
diag(GGM) <- 0

new_names = names
unconnected = which(rowMeans(GGM)==0)
if(length(unconnected)>0){
 GGM = GGM[-unconnected, -unconnected]
 rowMeans(GGM)
 new_names = names[-unconnected]
}

# Ensembl are mapped to external gene names
library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
mouse_gene_ids = new_names
all_new_gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
                  filters = 'ensembl_gene_id', values = mouse_gene_ids, mart = ensembl)


g1 <-
 graph_from_adjacency_matrix(
   abs(GGM),
   mode = c("undirected"),
   diag = FALSE,
   weighted = TRUE
 )
V(g1)$label <- as.matrix(substr(new_names, start=12, stop=length(new_names)))#all_new_gene$external_gene_name
colEdge <- c("blue", "red", "magenta")
line.thickness <- c(1, 1, 3)


#emf(file = paste0("FigureS6_MM_net.emf"), width = 6, height = 6 )

plot(
 g1,
 layout = layout_nicely(g1),#layout_as_tree(g1)
 edge.width = line.thickness[E(g1)$weight],#0.5 * E(g1)$weight,
 edge.color = colEdge[E(g1)$weight],
 vertex.label.color = "black",
 vertex.label.dist = - 0.55,
 vertex.label.cex = 0.75 ,
 vertex.color = c("black"),
 vertex.size = 1 ,
 edge.label.family = "Times"
)
#dev.off()

# legend(
#  x = -1,
#  y = -1,
#  c("Shrunk", "Unshrunk", "Overlap"),
#  pch = 21,
#  col = "#777777",
#  pt.bg = colEdge,
#  pt.cex = 0.8,
#  cex = 1,
#  bty = "n",
#  ncol = 1
# )                

#...................................................
# Number of edges. 
#sum(p.adjust(p = p.opt, method = 'BH' ) < 0.05)
#sum(p.adjust(p =p.UNSHRUNK, method = 'BH' )  < 0.05)

sum(p.opt < 0.001)
sum(p.UNSHRUNK  < 0.001)

# Number of nodes
sig.shrunk = names[ unique( sm.index(diag(p))[p.opt < 0.001] ) ]
#all_new_gene$external_gene_name[ unique( sm.index(diag(p))[p.opt < 0.001]  ) ]
length(sig.shrunk)
sig.unshrunk =  names[ unique( sm.index(diag(p))[ p.UNSHRUNK  < 0.001 ]  ) ]
#all_new_gene$external_gene_name[ unique( sm.index(diag(p))[ p.UNSHRUNK  < 0.001 ]  ) ]
length(sig.unshrunk)
# present in set 1 but not in set 2
setdiff(sig.unshrunk, sig.shrunk)
setdiff(sig.shrunk, sig.unshrunk)


# STRINGdb
u =get_STRING_species(version="10", species_name=NULL)
names(u)
which(get_STRING_species(version="10", species_name=NULL)$compact_name=='Mus musculus')
u[116,]
10090

dat<-data.frame(names)
dat.unshrunk<-data.frame(sig.unshrunk)
dat.shrunk<-data.frame(sig.shrunk)

dat.shrunk = as.data.frame(dat.shrunk[!is.na(dat.shrunk)])
dat.unshrunk= as.data.frame(dat.unshrunk[!is.na(dat.unshrunk)])

colnames(dat)<-"gene"
colnames(dat.unshrunk)<-"gene"
colnames(dat.shrunk)<-"gene"
head(dat.shrunk)

string_db <- STRINGdb$new( version="10", species=10090, 
                           score_threshold=0, input_directory="" )

full_mmus <- string_db$map(dat , "gene", removeUnmappedRows = TRUE )
backgroundV <- full_mmus$STRING_id 
string_db$set_background(backgroundV)
string_db <- STRINGdb$new( version="10",input_directory="", score_threshold=0, species=10090, backgroundV = backgroundV )

mmus_mapped.unshrunk <- string_db$map(dat.unshrunk , "gene", removeUnmappedRows = TRUE )
mmus_mapped.shrunk <- string_db$map(dat.shrunk , "gene", removeUnmappedRows = TRUE )

hits.unshrunk <- mmus_mapped.unshrunk$STRING_id
hits.shrunk <- mmus_mapped.shrunk$STRING_id

 getOption("SweaveHooks")[["fig"]]()
 string_db$plot_network( full_mmus$STRING_id )
 string_db$plot_network( hits.unshrunk )
 string_db$plot_network( hits.shrunk )
# 
# string_db$plot_ppi_enrichment( hits.unshrunk, quiet=TRUE )
# string_db$plot_ppi_enrichment( hits.shrunk, quiet=TRUE )
# 
# eh <- string_db$enrichment_heatmap( list( hits.unshrunk, hits.shrunk),
#                                     list("list1","list2"), title="My Lists" )

enrichmentGO.unshrunk <- string_db$get_enrichment( hits.unshrunk, category = "Process", methodMT = "fdr", iea = TRUE )
enrichmentGO.shrunk <- string_db$get_enrichment( hits.shrunk, category = "Process", methodMT = "fdr", iea = TRUE )

enrichmentKEGG.unshrunk <- string_db$get_enrichment( hits.unshrunk, category = "KEGG", methodMT = "fdr", iea = TRUE )
enrichmentKEGG.shrunk <- string_db$get_enrichment( hits.shrunk, category = "KEGG", methodMT = "fdr", iea = TRUE )

# present in set 1 but not in set 2
setdiff(enrichmentGO.unshrunk$term_description[enrichmentGO.unshrunk$pvalue_fdr<0.05],
        enrichmentGO.shrunk$term_description[enrichmentGO.shrunk$pvalue_fdr<0.05])
setdiff(enrichmentGO.shrunk$term_description[enrichmentGO.shrunk$pvalue_fdr<0.05],
        enrichmentGO.unshrunk$term_description[enrichmentGO.unshrunk$pvalue_fdr<0.05])

setdiff(enrichmentKEGG.unshrunk$term_description[enrichmentKEGG.unshrunk$pvalue_fdr<0.05], 
        enrichmentKEGG.shrunk$term_description[enrichmentKEGG.shrunk$pvalue_fdr<0.05])
setdiff(enrichmentKEGG.shrunk$term_description[enrichmentKEGG.shrunk$pvalue_fdr<0.05], 
        enrichmentKEGG.unshrunk$term_description[enrichmentKEGG.unshrunk$pvalue_fdr<0.05])

# unshrunk has 906 gos

#......................................................
# ..............Supplementary material..................
#........................................................


#......................................
# Figure S1: a linear correction is not enough
# ......................................
library(GeneNet)
library(stats4)
library(devEMF)

# Correlation matrix
m <- matrix(c(
  1,
  1 / 2,-1 / 4,-1 / 8,
  1 / 2,
  1,-3 / 4,-3 / 4,-1 / 4,-3 / 4,
  1,
  3 / 4,-1 / 8,-3 / 4,
  3 / 4,
  1
),
4,
4)
# The toy matrix has eigenvalues, 2-norm condition number, and the real partial correlation 
m
solve(m)*97

eigen(m)$values
kappa(m)
pcor0 <- sm2vec(cor2pcor(m))

# Reconstruct the shrunk partial correlations
lambda <- c(1:20) / 20
lambda =  lambda[-length(lambda)]
pcors <- lapply(lambda , function(z) {
  return(sm2vec(cor2pcor((1 - z) * m + z * diag(ncol(
    m
  )))))
  
})

pcors <-
  matrix(unlist(pcors), 0.5 * ncol(m) * (ncol(m) - 1), length(lambda))

pcors <-
  matrix(unlist(pcors), 0.5 * ncol(m) * (ncol(m) - 1), length(lambda))
pcors = pcors%*%diag(1/(1-lambda) )
pcor0 <-sm2vec(m)

emf(file = paste0("FigureS1.emf"), width = 5, height = 5 )

plot.ts(
  t(pcors),
  plot.type = "single",
  type = "l",
  lw = 4,
  lty = c(1:2),
  col = "grey40",# rainbow(0.5 * ncol(m) * (ncol(m) - 1)),
  ylab = "partial correlation",
  xlab = "LW - shrinkage",
  ylim = c(-1, 1),
  axes = F,
  cex = 1
)
points(
  rep(0.35, length(pcor0)),
  pcor0,
  pch = 20,
  cex = 2,#,
  bg="grey",
  col = "grey40"# rainbow(0.5 * ncol(m) * (ncol(m) - 1))
)
axis(2)
axis(1,
     labels = c(lambda),
     at = seq(
       from = 1,
       by = 1,
       length.out = length(lambda)
     ))
box()
dev.off()

#........................................
# Figure S2) and Supplementary Table S1
# Compare the order of the partial cor for different lambdas against the actual GGM
# let's check the top 10 most strong ones
#........................................

rm(list=ls())


library(igraph)
library(GeneNet)
library(ggplot2)     
library(stats4) 
source("functions_pvalues.R")
source("functions_Unshrink.R")


set.seed(123)

p<- 10
n<- 1000

etaA<-0.7
#top= floor(0.25*0.5*p*(p-1))

lambda<-   c(1:40)/40 

TrueNET <- ggm.simulate.pcor(p, etaA)
true.idx<- sort( abs(sm2vec(TrueNET)),decreasing = TRUE, index.return=TRUE)$ix

sim.data <- ggm.simulate.data( n , TrueNET)

opt.lambda<- attr(pcor.shrink(sim.data),"lambda")

V<-eigen(cor(sim.data), symmetric=TRUE)$vectors
D<-eigen(cor(sim.data), symmetric=TRUE)$values

which(D==0)
diag(1/((1-0)*D+0))

# with shrinkage
GGM <-lapply( lambda , function(L) { sm2vec(-cov2cor(V %*% diag(1/((1-L)*D+L)) %*%t(V)) )})
idx<- lapply(GGM, function(x) { sort( abs( x),decreasing = TRUE, index.return=TRUE)$ix})

# un shrinking
GGM_U <- unshrink_GGM(sim.data, PLOT = F)  
idx_U<-  sort( abs( GGM_U),decreasing = TRUE, index.return=TRUE)$ix

# table S1
indices<-data.frame(true.idx, idx_U, idx )
colnames(indices)<-c("true" ,"Un_shrunk", lambda)

# order
indices

#while the optimal lambda is small 
opt.lambda

# the order
top = length(GGM_U)
SPEARMAN <-lapply( lambda , function(L) { cor(true.idx, idx[[ signif(10*L, 1) ]], 
                                              method="spearman" )  })
SPEARMAN2 <-lapply( c(1:top), function(x) { cor(true.idx[1:x], idx[[ signif(10*opt.lambda, 1) ]][1:x], 
                                                method="spearman" )  })
SPEARMAN2_U <-lapply( c(1:top), function(x) { cor(true.idx[1:x], idx_U[1:x], 
                                                  method="spearman" )  })

plot(y = SPEARMAN , x = lambda, ylim=c(-1, 1) )
lines(y = rep( cor(true.idx, idx_U, method="spearman" ), length(lambda)), x= lambda)
abline(v= opt.lambda)

plot(y = SPEARMAN2 , x = c(1:top), ylim=c(-1, 1) )
lines(y = c(unlist(SPEARMAN2_U)) , x= c(1:top))
lines(y = abs(sm2vec(TrueNET)[true.idx]) , x= c(1:top), col='blue')
text(x=10, y=-0.5, "actual pcor", col='blue')

# 
# cor(true.idx, idx_U, method="spearman" )
# cor(true.idx, idx[[1]], method="spearman" )
# 
# cor(true.idx, idx_U, method="kendall" )
# cor(true.idx, idx[[1]], method="kendall" )

write.table(x= indices, file = "ordering_p10_n1000_eta_50%_optL_02.txt", sep= '\t', row.names = F, col.names = T)

#............. True GGM
gtrue<-graph_from_adjacency_matrix(TrueNET!=0, mode = c("undirected"), diag=FALSE )
plot(gtrue, layout = layout_in_circle(gtrue),edge.width= 20*abs(TrueNET[sm2vec(TrueNET)!=0]),
     vertex.label.color="black", vertex.label.dist=3, 
     vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	)
text(-1.2, -1.2,  expression(True_GGM),
     cex=1.5, pos=3, col="black") 

#..........top 25% edges for lambda GGM
top10_1 = lapply(idx, function(x) { x[1:top]})
nontop = lapply(idx, function(x) { x[-c(1:top)]})
topGGM<- lapply(c(1:length(lambda)), function(x) { 
  t<-GGM[[x]] 
  t[ c(nontop[[x]]) ]=0
  t<- vec2sm(t)
  return (t)}
)

sapply(c(1:length(lambda)), function(x) {
  print(x)
  g1<-graph_from_adjacency_matrix(topGGM[[x]]!=0, mode = c("undirected"), diag=FALSE )
  
  #win.metafile(paste0(x,"supp1.wmf"))
  
  plot(g1, layout = layout_in_circle(g1),edge.width= 20*abs(topGGM[[x]][topGGM[[x]]!=0]),
       vertex.label.color="black", vertex.label.dist=1.5, 
       vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	)
  text(-1.2, -1.2,  paste("lambda =", eval(lambda[x])),
       cex=1.5, pos=3, col="black") 
  
  #dev.off()
})               





#..................................................................     
# Figure S3 Heatmap: 
# Comparison of L1 distances to the actual partial correlation value
##................................................................
rm(list = ls())

library(GeneNet)
library(reshape)
library(ggplot2)
library(stats4)
library(Hmisc)


setwd("/home/victor/Desktop/")
source("functions_Unshrink.R")

set.seed(123)

# Initialize parameters
p<-c( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) # num of genes   
n<-c( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) # num of samples  

etaA <- 0  # proportion of  TP
#lambda<-0.3

all.lambda<-matrix(NA, length(n),length(p)) 
diff.shrunk.pos<-diff.shrunk.neg<-matrix(NA, length(n),length(p))
diff.un.shrunk.pos<-diff.un.shrunk.neg<-matrix(NA, length(n),length(p))

for (times in 1:length(p)){
  
  for (s.size in 1:length(n)){
    
    cat('p=', p[times] , '  n=', n[s.size] , '/n')
    
    true.pcor<-ggm.simulate.pcor(p[times], etaA)
    
    cat('create a pcor of', 0.5,"\n") 
    true.pcor[1, 2] = true.pcor[2, 1] = 0.5
    
    positives<-which(sm2vec(true.pcor)!=0)
    negatives<-which(sm2vec(true.pcor)==0)
    
    # Average the partial corr over several simulations
    
    for (i in 1:10){
      
      null.data<-ggm.simulate.data(n[s.size],true.pcor)
      
      while(attr(pcor.shrink(  null.data, verbose=FALSE  ), 'lambda') > 0.95){
        null.data<-ggm.simulate.data(n[s.size],true.pcor)
      }
      
      assign( paste0('GGM',i)  , pcor.shrink(  null.data, verbose=FALSE  ))
      assign( paste0('r',i)  , sm2vec(pcor.shrink(  null.data, verbose=FALSE  )))
      assign( paste0('UNGGM',i)  , unshrink_GGM (  null.data, PLOT=FALSE  ))
    }
    

    all.lambda[s.size, times]<- mean(attr(GGM1, 'lambda'),
                                    attr(GGM2, 'lambda'),
                                    attr(GGM3, 'lambda'),
                                    attr(GGM4, 'lambda'),
                                    attr(GGM5, 'lambda'),
                                    attr(GGM6, 'lambda'),
                                    attr(GGM7, 'lambda'),
                                    attr(GGM8, 'lambda'),
                                    attr(GGM9, 'lambda'),
                                    attr(GGM10, 'lambda'))
    

    # L1 distance  Unshrunk   
    diff.un.shrunk.pos[ s.size, times]<-(1/length(positives))* mean( sum( abs(c(UNGGM1 - sm2vec(true.pcor) )[positives] ) ),
                                                                    sum(abs( c(UNGGM2 - sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM3 - sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM4  - sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM5 - sm2vec(true.pcor))[positives])) ,
                                                                    sum(abs( c(UNGGM6 -sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM7 -sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM8  -sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM9 -sm2vec(true.pcor))[positives]) ),
                                                                    sum(abs( c(UNGGM10 - sm2vec(true.pcor))[positives])), na.rm = TRUE)  
    
    # diff.un.shrunk.neg[ s.size, times]<-(1/length(negatives))* mean(   sum(abs( c(UNGGM1 - sm2vec(true.pcor))[negatives] )) ,
    #                                                                 sum(abs(  c(UNGGM2 - sm2vec(true.pcor))[negatives]) ),
    #                                                                 sum(abs(  c(UNGGM3 - sm2vec(true.pcor))[negatives]) ),
    #                                                                 sum(abs(  c(UNGGM4  - sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM5 - sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM6 -sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM7 -sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM8  -sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM9 -sm2vec(true.pcor))[negatives])),
    #                                                                 sum(abs(  c(UNGGM10 - sm2vec(true.pcor))[negatives])), na.rm = TRUE)        
    # L1 distance  shrunk   
    diff.shrunk.pos[ s.size, times]<-(1/length(positives))* mean( sum(abs( c(r1 - sm2vec(true.pcor))[positives] )) ,
                                                                 sum(abs(  c(r2 - sm2vec(true.pcor))[positives]) ),
                                                                 sum(abs(  c(r3 - sm2vec(true.pcor))[positives]) ),
                                                                 sum(abs(  c(r4  - sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r5 - sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r6 -sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r7 -sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r8  -sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r9 -sm2vec(true.pcor))[positives])),
                                                                 sum(abs(  c(r10 - sm2vec(true.pcor))[positives])), na.rm = TRUE)  
    
    # diff.shrunk.neg[ s.size, times]<-(1/length(negatives))* mean(sum( c(r1 - sm2vec(true.pcor))[negatives] )) ,
    #                                                             sum(abs(  c(r2 - sm2vec(true.pcor))[negatives]) ),
    #                                                             sum(abs(  c(r3 - sm2vec(true.pcor))[negatives]) ),
    #                                                             sum(abs(  c(r4  - sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r5 - sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r6 -sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r7 -sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r8  -sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r9 -sm2vec(true.pcor))[negatives])),
    #                                                             sum(abs(  c(r10 - sm2vec(true.pcor))[negatives])), na.rm = TRUE)        
    
    
  }
}

save.image(file = "heatmap.R")
## Positives ...............

min.colour<- 0
max.colour<- 1
num.col<-5


# cols= p (number of genes) , rows are samples size
num_genes <- paste(rep("p", length(p)), p, sep="=") # rows
sample_size <- paste(rep("n", length(n)), n, sep="=") # columns
mc.ENF <- data.frame(genes = num_genes, 
                     t( diff.un.shrunk.pos), nrow = length(p), ncol =length(n) )

#mc.ENF$genes <- factor(mc.ENF$genes, levels = sort(unique(mc.ENF$genes)))

names(mc.ENF )[2:(length(n)+1)] <- sample_size 
mc.ENF<-mc.ENF[,-c(ncol(mc.ENF)-1, ncol(mc.ENF))]
df_heatmap.gnet <- melt(mc.ENF , id.vars = "genes")
names(df_heatmap.gnet)[2:3] <- c("sample", "P")
head(df_heatmap.gnet)
df_heatmap.gnet$genes <- factor(df_heatmap.gnet$genes, levels = unique(mc.ENF$genes))


fig = ggplot(df_heatmap.gnet, aes(x=sample, y=genes  )) +
  geom_tile(aes(fill = P), color = "white") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red" ,  midpoint = 0, limits=c(min.colour, max.colour))+
  geom_text(aes(label = round(df_heatmap.gnet$P, 3)),size=4) +
  ylab("") +
  xlab("") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size=14),
        axis.text.y = element_text(size=14) )+
  labs(fill = "|Un-Shrunk - 0.5|")  

library("devEMF")
emf(file = paste0("Heatmap_Un_shrunk.emf"), width = 10, height = 10 )
fig
dev.off()

# cols= p (number of genes) , rows are samples size
num_genes <- paste(rep("p", length(p)), p, sep="=") # rows
sample_size <- paste(rep("n", length(n)), n, sep="=") # columns
mc.ENF <- data.frame(genes = num_genes, 
                     t( diff.shrunk.pos), nrow = length(p), ncol =length(n) )

#mc.ENF$genes <- factor(mc.ENF$genes, levels = sort(unique(mc.ENF$genes)))

names(mc.ENF )[2:(length(n)+1)] <- sample_size 
mc.ENF<-mc.ENF[,-c(ncol(mc.ENF)-1, ncol(mc.ENF))]
df_heatmap.gnet <- melt(mc.ENF , id.vars = "genes")
names(df_heatmap.gnet)[2:3] <- c("sample", "P")
head(df_heatmap.gnet)
df_heatmap.gnet$genes <- factor(df_heatmap.gnet$genes, levels = unique(mc.ENF$genes))


fig = ggplot(df_heatmap.gnet, aes(x=sample, y=genes  )) +
  geom_tile(aes(fill = P), color = "white") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red" ,  midpoint = 0, limits=c(min.colour, max.colour))+
  geom_text(aes(label = round(df_heatmap.gnet$P, 3)),size=4) +
  ylab("") +
  xlab("") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size=14),
        axis.text.y = element_text(size=14) )+
  labs(fill = "|Shrunk - 0.5|")  


library("devEMF")
emf(file = paste0("Heatmap_Shrunk.emf"), width = 10, height = 10 )
fig
dev.off()



#.........................................
# Figure S4).- z-score and normality deviance 
#.........................................

# For small samples sizes (n<30) the standarized sample data 
# standarized data = (x-mean)/sd does not necessarily look normal
# it is common that the sd is poorly estimated and the z scores might show
# a different relationship
set.seed(123)
s.size = 100
s = rnorm(n=s.size, mean = 5, sd = 5)

 #  centering is done by subtracting the column means (omitting NAs) 
 # scaling is done by dividing the (centered) columns of x by their standard deviations
 mean(s) # far from 10?
 sd(s)# far from 100?
 
 norm.s = scale(s)
 cs = ecdf( x = norm.s)
 
 theoretical = pnorm(seq(-4, 4, by=0.1), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
 
 library("devEMF")
 emf(file = paste0("empirical_",s.size,".emf"), width = 5, height = 5 )
 # empirical function
 plot( cs, verticals = TRUE, do.points = TRUE , xlim=c(-4, 4), ylab="cumulative distribution")
 lines( seq(-4, 4, by=0.1) , theoretical ,"l", col="red", lwd=2)
 text(2, 0.5, paste0("n = ", s.size))
 dev.off()
 # predicted by empirical function
 plot( seq(-4, 4, by=0.1),cs(seq(-4, 4, by=0.1)), verticals = TRUE, do.points = TRUE)
 lines( seq(-4, 4, by=0.1) , theoretical, col="red" )
 
 # Remarks: The CLT holds the sample must be sufficiently large (n > 30). 
 # Two exceptions to this. If the population is normal, then the result holds for samples of any size 
 

 #..................................................................     
 # Figure S5:  ROC and PR curves 
 #..................................................................
 
 # ROC Curves show the trade-off between the true positive rate and false positive rate
 # varying the threshold.
 # Precision-Recall curves show the trade-off between the true positive rate (recall) and the positive predictive value (precision)
 # varying the threshold.
 # ROC curves are appropriate when the observations are balanced between each class, 
 # Precision-recall curves are appropriate for imbalanced datasets.
 #..................................................................
 rm(list=ls())
 
 
 library(ggplot2)
 library(igraph)
 library(GeneNet)
 library(stats4) 
 library(reshape)
 
 
 setwd("/home/victor/Desktop/")
 source("functions_Unshrink.R")
 
 
 set.seed(123)
 p<-10
 n1<- 50# c(1:10)*10
 etaA<- 0.3 # p=100 eta 0.005, p=30 eta 0.1
 
 TrueNET  <- ggm.simulate.pcor(p, etaA)#ggm.simulate.pcor.MODIF
 TrueNET
 positives.idx <- which(sm2vec(TrueNET)!=0)
 non.positives.idx <- which(sm2vec(TrueNET)==0)               
 
 # the simulated partial correlation
 plot(sm2vec(TrueNET)[positives.idx ] , pch=16 , ylab='partial correlations')
 text(x = 3, y = 0, labels = paste0('# pcors= ',length(positives.idx)) , cex = 1)
 
 
 cat('positives = ', sm2vec(TrueNET)[positives.idx])
 
 sim.data1 <-
   lapply(n1, function(x) {
     ggm.simulate.data(x , TrueNET)
   })
 
 GGM1.shrunk <-
   lapply(sim.data1, function(x) {
     sm2vec(pcor.shrink(x, verbose = FALSE))
   })
 
 opt.lambdas <-
   lapply(sim.data1, function(x) {
     return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
   })
 
 unshrunk1 <-
   lapply(sim.data1, function(x) {
     unshrink_GGM(x, l_0 = 0.01,PLOT = F , corrected = F)
   })
 
 
 simple_roc <- function(labels, scores){
   labels <- labels[order(scores, decreasing=TRUE)]
   data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
 }
 
 simple_PR <- function(labels, scores){
   labels <- labels[order(scores, decreasing=TRUE)]
   #TPR=cumsum(labels)/sum(labels)
   #FPR=cumsum(!labels)/sum(!labels)
   data.frame(TPR=cumsum(labels)/sum(labels), PPV=cumsum(labels)/sum(labels)/((cumsum(labels)/sum(labels))+(cumsum(!labels)/sum(!labels))), labels)
 }
 
 simple_auc <- function(TPR, FPR){
   # inputs already sorted, best scores first 
   dFPR <- c(diff(FPR), 0)
   dTPR <- c(diff(TPR), 0)
   sum(TPR * dFPR) + sum(dTPR * dFPR)/2
 }
 
 simp_roc1 <- simple_roc(c(sm2vec(TrueNET)!=0), abs(GGM1.shrunk[[1]]) )
 Usimp_roc1 <- simple_roc(c(sm2vec(TrueNET)!=0), abs(unshrunk1[[1]]) )
 simp_pr1 <- simple_PR(c(sm2vec(TrueNET)!=0), abs(GGM1.shrunk[[1]]) )
 Usimp_pr1 <- simple_PR(c(sm2vec(TrueNET)!=0), abs(unshrunk1[[1]]) )
 
 #library("devEMF")
 #emf(file = paste0("ROC_p",p, "n_", n1," etA_",etaA,".emf"), width = 4, height = 4 )
 plot(simp_roc1$FPR, simp_roc1$TPR, col="black", type = "b", cex=1.5, ylim=c(0, 1),xlim=c(0, 1),
      ylab="TPR", xlab="FPR")
 lines(Usimp_roc1$FPR, Usimp_roc1$TPR, col="blue", pch=20, type = "b")
 abline(a=0, b=1)          
 text(x = c(0.75),y=0.25, 
      labels = paste0('shrunk AUC=',signif(simple_auc(simp_roc1$TPR, simp_roc1$FPR), 3) ))
 text(x = c(0.75),y=0.15, 
      labels = paste0('un-shrunk AUC=',signif(simple_auc(Usimp_roc1$TPR, Usimp_roc1$FPR), 3) ))
 text(x = c(0.5),y=0.05, 
      labels = paste0('p=',p, " n=", n1," etA =",etaA) )   
 dev.off() 
 
 
 library("devEMF")
# emf(file = paste0("PR_p",p, "n_", n1," etA_",etaA,".emf"), width = 4, height = 4 )
 plot(simp_pr1$PPV, simp_pr1$TPR, col="black", type = "b", cex=1.5, ylim=c(0, 1), xlim=c(0, 1),
      ylab="PPV", xlab="TPR")
 lines(Usimp_pr1$PPV, Usimp_pr1$TPR, col="blue", pch=20, type = "b")
 abline(a=0, b=1)          
 text(x = c(0.5),y=0.25, 
      labels = paste0('shrunk AUC =',signif(simple_auc(simp_pr1$PPV, simp_pr1$TPR), 3) ))
 text(x = c(0.5),y=0.15, 
      labels = paste0('un-shrunk AUC =',signif(simple_auc(Usimp_pr1$PPV, Usimp_pr1$TPR), 3) ))   
 text(x = c(0.5),y=0.05, 
      labels = paste0('p=',p, " n=", n1," etA =",etaA) )   
 dev.off()   
 #f1.unshrunk = 2 * Usimp_pr1$PPV * Usimp_pr1$TPR/ (Usimp_pr1$PPV+ Usimp_pr1$TPR) 
 #f1.shrunk = 2 * simp_pr1$PPV * simp_pr1$TPR/ (simp_pr1$PPV+ simp_pr1$TPR) 
 #plot(f1.unshrunk , f1.shrunk)
 
 #..................................................................     
 # Supp Figure 3  : 
 #..................................................................
 rm(list=ls())
 
 
 library(ggplot2)
 library(igraph)
 library(GeneNet)
 library(stats4) 
 library(RColorBrewer)#library(pracma)
 library(reshape)
 
 
 setwd("/home/victor/Desktop/")
 source("/media/victor/VICTOR/Semester7(15032017)/RUG/Shrinkage/functions_pvalues.R")
 source("/media/victor/VICTOR/Semester7(15032017)/RUG/Shrinkage Symbolical/Rcode/functions_Unshrink.R")
 
 
 set.seed(123)
 
 
 p<-10
 n1<- c(10, 30, 1000) 
 etaA<- 0.1 # p=100 eta 0.005, p=30 eta 0.1
 TrueNET  <- ggm.simulate.pcor(p, etaA)#ggm.simulate.pcor.MODIF
 TrueNET
 positives.idx <- which(sm2vec(TrueNET)!=0)
 non.positives.idx <- which(sm2vec(TrueNET)==0)               
 sm2vec(TrueNET)[positives.idx]
 
 sim.data1 <-
   lapply(n1, function(x) {
     ggm.simulate.data(x , TrueNET)
   })
 
 GGM1.shrunk <-
   lapply(sim.data1, function(x) {
     sm2vec(pcor.shrink(x, verbose = FALSE))
   })
 
 opt.lambdas <-
   lapply(sim.data1, function(x) {
     return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
   })
 
 unshrunk1 <-
   lapply(sim.data1, function(x) {
     unshrink_GGM (x,  PLOT = F)
   })
 
 
 one = data.frame('Method'= c(rep('Shrunk', length(GGM1.shrunk[[1]][positives.idx ])), rep('Unshrunk', length(unshrunk1 [[1]][positives.idx ])),rep('Actual', length(sm2vec(TrueNET)[positives.idx ])  )) , 
                  'pcor' = c(GGM1.shrunk[[1]][positives.idx ] ,unshrunk1[[1]][positives.idx ], sm2vec(TrueNET)[positives.idx ]) ,
                  'pair'= c(rep(names(unshrunk1[[1]])[positives.idx ], 3) ))
 two = data.frame('Method'= c(rep('Shrunk', length(GGM1.shrunk[[2]][positives.idx ])), rep('Unshrunk', length(unshrunk1[[2]][positives.idx ])),rep('Actual', length(sm2vec(TrueNET)[positives.idx ])  )) , 
                  'pcor' = c(GGM1.shrunk[[2]][positives.idx ] ,unshrunk1[[2]][positives.idx ], sm2vec(TrueNET)[positives.idx ]),
                  'pair'= c(rep(names(unshrunk1[[2]])[positives.idx ], 3) ))
 three = data.frame('Method'= c(rep('Shrunk', length(GGM1.shrunk[[3]][positives.idx ])), rep('Unshrunk', length(unshrunk1[[3]][positives.idx ])),rep('Actual', length(sm2vec(TrueNET)[positives.idx ])  )) , 
                    'pcor' = c(GGM1.shrunk[[3]][positives.idx ] ,unshrunk1[[3]][positives.idx ], sm2vec(TrueNET)[positives.idx ]),
                    'pair'= c(rep(names(unshrunk1[[2]])[positives.idx ], 3) ))
 
 one$Method = factor(one$Method, levels=unique(one$Method))
 two$Method = factor(two$Method, levels=unique(two$Method))
 three$Method = factor(three$Method, levels=unique(three$Method))
 
 #  Parallel coordinates 
 ggplot(data= one,  aes(x= Method, y= pcor, color=Method)) + 
   geom_point(size= 3.5 ) +
   geom_line(data= one,  aes(x= Method, y= pcor, group = pair ), lwd=1)+
   theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
         axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
         axis.line = element_line(size = 1, colour = "grey80"),
         #axis.title =     element_blank(),
         axis.title.x = element_text(color = "grey20", size = 14),
         axis.title.y = element_text(color = "grey20", size = 14),
         panel.border = element_blank(),
         panel.background = element_blank(),
         panel.grid = element_line(colour = "grey"),
         legend.text=element_text(size=14),
         legend.title=element_text(size=14),
         legend.position = "top",
         legend.key = element_rect(fill = "white", colour = "black"))  +
   scale_color_manual(values= c( rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75) , rgb(0, 0, 0, 1) ))+
   scale_fill_manual(values= c( rgb(0, 0, 1, 0.5), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )) +
   scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+
   annotate(geom="text", x=1, y=0.85, label=paste0(expression("\u03BB")," = ", signif(opt.lambdas[[1]],2) ),
            color="black", size=5)+
   annotate(geom="text", x=1, y=0.75, label=paste0("n = ",n1[1]),
            color="black", size=5)+
   geom_hline(yintercept=0, linetype="dashed", color = "black")
 
 
 
 ggplot(data= two,  aes(x= Method, y= pcor, color=Method)) + 
   geom_point(size= 3.5 ) +
   geom_line(data= two,  aes(x= Method, y= pcor, group = pair ), lwd=1)  +
   theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
         axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
         axis.line = element_line(size = 1, colour = "grey80"),
         #axis.title =     element_blank(),
         axis.title.x = element_text(color = "grey20", size = 14),
         axis.title.y = element_text(color = "grey20", size = 14),
         panel.border = element_blank(),
         panel.background = element_blank(),
         panel.grid = element_line(colour = "grey"),
         legend.text=element_text(size=14),
         legend.title=element_text(size=14),
         legend.position = "top",
         legend.key = element_rect(fill = "white", colour = "black"))  +
   scale_color_manual(values= c( rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75) , rgb(0, 0, 0, 1) ))+
   scale_fill_manual(values= c( rgb(0, 0, 1, 0.5), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )) +
   scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+
   annotate(geom="text", x=1, y=0.85, label=paste0(expression("\u03BB")," = ", signif(opt.lambdas[[2]],2) ),
            color="black", size=5)+
   annotate(geom="text", x=1, y=0.75, label=paste0("n = ",n1[2]),
            color="black", size=5)+
   geom_hline(yintercept=0, linetype="dashed", color = "black")
 
 
 ggplot(data= three,  aes(x= Method, y= pcor, color=Method)) + 
   geom_point(size= 3.5 ) +
   geom_line(data= three,  aes(x= Method, y= pcor, group = pair ), lwd=1)  +
   theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
         axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
         axis.line = element_line(size = 1, colour = "grey80"),
         #axis.title =     element_blank(),
         axis.title.x = element_text(color = "grey20", size = 14),
         axis.title.y = element_text(color = "grey20", size = 14),
         panel.border = element_blank(),
         panel.background = element_blank(),
         panel.grid = element_line(colour = "grey"),
         legend.text=element_text(size=14),
         legend.title=element_text(size=14),
         legend.position = "top",
         legend.key = element_rect(fill = "white", colour = "black"))  +
   scale_color_manual(values= c( rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75) , rgb(0, 0, 0, 1) ))+
   scale_fill_manual(values= c( rgb(0, 0, 1, 0.5), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )) +
   scale_y_continuous(limits =c(-1, 1) ,breaks = round(seq(-1, 1, by = 0.2),1))+
   annotate(geom="text", x=1, y=0.85, label=paste0(expression("\u03BB")," = ", signif(opt.lambdas[[3]],2) ),
            color="black", size=5)+
   annotate(geom="text", x=1, y=0.75, label=paste0("n = ",n1[3]),
            color="black", size=5)+
   geom_hline(yintercept=0, linetype="dashed", color = "black")
 
 
 ## Session Info
 sessionInfo(package = NULL)
 
# ###################################################################
# # Supplementary material
# ###################################################################
# #........................................
# # Figure S1) and Supplementary Table S1
# #
# # Compare the order of the partial cor for different lambdas against the actual GGM
# # let's check the top 10% most strong ones
# #........................................
# 
# rm(list=ls())
# 
# 
# library(igraph)
# library(GeneNet)
# library(ggplot2)     
# library(stats4) 
# source("G:/Semester7(15032017)/RUG/Shrinkage/functions_pvalues.R")
# source("functions_Unshrink.R")
# 
# 
# set.seed(123)
# 
# p<- 10
# n<- 30
# 
# etaA<-0.5
# top= floor(0.25*0.5*p*(p-1))
# 
# lambda<-   c(1:10)/10 
# 
# TrueNET <- ggm.simulate.pcor(p, etaA)
# true.idx<- sort( abs(sm2vec(TrueNET)),decreasing = TRUE, index.return=TRUE)$ix
# 
# sim.data <- ggm.simulate.data( n , TrueNET)
# 
# opt.lambda<- attr(pcor.shrink(sim.data),"lambda")
# 
# V<-eigen(cor(sim.data), symmetric=TRUE)$vectors
# D<-eigen(cor(sim.data), symmetric=TRUE)$values
# 
# which(D==0)
# diag(1/((1-0)*D+0))
# 
# GGM <-lapply( lambda , function(L) { sm2vec(-cov2cor(V %*% diag(1/((1-L)*D+L)) %*%t(V)) )})
# idx<- lapply(GGM, function(x) { sort( abs( x),decreasing = TRUE, index.return=TRUE)$ix})
# 
# indices<-data.frame(true.idx, idx )
# colnames(indices)<-c("true" , lambda)
# 
# indices
# 
# #............. True GGM
# gtrue<-graph_from_adjacency_matrix(TrueNET!=0, mode = c("undirected"), diag=FALSE )
# plot(gtrue, layout = layout_in_circle(gtrue),edge.width= 20*abs(TrueNET[sm2vec(TrueNET)!=0]),
#      vertex.label.color="black", vertex.label.dist=3, 
#      vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	)
# text(-1.2, -1.2,  expression(True_GGM),
#      cex=1.5, pos=3, col="black") 
# 
# #..........top 25% edges for lambda GGM
# top10_1 = lapply(idx, function(x) { x[1:top]})
# nontop = lapply(idx, function(x) { x[-c(1:top)]})
# topGGM<- lapply(c(1:length(lambda)), function(x) { 
#   t<-GGM[[x]] 
#   t[ c(nontop[[x]]) ]=0
#   t<- vec2sm(t)
#   return (t)}
# )
# 
# sapply(c(1:length(lambda)), function(x) {
#   print(x)
#   g1<-graph_from_adjacency_matrix(topGGM[[x]]!=0, mode = c("undirected"), diag=FALSE )
#   #win.metafile(paste0(x,"top25.wmf"))
#   plot(g1, layout = layout_in_circle(g1),edge.width= 20*abs(topGGM[[x]][topGGM[[x]]!=0]),
#        vertex.label.color="black", vertex.label.dist=1.5, 
#        vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	)
#   text(-1.2, -1.2,  paste("lambda =", eval(lambda[x])),
#        cex=1.5, pos=3, col="black") 
#   #dev.off()
# })
# 
# nooverlap<-lapply(topGGM, function(x){ 
#   diag(x)=0
#   z<-x!=0
#   y<-TrueNET!=0
#   return (sum(sm2vec(x+y==1))/(0.5*p*(p-1)) ) })
# 
# 
# df<-data.frame(nooverlap)
# colnames(df)<-c("% different edges")
# df
# plot(lambda, df, ylim=c(0, 1))
# abline(v=opt.lambda, col="red")
# 
# 
# # confusion matrix edges   TABLE( Row, Cols)
# edges <-lapply( topGGM , function(x) { which(sm2vec(x!=0)) })
# confMatrix<-data.frame(sm2vec(TrueNET!=0), sm2vec(topGGM[[1]])!=0, sm2vec(topGGM[[2]])!=0)
# colnames(confMatrix)=c("Real","1","2")
# 
# # True vs GGM1   
# table(confMatrix$`1`, confMatrix$Real)  
# # True vs GGM2
# table(confMatrix$`2`, confMatrix$Real)
# # GGM1 vs GGM2
# table(confMatrix$`2`,confMatrix$`1`)
# 
# # confusion matrix connected genes TABLE( Row, Cols)
# connTrue<- unique( sm.index( diag(p))[which(sm2vec(TrueNET!=0))] )
# conn <-lapply( topGGM , function(x) { unique( sm.index( diag(p))[which(sm2vec(x!=0))] )})
# 
# # GGM1 vs GGM2   
# setdiff(conn[[1]],conn[[5]])
# 
# # change order 
# ttt<-matrix(unlist(GGM), 0.5*p*(p-1),10)
# ttt <- data.frame(t(ttt))
# rownames(ttt)<-lambda
# # expand with monotnic function
# #ttt<-  100*ttt
# 
# #win.metafile("3_Paths_ChangeOrder.wmf")
# plot.ts(ttt,
#         plot.type = c("single"), 
#         lw=2, lty = 1:2, col=rainbow(0.5*p*(p-1)),
#         ylab="partial correlation",xlab="shrinkage",
#         axes=F, cex=1)
# 
# axis(2)
# axis(1, labels=c(lambda), at=seq(from=1, by=1, length.out=length(lambda))) 
# box() 
# #dev.off()
# 
# #...................................
# # pcorr order change scatter 
# 
# id<-sm.index(matrix(0, p ,p), diag = FALSE)
# df<-data.frame(id[,1],id[,2])
# edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
# 
# scat.pcor<-data.frame(edge, GGM[[1]],GGM[[3]])
# colnames(scat.pcor)<-c("Edges","one","three")
# scat.pcor$Edges<-as.character(scat.pcor$Edges)
# 
# #win.metafile("scattertop25.wmf")
# ggplot(scat.pcor, aes(x=one, y=three))  +  geom_point(color = "red", size = 1) +
#   theme_classic(base_size = 12, base_family ="")+
#   geom_text(aes(label=scat.pcor$Edges),check_overlap = TRUE, size=5, nudge_y = 0.05, nudge_x = 0.01)+labs(x = "1",y = "2")  +
#   xlim(min(-1, -1 ,na.rm = TRUE), max(1, 1, na.rm = TRUE))+
#   ylim(min(-1, -1, na.rm = TRUE), max(1, 1, na.rm = TRUE))
# dev.off()
# 
# 
# #
# #Figure 4
# #...................................
# # p- values order change
# 
# alpha=0.05
# p.vals<-mapply(function(x, y){-log10(p.shrunk(x, p, n, y))},x=GGM[1:9],y=as.list(lambda[1:9])  ) 
# 
# id<-sm.index(matrix(0, p ,p), diag = FALSE)
# df<-data.frame(id[,1],id[,2])
# edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
# #
# scat.pval<-data.frame(edge, p.vals[,1],p.vals[,3])
# colnames(scat.pval)<-c("Edges","one","three")
# scat.pval$Edges<-as.character(scat.pval$Edges)
# scat.pval$Edges[ which(scat.pval$one < -log10(alpha) & scat.pval$three < -log10(alpha)) ]<-NA
# 
# #win.metafile("scattertop25.wmf")
# ggplot(scat.pval, aes(x=one, y=three))  +  geom_point(color = "red", size = 0.5) +
#   theme_classic(base_size = 12, base_family ="")+
#   geom_text(aes(label=scat.pval$Edges),size=5, nudge_y = 0.050, nudge_x = 0.05)+labs(x = "1",y = "2")+
#   geom_hline(yintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(alpha), linetype="dashed", color = "grey")
# dev.off()
# 
# #........................................
# # 4) SIM data: Comparing networks obtained with different sample sizes (diff optimal lambda)
# #........................................
# 
# rm(list=ls())
# 
# 
# library(igraph)
# library(GeneNet)
# library(stats4) 
# source("F:/Semester7(15032017)/RUG/Shrinkage/functions_pvalues.R")
# source("F:/Semester7(15032017)/RUG/Shrinkage Symbolical/Rcode/functions_Unshrink.R")
# 
# alpha=0.05
# 
# p<- 30
# n1<- 30
# n2<- floor(n1-0.25*n1)
# n2
# 
# etaA<- 0.5
# 
# top<-floor(0.1*0.5*p*(p-1))
# top
# 
# TrueNET <- ggm.simulate.pcor(p, etaA)
# true.idx<- sort( abs(sm2vec(TrueNET)),decreasing = TRUE, index.return=TRUE)$ix
# 
# sim.data1 <- ggm.simulate.data( n1 , TrueNET)
# sim.data2 <- ggm.simulate.data( n2, TrueNET)
# 
# # shrunk unshrunk
# 
# GGM1 <-pcor.shrink(sim.data1)
# GGM2 <-pcor.shrink(sim.data2)
# # ...............................
# # optimal lambda => not comparable (top, cohen significant)
# #.............................     
# 
# # cohen criteria 0.5
# plot(sm2vec(GGM1)/(1-attr(pcor.shrink(sim.data1),"lambda")),
#      sm2vec(GGM2)/(1-attr(pcor.shrink(sim.data2),"lambda")), xlim=c(-1, 1),
#      pch=20, ylim=c(-1, 1), xlab="shrunk 1",ylab="shrunk 2" )
# abline(v=c(-0.5, 0.5), col="grey",lw=2)
# abline(h=c(-0.5, 0.5), col="grey",lw=2)
# abline(a=0, b=1, col="red",lw=2)
# 
# 
# # Top criteria
# idx1<- sort( abs(sm2vec(GGM1)),decreasing = TRUE, index.return=TRUE)$ix
# idx2<- sort( abs(sm2vec(GGM2)),decreasing = TRUE, index.return=TRUE)$ix
# 
# indices<-data.frame(true.idx, idx1, idx2 )
# colnames(indices)<-c("true" , attr(pcor.shrink(sim.data1),"lambda"),attr(pcor.shrink(sim.data2),"lambda"))
# indices
# 
# 
# top10_1 = idx1[c(1:top)]
# top10_2 = idx2[1:top]
# nontop1 = idx1[-c(1:top)]
# nontop2 = idx2[-c(1:top)]
# 
# GGM1<-sm2vec(GGM1)
# GGM2<-sm2vec(GGM2)
# GGM1[nontop1]=0
# GGM2[nontop2]=0
# GGM1<-vec2sm(GGM1)
# GGM2<-vec2sm(GGM2)
# diag(GGM1)<-0
# diag(GGM2)<-0
# 
# 
# g1<-graph_from_adjacency_matrix(GGM1!=0, mode = c("undirected"), diag=FALSE )
# plot(g1, layout = layout_in_circle(g1),edge.width= 20*abs(GGM1[GGM1!=0]),
#      vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	,
#      xlab="shrunk 1")
# text(-1.2, -1.2,  paste("lambda =", 
#                         round(eval(attr(pcor.shrink(sim.data1),"lambda")),2)),
#      cex=1, pos=3, col="black") 
# 
# 
# g2<-graph_from_adjacency_matrix(GGM2!=0, mode = c("undirected"), diag=FALSE )
# plot(g2, layout = layout_in_circle(g2),edge.width= 20*abs(GGM2[GGM2!=0]),
#      vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c("black"), vertex.size=8 , edge.label.family="Times"	,
#      xlab="shrunk 2")
# text(-1.2, -1.2,  paste("lambda =", 
#                         round(eval(attr(pcor.shrink(sim.data2),"lambda")),2)),
#      cex=1, pos=3, col="black") 
# 
# #...............
# temp1<- GGM1!=0
# temp2<-  GGM2!=0
# nooverlap<-temp1+temp2==1
# 
# df<-data.frame(sum(sm2vec(nooverlap))/(0.5*p*(p-1)))
# colnames(df)<-c("% different edges")
# df
# 
# g3<-graph_from_adjacency_matrix(nooverlap, mode = c("undirected"), diag=FALSE )
# plot(g2, layout = layout_in_circle(g2),edge.width= 2*abs(nooverlap[nooverlap!=0]),
#      vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c("black"), vertex.size=8 , edge.label.family="Times",
#      xlab="different edges")
# 
# 
# 
# # ...............................
# # optimal lambda => not comparable (TOP)
# #.............................                    
# # confusion matrix edges   TABLE( Row, Cols)
# confMatrix<-data.frame(sm2vec(TrueNET!=0), sm2vec(GGM1!=0), sm2vec(GGM2!=0))
# colnames(confMatrix)=c("Real","1","2")
# 
# # True vs GGM1   
# table(confMatrix$`1`, confMatrix$Real)  
# # True vs GGM2
# table(confMatrix$`2`, confMatrix$Real)
# # GGM1 vs GGM2
# table(confMatrix$`2`,confMatrix$`1`)
# 
# # confusion matrix connected genes TABLE( Row, Cols)
# connTrue<- unique( sm.index( diag(p))[which(sm2vec(TrueNET!=0))] )
# conn1<- unique( sm.index( diag(p))[which(sm2vec(GGM1!=0))] )
# conn2<- unique( sm.index( diag(p))[which(sm2vec(GGM2!=0))] )
# 
# # GGM1 vs GGM2   
# setdiff(conn1, conn2) 
# setdiff(which(sm2vec(GGM1!=0)),which(sm2vec(GGM2!=0)))
# #..............................................
# 
# # Comparison Optimals vs TRUE
# idx.TRUE<- sort( abs( sm2vec(TrueNET)),decreasing = TRUE, index.return=TRUE)$ix
# 
# cumul.diff.shrunk1<-sapply(c(1:length(idx1)),function(x){
#   length(setdiff(idx1[1:x],idx.TRUE[1:x]))
# })
# cumul.diff.shrunk2<-sapply(c(1:length(idx2)),function(x){
#   length(setdiff(idx2[1:x],idx.TRUE[1:x]))
# })
# 
# plot(cumul.diff.shrunk1-cumul.diff.shrunk2,
#      c(1:length(idx1))/c(length(idx1)),col="red",pch=20,
#      xlim=c(-10, 10),xlab=c("number diff edges GGM1 - GGM2 wrt true"), ylab=c("top edges %"))
# abline(v=0)
# 
# # p values shrunks
# 
# p.shrunk1<-p.shrunk(sm2vec(pcor.shrink(sim.data1)),p, n1, lambda=attr(pcor.shrink(sim.data1),"lambda"))
# p.shrunk2<-p.shrunk(sm2vec(pcor.shrink(sim.data2)),p, n2, lambda=attr(pcor.shrink(sim.data2),"lambda"))
# id<-sm.index(matrix(0, p ,p), diag = FALSE)
# df<-data.frame(id[,1],id[,2])
# edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
# scat.pval<-data.frame(edge,-log10(p.shrunk1),-log10(p.shrunk2))
# colnames(scat.pval)<-c("Edges","one","three")
# scat.pval$Edges<-as.character(scat.pval$Edges)
# scat.pval$Edges[ which(scat.pval$one < -log10(alpha) & scat.pval$three < -log10(alpha)) ]<-NA
# 
# #win.metafile("scattertop25.wmf")
# ggplot(scat.pval, aes(x=one, y=three))  +  geom_point(color = "red", size = 0.5) +
#   theme_classic(base_size = 12, base_family ="")+
#   geom_text(aes(label=scat.pval$Edges),size=5, nudge_y = 0.050, nudge_x = 0.05)+labs(x = "1",y = "2")+
#   geom_hline(yintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(alpha), linetype="dashed", color = "grey")
# #dev.off()
# sum(p.shrunk1<=alpha)/length(p.shrunk1)
# sum(p.shrunk2<=alpha)/length(p.shrunk2)
# 
# 

#..................................................................     
# Supp Figure 3. Unshrunk are comparable               
#..................................................................  
library(plotly)

N <- DF$n
P <- DF$p
C <- DF$correction
#Graphing your 3d scatterplot using plotly's scatter3d type:

#plot_ly(x=X, y=pressure, z=dtime, type="scatter3d", mode="markers", color=pressure )

plot_ly(x=N, y=P, z=C, type="scatter3d", mode="markers", color=N)
#......................................
N <- as.factor(sort(rep(log(n), length(n))))
P <- as.factor(rep(log(p), length(p)))
C <- log(DF$correction)

plot_ly(x=N, y=P, z=C, type="scatter3d", mode="markers", color=N)

plot( log(p) , log(DF$correction)[1:20],type="l")
lines(log(p), log(n[1])-log(p),col="red")

#.............................................................
# ROC
p<-100
n1<-10
TrueNET  <- ggm.simulate.pcor.MODIF(p, 20)
sim.data1 <- ggm.simulate.data( n1, TrueNET)

GGM1 <-sm2vec(pcor.shrink(sim.data1, verbose = FALSE))
unshrunk1<-unshrink_GGM(sim.data1)
max(abs(unlist(unshrunk1)))

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

simp_roc1 <- simple_roc(c(sm2vec(TrueNET!=0)), abs(GGM1))
Usimp_roc1 <- simple_roc(c(sm2vec(TrueNET!=0)), abs(unshrunk1))


plot(simp_roc1$FPR, simp_roc1$TPR, col="black")
points(Usimp_roc1$FPR, Usimp_roc1$TPR, col="blue", pch=20)
abline(a=0, b=1)

library(pROC)
roc_obj.shrunk <- roc(c(sm2vec(TrueNET!=0)), abs(GGM1))
auc(roc_obj.shrunk)
roc_obj.unshrunk <- roc(c(sm2vec(TrueNET!=0)), abs(unshrunk1))
auc(roc_obj.unshrunk)
# 
# 
# 
# library(corrplot)
# 
# corrplot(L2.performance, method="number", 
#          is.corr = FALSE, addgrid.col = "black", order = "original")
# 
# corrplot(non.L2.performance, method="number", 
#          is.corr = FALSE, addgrid.col = "black", order = "original")
# 
# corrplot(pval.performance, method="number", 
#          is.corr = FALSE, addgrid.col = "black", order = "original")
# 
# corrplot(non.pval.performance, method="number", 
#          is.corr = FALSE, addgrid.col = "black", order = "original")
# differnt sample size vs true
# dev.off()
# plot(corr.unshrunk1, corr.unshrunk2, col="red",xlim=c(-1, 1),ylim=c(-1, 1),ylab="2",xlab="1")
# points(corr.shrunk1, corr.shrunk2, col="black",xlim=c(-1, 1),ylim=c(-1, 1))
# abline(a=0, b=1) # better= near diagonal
# 
# data.frame(c(unlist(corr.shrunk1)-c(unlist(corr.shrunk2))),
#            c(unlist(corr.unshrunk1)-unlist(corr.unshrunk2)))
# 
# boxplot(c(unlist(corr.shrunk1)),c(unlist(corr.shrunk2)),
#         c(unlist(corr.unshrunk1)),c(unlist(corr.unshrunk2)),col=c("white","white","red","red"))  
# 
#boxplot(sm2vec(TrueNET),c(unlist(GGM1)),c(unlist(GGM2)),
#         c(unlist(unshrunk1)),c(unlist(unshrunk2)))  


# differnt sample size across methods
# corr.UNvs<-mapply(function(X, Y){cor(X[positives.idx],Y[positives.idx])},Y=unshrunk2, X=unshrunk1)
# corr.vs<-mapply(function(X, Y){cor(X[positives.idx],Y[positives.idx])},Y=GGM2, X=GGM1)
# 
# plot(corr.vs, corr.UNvs, col="red",xlim=c(-1, 1),ylim=c(-1, 1))
# abline(a=0, b=1) # better= away diagonal
# 
# data.frame(corr.vs, corr.UNvs)
# hist(corr.UNvs, 5, col="red")
# hist(corr.vs, add=T, col="blue")
# 
#
# pval.TRUE<-p.shrunk(sm2vec(TrueNET),p, n1, lambda=0.0000001)
# pval.UNshrunk1<-lapply(unshrunk1, function(x){p.shrunk(x, p, n1, lambda=0.00000001)})
# pval.UNshrunk2<-lapply(unshrunk2, function(x){p.shrunk(x, p,(2)*n1, lambda=0.00000001)})
# 
# plot(pval.UNshrunk1[[1]], pval.UNshrunk2[[1]],col="red")
# data.frame(pval.UNshrunk1[[1]], pval.UNshrunk2[[1]])
# 
# plot(pval.TRUE, pval.UNshrunk2[[1]],col="black")
# points(pval.TRUE, pval.UNshrunk2[[1]],col="red")

# pval.shrunk1<-mapply(function(X, Y){p.shrunk(X, p, n1, lambda=Y)},
#                      X=GGM1, Y=opt.lambdas)

# plot(pval.TRUE[positives.idx], pval.shrunk1[,1][positives.idx])
# points(pval.TRUE[positives.idx], pval.UNshrunk2[[1]][positives.idx],col="red")

#GGM1 <-lapply( sim.data1, function(x){ sm2vec(pcor.shrink(x))})
#opt.lambdas <-lapply( sim.data1, function(x){ return(attr(pcor.shrink(x),"lambda"))})


#pval.TRUE<-p.shrunk(sm2vec(TrueNET),p, n1, lambda=0.0000001)

#lapply(GGM1, opt.lambdas, function(x, y){p.shrunk(x, p, n1, lambda=y)})

# plot(pval.TRUE[positives.idx], pval.shrunk1[,1][positives.idx])
# points(pval.TRUE[positives.idx], pval.UNshrunk1[[1]][positives.idx],col="red")

# rank correlation of magnitudes
#  corr.shrunk1<-lapply(GGM1, function(x){ cor(x[positives.idx]^2, sm2vec(TrueNET)[positives.idx]^2, method = "spearman")})
#  corr.unshrunk1<-lapply(unshrunk1, function(x){ cor(x[positives.idx]^2, sm2vec(TrueNET)[positives.idx]^2, method = "spearman")})
#  
#      plot(corr.shrunk1, corr.unshrunk1, col="red",xlim=c(0, 1),ylim=c(0, 1))
#      abline(a=0, b=1)
#      
# # rank correlation of pvals
#      corr.shrunk1<-apply(pval.shrunk1, 2, function(x){ cor(x[positives.idx],pval.TRUE[positives.idx],method = "spearman")})
#      corr.unshrunk1<-lapply(pval.UNshrunk1, function(x){ cor(x[positives.idx],pval.TRUE[positives.idx],method = "spearman")})
#      
#      plot(unlist(corr.shrunk1), unlist(corr.unshrunk1),col="red",xlim=c(0, 1),ylim=c(0, 1))
#      abline(a=0, b=1)
#      
#      hist(unlist(corr.unshrunk1)-corr.shrunk1, col="red",10)

# 
# t.UN<- t.test(matrix(unlist(corr.unshrunk1), nrow=length(corr.unshrunk1), byrow=T) ,
#         matrix(unlist(corr.shrunk1), nrow=length(corr.shrunk1), byrow=T), 
#         alternative = "greater")
#  
#  t.UN$p.value<0.05



#                 #..........................................
#                 level=0.005
#                 # confusion matrix edges   TABLE( Row, Cols)
#                 confMatrix<-data.frame(sm2vec(TrueNET!=0), 
#                                        abs(unshrunk1)>= level, abs(unshrunk2)>= level)
#                 colnames(confMatrix)=c("Real","1","2")
#                 # GGM1 vs GGM2
#                 table(confMatrix$`2`,confMatrix$`1`)
#                 
#                 
#                 # confusion matrix edges   TABLE( Row, Cols)
#                 confMatrix2<-data.frame(sm2vec(TrueNET!=0), 
#                                        abs(sm2vec(GGM1))>= level, abs(sm2vec(GGM2))>= level)
#                 colnames(confMatrix2)=c("Real","1","2")
#                 # GGM1 vs GGM2
#                 table(confMatrix2$`2`,confMatrix2$`1`)
#                 #..........................................
#                 
# 
# 
#                   # compare coefficients
#                   plot(unshrunk1, unshrunk2, ylim=c(-1, 1),
#                        pch=20, xlim= c(-1, 1) )
#                   
#                   points(sm2vec(GGM1)/(1-attr(pcor.shrink(sim.data1),"lambda")),
#                          sm2vec(GGM2)/(1-attr(pcor.shrink(sim.data2),"lambda")),col=rgb(.8, 0,.8, 0.2))
#                   abline(v=c(-0.5, 0.5), col="grey",lw=2)
#                   abline(h=c(-0.5, 0.5), col="grey",lw=2)
#                   abline(a=0, b=1, col="red",lw=2)
#                   
# # compare pvalues
# p.UNshrunk1<-p.shrunk(unshrunk1, p, n1, lambda=0.01)
# p.UNshrunk2<-p.shrunk(unshrunk2, p, n2, lambda=0.01)
# p.cohen.mid<-p.shrunk(0.3, p, n2, lambda=0.01)
# p.cohen.low<-p.shrunk(0.1, p, n2, lambda=0.01)
# 
# id<-sm.index(matrix(0, p ,p), diag = FALSE)
# df<-data.frame(id[,1],id[,2])
# edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
# scat.pval<-data.frame(edge,-log10(p.UNshrunk1),-log10(p.UNshrunk2))
# colnames(scat.pval)<-c("Edges","one","three")
# scat.pval$Edges<-as.character(scat.pval$Edges)
# scat.pval$Edges[ which(scat.pval$one < -log10(alpha) & scat.pval$three < -log10(alpha)) ]<-NA
# 
# #win.metafile("scattertop25.wmf")
# colref=c("black","red")
# ggplot(scat.pval, aes(x=one, y=three),
#        ylim=c(0, max(scat.pval)), xlim=c(0, max(scat.pval)))  +  
#   geom_point(color = colref[1+is.na(scat.pval$Edges)], size = (2-is.na(scat.pval$Edges))) +
#   theme_classic(base_size = 12, base_family ="")+
#   #geom_text(aes(label=scat.pval$Edges),size=5, nudge_y = 0.050, nudge_x = 0.05)+
#   labs(x = "1",y = "2")+
#   geom_hline(yintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_hline(yintercept=-log10(p.cohen.mid), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(p.cohen.mid), linetype="dashed", color = "grey")+
#   geom_hline(yintercept=-log10(p.cohen.low), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(p.cohen.low), linetype="dashed", color = "grey")
# #+
#  # geom_text(sum(p.UNshrunk1<=alpha)/length(p.UNshrunk1) )
#   #geom_text(sum(p.UNshrunk2<=alpha)/length(p.UNshrunk2)  )
# #...................................................
# # compare pvalues
# p.shrunk1<-p.shrunk(sm2vec(abs(GGM1)),p, n1, lambda=attr(pcor.shrink(sim.data1),"lambda"))
# p.shrunk2<-p.shrunk(sm2vec(abs(GGM2)),p, n2, lambda=attr(pcor.shrink(sim.data2),"lambda"))
# p.cohen.mid<-p.shrunk(0.3, p, n2, lambda=0.01)
# p.cohen.low<-p.shrunk(0.1, p, n2, lambda=0.01)
# 
# id<-sm.index(matrix(0, p ,p), diag = FALSE)
# df<-data.frame(id[,1],id[,2])
# edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
# scat.pval<-data.frame(edge,-log10(p.shrunk1),-log10(p.shrunk2))
# colnames(scat.pval)<-c("Edges","one","three")
# scat.pval$Edges<-as.character(scat.pval$Edges)
# scat.pval$Edges[ which(scat.pval$one < -log10(alpha) & scat.pval$three < -log10(alpha)) ]<-NA
# 
# colref=c("black","red")
# ggplot(scat.pval, aes(x=one, y=three),
#        ylim=c(0, max(scat.pval)), xlim=c(0, max(scat.pval)))  +  
#   geom_point(color = colref[1+is.na(scat.pval$Edges)], size = (2-is.na(scat.pval$Edges))) +
#   theme_classic(base_size = 12, base_family ="")+
#   #geom_text(aes(label=scat.pval$Edges),size=5, nudge_y = 0.050, nudge_x = 0.05)+
#   labs(x = "1",y = "2")+
#   geom_hline(yintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(alpha), linetype="dashed", color = "grey")+
#   geom_hline(yintercept=-log10(p.cohen.mid), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(p.cohen.mid), linetype="dashed", color = "grey")+
#   geom_hline(yintercept=-log10(p.cohen.low), linetype="dashed", color = "grey")+
#   geom_vline(xintercept=-log10(p.cohen.low), linetype="dashed", color = "grey")
# # geom_text(sum(p.UNshrunk1<=alpha)/length(p.UNshrunk1) )
# #geom_text(sum(p.UNshrunk2<=alpha)/length(p.UNshrunk2)  )

#............................................................
# #................................................................
# # P values
# # Network comparison for the actual positives and actual negatives
# #..................................................................
# rm(list=ls())
#                
#                
# library(ggplot2)
# library(igraph)
# library(GeneNet)
# library(stats4) 
# library(pracma)
# library(reshape)
#                
# setwd("/home/victor/Desktop/")
# source("/media/victor/VICTOR/Semester7(15032017)/RUG/Shrinkage/functions_pvalues.R")
# source("/home/victor/Desktop/functions_Unshrink.R")
#                
#                
# set.seed(123)
#                
# # choose etaA s.t. TrueNET has pcors greater than 0.5
# p<- 100
# n1<- c(1)*10
# etaA<- 0.01 
# TrueNET_2  <- ggm.simulate.pcor(p, etaA)#ggm.simulate.pcor.MODIF
# TrueNET <- matrix(0, ncol = p, nrow=p )
# #true.idx <- sort( abs(sm2vec(TrueNET)),decreasing = TRUE, index.return=TRUE)$ix
# positives.idx_2 <- which(TrueNET_2!=0, arr.ind = TRUE)
# TrueNET_2[positives.idx_2]
# 
# TrueNET[positives.idx_2[,1] , positives.idx_2[,2]] <- 
#               TrueNET_2[positives.idx_2[,1] , positives.idx_2[,2]]
# det(TrueNET)
# TrueNET
# positives.idx <- which(sm2vec(TrueNET)!=0)
# non.positives.idx <- which(sm2vec(TrueNET)==0)               
# #
# shrunk.pos <- c()
# unshrunk.pos <- c()
# shrunk.neg <- c()
# unshrunk.neg <- c()
# shrunk.pos.se <- c()
# unshrunk.pos.se <- c()
# shrunk.neg.se <- c()
# unshrunk.neg.se <- c()
# pval.shrunk<- c()
# pval.unshrunk<- c()
# opt.lambdas.mean<-c()
#            
# times = 3
#                
# for (samples in 1:length(n1)) {
#     cat("p=", p , "n=", n1[samples], "\n")
#                  
#    sim.data1 <-
#      lapply(c(1:times), function(x) {
#            ggm.simulate.data(n1[samples] , TrueNET)
#             })
#                  
#                  GGM1.shrunk <-
#                    lapply(sim.data1, function(x) {
#                      sm2vec(pcor.shrink(x, verbose = FALSE))
#                    })
#                  
#                  opt.lambdas <-
#                    lapply(sim.data1, function(x) {
#                      return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
#                    })
#                  
#                  
#                  while (sum(unlist(opt.lambdas) == 1)) {
#                        sim.data1 <-
#                          lapply(c(1:times), function(x) {
#                            ggm.simulate.data(n1[samples] , TrueNET)
#                          })
#                        
#                        GGM1.shrunk <-
#                          lapply(sim.data1, function(x) {
#                            sm2vec(pcor.shrink(x, verbose = FALSE))
#                          })
#                        opt.lambdas <-
#                          lapply(sim.data1, function(x) {
#                            return(attr(pcor.shrink(x, verbose = FALSE), "lambda"))
#                          })
#                        
#                  }
#                  
#                  unshrunk1 <-
#                    lapply(sim.data1, function(x) {
#                      unshrink_GGM (x, tol = 5, PLOT = F)
#                    })
#                  max(abs(unlist(unshrunk1)))
#                  
#                  
#                  # minimum lambda
#                  minL <- lapply(sim.data1, function(x) {
#                    D <- eigen(cor(x), symmetric = TRUE)$values
#                    tol = 1e-2
#                    minL = tol + (min(D)) / (min(D) - 1)
#                    
#                    return (minL)
#                  })
#                  
#                  
#                    pval.TRUE<-p.shrunk(sm2vec(TrueNET)[positives.idx],p, n1[samples],lambda=0)
#                  # #
#                  # # pval.UNshrunk<-lapply(unshrunk1, function(x){p.shrunk(x, p, n1[samples],lambda = 1e-5  )})
#                  # #
#                    pval.UNshrunk_1<-mapply( function(x, y){p.shrunk(x, p, n1[samples],lambda=y)} ,
#                                           x=unshrunk1, y=minL)
#                  # #
#                   pval.shrunk_1<-mapply( function(x, y){p.shrunk(x, p, n1[samples],lambda=y)} ,
#                                         x=GGM1[positives.idx],y=opt.lambdas)
#                  #
#                  #
#                  
#                   pval.UNshrunk<-split(pval.UNshrunk, rep(1:ncol(pval.UNshrunk), each = nrow(pval.UNshrunk)))
#                   pval.shrunk<-split(pval.shrunk, rep(1:ncol(pval.shrunk), each = nrow(pval.shrunk)))
#                  
#                  # new mean
#                  
#                  shrunk.pos <-
#                    cbind(shrunk.pos, rowMeans(sapply(pval.shrunk, function(x) {
#                      x[positives.idx]
#                    })))
#                  unshrunk.pos <-
#                    cbind(unshrunk.pos, rowMeans(sapply(pval.UNshrunk, function(x) {
#                      x[positives.idx]
#                    })))
#                  
#                  shrunk.neg <-
#                    cbind(shrunk.neg, rowMeans(sapply(pval.shrunk, function(x) {
#                      x[non.positives.idx]
#                    })))
#                  shrunk.neg <-
#                    cbind(shrunk.neg, rowMeans(sapply(pval.UNshrunk, function(x) {
#                      x[non.positives.idx]
#                    })))
# 
#                  opt.lambdas.mean <-
#                    cbind(opt.lambdas.mean, mean(unlist(opt.lambdas)))
#                  
#                  # new se
#                  shrunk.pos.se <- cbind(shrunk.pos.se, apply(sapply(pval.shrunk, function(x) {
#                    x[positives.idx]
#                  }) , 1 ,
#                  function(x) {
#                    sd(x) / sqrt(length(positives.idx))
#                  }))
#                  unshrunk.pos.se <- cbind(unshrunk.pos.se, apply(sapply(pval.UNshrunk, function(x) {
#                    x[positives.idx]
#                  }) , 1 ,
#                  function(x) {
#                    sd(x) / sqrt(length(positives.idx))
#                  }))
#                  shrunk.neg.se <- cbind(shrunk.neg.se, apply(sapply(pval.shrunk, function(x) {
#                    x[non.positives.idx]
#                  }) , 1 ,
#                  function(x) {
#                    sd(x) / sqrt(length(non.positives.idx))
#                  }))
#                  unshrunk.neg.se <- cbind(unshrunk.neg.se, apply(sapply(pval.UNshrunk, function(x) {
#                    x[non.positives.idx]
#                  }) , 1 ,
#                  function(x) {
#                    sd(x) / sqrt(length(non.positives.idx))
#                  }))
#                  # pval.shrunk<- cbind(pval.shrunk, rowMeans(pval.shrunk_1) )
#                  # pval.unshrunk<-cbind(pval.unshrunk, rowMeans(pval.UNshrunk_1) )
#                  #
# }
# ## distance plots
#                
# library(ggplot2)
# library(reshape2)
# #
# #dev.off()
# 
# 
# 
# #for(w in c(1:length(n1))){
# w=1
# n1[w]
# 
# spider.top<- c(rep(c(as.character(signif(sm2vec(TrueNET)[positives.idx],digits = 3) ), "0")
#                     ,1))
# df = data.frame(Group = spider.top,
#                 "Shrunk" = c(shrunk.pos[,w], mean(shrunk.neg[,w])) ,
#                 "Un-Shrunk" = c(unshrunk.pos[,w],mean(unshrunk.neg[,w])) ,
#                 "Actual" =c(sm2vec(TrueNET)[positives.idx],0) )
# 
# df.se = data.frame(Group = spider.top,
#                    "Shrunk" = c(shrunk.pos.se[,w], mean(shrunk.neg.se[,w])) ,
#                    "Un-Shrunk" = c(unshrunk.pos.se[,w],mean(unshrunk.neg.se[,w])) ,
#                    "Actual" =rep(0, length(c(sm2vec(TrueNET)[positives.idx],0)) ))
# 
# df.m <- melt(df, 
#              id.vars= c("Group"), 
#              measure.vars= c("Shrunk", "Un.Shrunk","Actual"),
#              variable.name= "Color",
#              value.name=    "val"
# )
# df.mse <- melt(df.se, 
#                id.vars= c("Group"), 
#                measure.vars= c("Shrunk", "Un.Shrunk","Actual"),
#                variable.name= "Color",
#                value.name=    "val"
# )
# 
# 
# df.m<-df.m[order(as.numeric(as.character(df.m$Group))), ]
# df.m$Group<-as.factor(df.m$Group)
# df.mse<-df.mse[order(as.numeric(as.character(df.mse$Group))), ]
# df.mse$Group<-as.factor(df.mse$Group)
# 
# df_fin<-cbind(df.m, df.mse$val)
# colnames(df_fin)<-c("Label",  "Method","ave","SE")
# df_fin$ave<- df_fin$ave
# df_fin$SE<-2*df_fin$SE
# 
#  # make V1 an ordered factor
# df_fin$Label <- factor(df_fin$Label, levels = unique(df_fin$Label))
# 
# fig<-ggplot(data=df_fin,  aes(x= Label, y=ave, group= Method, colour=Method, fill=NA)) + 
#   geom_point(size=2.5 ) +
#   geom_errorbar(aes(ymin= ave - SE, ymax= ave + SE , color=Method), width=0.1, size=1.5)+
#   ylim(-1, 1) + 
#   #scale_x_discrete(breaks = round(seq(-1, 1, by = 0.01),4)) +
#   theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
#         axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#         axis.line = element_line(size = 2, colour = "grey80"),
#         axis.title =     element_blank(),
#         #panel.grid.major = element_blank(),
#         #panel.grid.minor = element_blank(colour = "black"),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_line(colour = "grey"),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         legend.position = "top",
#         legend.key = element_rect(fill = "white", colour = "black"))  +
#         scale_color_manual(values= c( rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75) , rgb(0, 0, 0, 1) ))+
#         scale_fill_manual(values= c( rgb(0, 0, 1, 0.5), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )) 
# #+ coord_polar()
# 
# #tiff(paste0(paste0("p15_n_",n1[w]),".tiff"),width = 600, height = 450, pointsize = 22 )
# fig+ annotate(geom="text", x=2.5, y=0.85, label=paste0("lambda = ", signif(opt.lambdas.mean[w],2) ),
#                  color="black", size=5)+
#         annotate(geom="text", x=2.5, y=0.75, label=paste0("n = ",n1[w]),
#                  color="black", size=5)+
#    geom_hline(yintercept=0, linetype="dashed", color = "black")


#dev.off()
#}


# # ## Spider
# library(fmsb)
# par(mar=rep(0.8, 4))
# par(mfrow=c(2, 3))
# w=1
# for(w in c(1:length(n1))){
#     n1[w]
# 
# spider<- rbind(c(shrunk.pos[,w], mean(shrunk.neg[,w])),
#                c(unshrunk.pos[,w],mean(unshrunk.neg[,w])),
#                c(sm2vec(TrueNET)[positives.idx],0) )
# spider<-abs(spider)
# spider<-as.data.frame(spider)
# colnames(spider) <- c(as.character(signif(abs(sm2vec(TrueNET)[positives.idx]),digits = 3)) , "0")
# rownames(spider) <- c("Shrunk", "Un-shrunk","Actual")
# spider<-spider[,order(spider[3,])]
# 
# spider <- rbind(rep(1, ncol(spider)) , rep(0, ncol(spider)) , spider)
# 
# 
# colors_border=c( rgb(0.2, 0.5, 0.5, 0.9), rgb(0.8, 0.2, 0.5, 0.9) , rgb(0.7, 0.5, 0.1, 0.9) )
# colors_in=c( rgb(0.5, 0.5, 0.5, 0.1), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )
# 
# radarchart(spider,
#            #custom polygon
#            pcol=colors_border , pfcol=colors_in ,
#            pty=16, plwd =0.1,
#            #custom the grid
#            vlcex=1.75, axistype=1,
#            caxislabels=seq(0, 20, 5)/20, cglcol="grey",
#            axislabcol="black",  cglwd=2, cglty=1)
# legend(x=1.2, y=0.95, legend = rownames(spider[-c(1, 2),]), bty = "n", pch=20 , col=c(rgb(0.2, 0.5, 0.5, 0.9), rgb(0.8, 0.2, 0.5, 0.9) , rgb(0.7, 0.5, 0.1, 0.9) ), text.col = "black", cex=1.25, pt.cex=2)
# text(x=-.85, y= 1.25, labels = paste('lambda=',10), cex=1.5)
# text(x=-.85, y=-1.25, labels = paste('n=',n1[w]), cex=1.5)
# 
#     }






#..... PCA on cor

# The eigenvectors represent the principal components of S.
# The eigenvalues of S are used to find the proportion of the total variance explained by the components.
cumsum(D) / sum(D)
plot(D, pch = 20, ylab = "eigenvalues")
df <- data.frame(V[, 1], V[, 2])
colnames(df) <- c("PC1", "PC2")
plot(PC2 ~ PC1, data = df, cex = 0.25)
with(df, text(PC2 ~ PC1, labels = names), pos = 2)

cbind(V[, 1], V[, 2]) > 0.1

plot(
  unshrunk.pcor,
  opt.pcor,
  ylim = c(-1, 1),
  pch = 20,
  xlim = c(-1, 1)
)
abline(v = c(-0.5, 0.5),
       col = "grey",
       lw = 3)
abline(h = c(-0.5, 0.5),
       col = "grey",
       lw = 3)
abline(v = c(-0.3, 0.3),
       col = "light grey",
       lw = 2)
abline(h = c(-0.3, 0.3),
       col = "light  grey",
       lw = 2)
abline(a = 0,
       b = 1,
       col = "red",
       lw = 2)
# high corr, weak pcorr
hist(sm2vec(cor(ecoli)), main="histogram of correlations")

hist(unshrunk.pcor, main="histogram of partial correlations",
     col = rgb(0.1, 0.1, 0.1, 0.5), breaks = 15, xlab="pcor")
hist(opt.pcor, add = T, col = rgb(1, 1, 1, 0.5), breaks = 15)

# corr vs pcor
plot(
  sm2vec(pcor.shrink(ecoli)),
  sm2vec(cor(ecoli)),
  xlim = c(-(1), (1)),
  ylim = c(-1, 1),
  pch = 20,
  ylab = c("corr"),
  xlab = c("shrink pcorr")
)
points(unshrunk.pcor, sm2vec(cor(ecoli)), col = rgb(1, 0, 0, 0.2))
abline(v = 0.5 * c(-1, 1),
       col = "grey",
       lty = "dashed")
abline(h = 0.5 * c(-1, 1),
       col = "grey",
       lty = "dashed")

# # PCA on pcor
# GGM1 <- vec2sm(unshrunk.pcor)
# diag(GGM1) <- 1
# V2 <- eigen(GGM1, symmetric = TRUE)$vectors
# D2 <- eigen(GGM1, symmetric = TRUE)$values
# cumsum(D2) / sum(D2)
# plot(D2, pch = 20, ylab = "eigenvalues")
# df <- data.frame(V2[, 1], V2[, 2])
# colnames(df) <- c("PC1", "PC2")
# plot(PC2 ~ PC1, data = df, cex = 0.25)
# with(df, text(PC2 ~ PC1, labels = names), pos = 2)
# cbind(V2[, 1], V2[, 2]) > 0.5

#... TOP Criteria
idx.unshrunk <-
  sort(abs(unshrunk.pcor),
       decreasing = TRUE,
       index.return = TRUE)$ix
idx.shrunk <-
  sort(abs(opt.pcor),
       decreasing = TRUE,
       index.return = TRUE)$ix

counts <-
  table(idx.unshrunk[1:top] %in% idx.shrunk[1:top])
rownames(counts) <- c("both not top", "both top")
barplot(
  counts,
  xlab = "?",
  col = c("darkblue", "red"),
  legend = rownames(counts)
)

#top.both<-names[ unique(sm.index(diag(p))[setdiff(which(cohen.unshrunk==T),
#                                        which(cohen.shrunk==T))])]



#... Cohen Criteria
cohen.unshrunk <- c(abs(unshrunk.pcor) >= 0.5)
cohen.shrunk <- c(abs(opt.pcor ) >= 0.5)

counts <- table(cohen.unshrunk, cohen.shrunk)
rownames(counts) <-
  c("both not strong", "both strong")
barplot(
  counts,
  xlab = "?",
  col = c("darkblue", "red"),
  legend = rownames(counts)
)
#cohen.both<-names[ unique(sm.index(diag(p))[setdiff(which(cohen.unshrunk==T),
#                                        which(cohen.shrunk==T))])]

## top edges
test.results.UNSHRUNK <- network.test.edges(vec2sm(unshrunk.pcor),  plot = FALSE)
test.results.SHRUNK <- network.test.edges(vec2sm(opt.pcor),  plot = FALSE)
net.UNSHRUNK <- extract.network(test.results.UNSHRUNK, method.ggm="number", cutoff.ggm=15)
net.SHRUNK <- extract.network(test.results.SHRUNK, method.ggm="number", cutoff.ggm=15)

data.frame(net.UNSHRUNK[,c(2, 3)], net.SHRUNK[,c(2, 3)])
# geom_hline(
#   yintercept = -log10(alpha),
#   size = 1 ,
#   linetype = "dashed",
#   color = "grey") +
# geom_hline(
#   yintercept = -log10(alpha),
#   size = 1 ,
#   linetype = "dashed",
#   color = "grey") +
# geom_vline(
#   xintercept = -log10(alpha),
#   size = 1 ,
#   linetype = "dashed",
#   color = "grey")


#......................
# hist(p.UNSHRUNK)
# hist(p.opt, add = T, col = rgb(1, 0, 0, 0.5))
# 
# hist(p.adjust(p.UNSHRUNK, 'BH', length(p.UNSHRUNK)))
# hist(p.adjust(p.opt,  'BH', length(p.opt)), add = T, col = rgb(1, 0, 0, 0.5))
# 
# counts <- table(p.UNSHRUNK >= alpha, p.opt >=alpha)
# rownames(counts) <-
#   c("both not strong", "both strong")
# barplot(
#   counts,
#   xlab = "?",
#   col = c("darkblue", "red"),
#   legend = rownames(counts)
# )
# 
# #.....................................
# #shared edges
# which(p.UNSHRUNK <= alpha & p.opt <= alpha)
# 
# # only edges in unshrunk
# which(p.UNSHRUNK <= alpha & p.opt > alpha)
# edges.only.unshrunk<-paste0(names[sm.index(diag(p))[which(p.UNSHRUNK <= alpha &
#                                        p.opt > alpha), 1]], "-",
#        names[sm.index(diag(p))[which(p.UNSHRUNK <= alpha &
#                                        p.opt > alpha), 2]])
# #write.table(x = data.frame(edges.only.unshrunk, 
#                      p.UNSHRUNK[which(p.UNSHRUNK <= alpha & p.opt > alpha)], 
#                      p.opt[which(p.UNSHRUNK <= alpha & p.opt > alpha)] ) ,
#       file = paste0(0, "only_edges_ecoli_Unshrunk.txt"),sep ="\t",  row.names = FALSE, quote = FALSE)
# #...............................................................
# write(x = names[unique(sm.index(diag(p))[which(p.UNSHRUNK <=
#                                                  alpha & p.opt > alpha)])] ,
#      file = paste0(0, "more connected_ecoli_Unshrunk.txt"))

# only edges in shrunk
#which(p.UNSHRUNK > alpha & p.opt <= alpha)
#...............................................................
#venn diagram

# vennDiagram(
#   data.frame(p.UNSHRUNK <=alpha, p.opt <= alpha),
#   include = "both",
#   names = c("unshrunk", "shrunk"),
#   cex = 1,
#   counts.col = "black")


#write(x = sig.unshrunk ,
#     file = paste0(0, "ecoli_Unshrunk.txt"))
#write(x = sig.shrunk ,
#       file = paste0(opt.lambda, "ecoli_Shrunk.txt"))

#write.table(x = enrichmentGO.unshrunk,
#            file = paste0(0, "GOecoli_Unshrunk.txt"),sep ="\t",  row.names = FALSE, quote = FALSE)
# write.table(x = enrichmentGO.shrunk,
#      file = paste0(0, "GOecoli_shrunk.txt"),sep ="\t",  row.names = FALSE, quote = FALSE)
#write.table(x = enrichmentKEGG.unshrunk,
#       file = paste0(0, "KEGGecoli_Unshrunk.txt"),sep ="\t",  row.names = FALSE, quote = FALSE)
#write.table(x = enrichmentKEGG.shrunk,
#       file = paste0(0, "KEGGecoli_shrunk.txt"),sep ="\t", row.names = FALSE, quote = FALSE)

#..... PCA on cor

# The eigenvectors represent the principal components of S. 
# The eigenvalues of S are used to find the proportion of the total variance explained by the components.

cumsum(D)/sum(D)
plot(D, pch=20, ylab="eigenvalues")
df<-data.frame(V[,1],V[,2])
colnames(df)<-c("PC1","PC2")
plot(PC2~PC1, data = df, cex=0.25)
with(df, text(PC2~PC1, labels =names),pos = 2, size=0.5)
cbind(V[,1],V[,2])>0.5

# on pcor
GGM1<-vec2sm(unshrunk.pcor)
diag(GGM1)<-1
V2<-eigen(GGM1, symmetric=TRUE)$vectors
D2<-eigen(GGM1, symmetric=TRUE)$values
cumsum(D2)/sum(D2)
plot(D2, pch=20, ylab="eigenvalues")
#negative!
df<-data.frame(V2[,1],V2[,2])
colnames(df)<-c("PC1","PC2")
plot(PC2~PC1, data = df, cex=0.25)
with(df, text(PC2~PC1, labels =names),pos = 2, size=0.5)

cbind(V2[,1],V2[,2])>0.5

# high corr, weak pcorr
hist(sm2vec(cor(tcell.10)))
hist(opt.pcor/(1-opt.lambda),add=T, col=rgb(1, 0, 0, 0.5))
hist(unshrunk.pcor, add=T, col=rgb(1, 1, 0, 0.5))

# coor vs pcor
plot(sm2vec(pcor.shrink(tcell.10))/(1-opt.lambda),sm2vec(cor(tcell.10)),
     xlim=c(-(1),(1)),ylim=c(-1, 1), pch=20,
     ylab=c("corr"), xlab=c("shrink pcorr"))
points(unshrunk.pcor, sm2vec(cor(tcell.10)),col=rgb(1, 0, 0, 0.2))
abline(v=0.5*c(-1, 1), col="grey",lty = "dashed")
abline(h=0.5*c(-1, 1), col="grey", lty = "dashed")



counts <-table(idx.unshrunk[1:top]%in% idx.shrunk[1:top])
rownames(counts)<- c("both not top", "both top")
barplot(counts, xlab="?", col=c("darkblue","red"),
        legend = rownames(counts))

#top.both<-names[ unique(sm.index(diag(p))[setdiff(which(cohen.unshrunk==T),
#                                        which(cohen.shrunk==T))])]



#... Cohen Criteria
# cohen.unshrunk<- c(abs(unshrunk.pcor)>=0.3)
# cohen.shrunk<-c(abs(opt.pcor)>=0.3)
# counts <-table(cohen.unshrunk, cohen.shrunk)
# rownames(counts)<- c("both not strong", "both strong")
# barplot(counts, xlab="?", col=c("darkblue","red"),
#         legend = rownames(counts))
#cohen.both<-names[ unique(sm.index(diag(p))[setdiff(which(cohen.unshrunk==T),
#                                        which(cohen.shrunk==T))])]




# percentage of sig edges and number of genes
# sig.unshrunk<-names[ unique(sm.index(diag(p))[which(p.UNSHRUNK<=alpha)]) ]
# sig.shrunk<-names[ unique(sm.index(diag(p))[which(p.opt<=alpha)] )]
# 
# setdiff(sig.unshrunk, sig.shrunk)          
# Significant at 0.01 are small shrunk effects
r.pval.shrunk = 0.0795
p.montecarlo(r = 0.0795, number = 1000, p=102, n = 9, lambda = opt.lambda)
p.shrunk( r = 0.0795, p=102, n = 9, lambda = opt.lambda)

r.pval.unshrunk = 0.0815
p.montecarlo(r = 0.0815, number = 1000, p=102, n = 9, lambda = 4.846791e-14)
p.shrunk( r = 0.0815, p=102, n = 9, lambda = 4.846791e-14)

# Parallel coordinates: the order changes
idx.shrunk = order( abs(opt.pcor) , decreasing = TRUE)
idx.unshrunk = order( abs(unshrunk.pcor), decreasing = TRUE)

# L2 distance
cor(opt.pcor, unshrunk.pcor)
sum(opt.pcor-unshrunk.pcor)^2
cor(idx.shrunk , idx.unshrunk, method = "kendall")
head(data.frame('shrunk'= idx.shrunk , 'unshrunk'=idx.unshrunk) , n = 10)                        


# Greatest changes
changes = order( abs(opt.pcor
                     - unshrunk.pcor ) , 
                 decreasing = TRUE)[1:20]

ecoli.parcoord = data.frame('Method'= c(rep('Shrunk', length(opt.pcor[changes] )), 
                                        rep('Unshrunk', length(unshrunk.pcor[changes]))) , 
                            'pcor' = c(opt.pcor[changes] ,unshrunk.pcor[changes]) ,
                            'pair'= c( rep(1:length(opt.pcor[changes]) , 2) ))
ecoli.parcoord$Method = factor(ecoli.parcoord$Method, levels=unique(ecoli.parcoord$Method))


ggplot(data= ecoli.parcoord , aes(x= Method, y= pcor, color=Method)) + 
  geom_point(size= 2.5 ) +
  geom_line(data= ecoli.parcoord ,  aes(x= Method, y= pcor, group = pair ), lwd=1)+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.line = element_line(size = 1, colour = "grey80"),
        #axis.title =     element_blank(),
        axis.title.x = element_text(color = "grey20", size = 14),
        axis.title.y = element_text(color = "grey20", size = 14),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "black"))  +
  scale_color_manual(values= c( rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75) , rgb(0, 0, 0, 1) ))+
  scale_fill_manual(values= c( rgb(0, 0, 1, 0.5), rgb(0.8, 0.2, 0.5, 0.1) , rgb(1, 1, 1, 0.1) )) +
  scale_y_continuous(limits =c(min(ecoli.parcoord$pcor),max(ecoli.parcoord$pcor)) )+
  annotate(geom="text", x=1, y=0.85, label=paste0(expression("\u03BB")," = ", signif(opt.lambda, 2) ),
           color="black", size=5)+
  annotate(geom="text", x=1, y=0.75, label=paste0("n = ",n),
           color="black", size=5)+
  geom_hline(yintercept=0, linetype="dashed", color = "black")


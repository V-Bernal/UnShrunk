#------------------
# - Title/Project: Rebuttal to the reviewers / The "un-shrunk" partial correlation in Gaussian Graphical Models 
# - Author            : Victor Bernal*, Rainer Bischoff, Victor Guryev, Peter Horvatovich, Marco Grzegorczyk.
# - Date created      : 24 MAY 2021
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
#-----------------
# Notes: 
# - d a is a n x p matrix of data 
# - the p columns of d are the random variables 
# - the n rows of d are samples
#-----------------
# References
# [1] Sch?fer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
#------------------

#----------------------------------------
# Reviewer 1
#----------------------------------------

#----------------------------------------
# Figure 1a.
# The order of the partial correlations changes with lambda: Toy example
#----------------------------------------


rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "huge", 
              "stats4", 
              "devEMF")

lapply(Packages, library, character.only = TRUE)

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

# The toy matrix has eigenvalues, 2-norm condition number,and the real partial correlation
m
solve(m)*97

eigen(m)$values
kappa(m)
pcor0<-sm2vec(cor2pcor(m))

# Reconstruct the shrunk partial correlations
lambda<-c(1:40)/40

pcors <- 
  lapply( lambda , function( z ){
                          return( sm2vec( cor2pcor( (1-z) * m + z * diag( ncol(m) ))) )
    }
  )

pcors<-
    matrix(unlist(pcors),0.5*ncol(m)*(ncol(m)-1),length(lambda))


#emf( file = paste0( "Rev_Figure1_a_LW.emf" ), width = 5 , height = 5 )

plot.ts( t(pcors),
          plot.type="single",
          type="l",
          lw=4,
          lty=1,
          col=rainbow( 0.5 * ncol(m) * (ncol(m)-1) ),#"grey40",#rainbow(0.5*ncol(m)*(ncol(m)-1)),
          ylab="partialcorrelation",
          xlab="LW-shrinkage",
          ylim=c(-1,1),
          axes=F,
          cex = 1 )
  points( rep(0.35,length(pcor0)),
            pcor0,
            pch=20,
            cex=2,#,
            bg="grey",
            col=rainbow(0.5*ncol(m)*(ncol(m)-1))#"grey40"# 
          )
  axis(2)
  axis(1, 
       labels = c(lambda),
       at = seq( from=1, by=1,length.out = length(lambda))
      )
box()

#dev.off()

#----------------------------------------
# Compare with Glasso
#----------------------------------------
pcors <- huge(m,method="glasso",lambda=sort(lambda))
pcors <- sapply(pcors$icov,function(x){-sm2vec(cov2cor(x))})
pcors<-
    matrix(unlist(pcors),0.5*ncol(m)*(ncol(m)-1),length(lambda))

#emf(file=paste0("Rev_Figure1_a_LASSO.emf"),width=5,height=5)

plot.ts( t(pcors),
          plot.type="single",
          type="l",
          lw=4,
          lty=1,#,1:2,
          col=rainbow(0.5*ncol(m)*(ncol(m)-1)),#"grey40",
          ylab="partialcorrelation",
          xlab="glassoshrinkage",
          ylim=c(-1,1),
          axes=F,
          cex=1
        )
points( rep(0.35,length(pcor0)),
        pcor0,
        pch=20,
        cex=2,#,
        bg="grey",
        col=rainbow(0.5*ncol(m)*(ncol(m)-1))
        )
  axis(2)
  axis(1,labels=c(lambda), 
       at=seq( from=1, by=1, length.out=length(lambda) ))
box()

#dev.off()


#----------------------------------------
# Correlations
#----------------------------------------

cor0 <- sm2vec(cov2cor(m))

cors <- lapply(lambda,function(z){
  
  return( sm2vec(cov2cor(( 1 - z ) * m + z * diag( ncol(m) )) ) )
  
  }
  
  )

cors<-
  matrix( unlist(cors),0.5 * ncol(m) * (ncol(m)-1), length(lambda) )


#emf(file=paste0("Rev_Figure1_a_LW.emf"),width=5,height=5)

plot.ts( t(cors),
          plot.type="single",
          type="l",
          lw=4,
          lty=1,
          col=rainbow(0.5*ncol(m)*(ncol(m)-1)),#"grey40",#rainbow(0.5*ncol(m)*(ncol(m)-1)),
          ylab="correlation",
          xlab="LW-shrinkage",
          ylim=c(-1,1),
          axes=F,
          cex=1
        )
points( rep(0.35,length(cor0)),
        cor0,
        pch=20,
        cex=2,
        bg="grey",
        col=rainbow(0.5*ncol(m)*(ncol(m)-1))#"grey40"#
        )
  axis(2)
  axis(1,
  labels=c(lambda),
  at=seq(
  from=1,
  by=1,
  length.out=length(lambda)
  ))
box()

#dev.off()

#----------------------------------------
# Compare with Glasso
#----------------------------------------

cors <- huge(x = m, method = "glasso", lambda = sort(lambda) )
cors <- sapply(X = cors$icov, FUN = function(x){ -sm2vec(cor2pcor(x)) })
cors<-
  matrix(data = unlist(cors), nrow = 0.5*ncol(m)*(ncol(m)-1), ncol = length(lambda))


#emf(file=paste0("Rev_Figure1_a_LASSO.emf"),width=5,height=5)

plot.ts(
  t(cors),
  plot.type="single",
  type="l",
  lw=4,
  lty=1,#,1:2,
  col=rainbow(0.5*ncol(m)*(ncol(m)-1)),#"grey40",
  ylab="correlation",
  xlab="glassoshrinkage",
  ylim=c(-1,1),
  axes=F,
  cex=1
)
points(
  rep(0.35,length(cor0)),
  cor0,
  pch=20,
  cex=2,#,
  bg="grey",
  col=rainbow(0.5*ncol(m)*(ncol(m)-1))
)
  axis(2)
  axis(1,
  labels=c(lambda),
  at=seq(
  from=1,
  by=1,
  length.out=length(lambda)
  ))
box()

#dev.off()

#----------------------------------------
# Figure S5: ROC and PR curves
#----------------------------------------

# ROC curves show the trade-off between the TPR and FPR
# Precision Recall curves show the trade-off between the TPR (recall) and PPV (precision)
# ROC curves are appropriate when the observations are balanced between each class,
# PR curves are appropriate for imbalanced datasets.

rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "huge", 
              "stats4", 
              "devEMF")

lapply(Packages, library, character.only = TRUE)

source(file = "functions_Unshrink.R")

set.seed(1)

p<-100
n1<-10
etaA<-0.05

TrueNET<-ggm.simulate.pcor(p,etaA)
positives.idx<-which(sm2vec(TrueNET)!=0)
non.positives.idx<-which(sm2vec(TrueNET)==0)

# the simulated partial correlations
plot(sm2vec(TrueNET)[positives.idx],pch=16,ylab='partial correlations')
text(x=3,y=0,labels=paste0('pcors=',length(positives.idx)), cex = 1)
cat('positives=',sm2vec(TrueNET)[positives.idx])

#-----------------------
# Simulate data
sim.data1<-
  lapply(n1,function(x){
  ggm.simulate.data(x,TrueNET)
})

# Shrunk partial correlations
GGM1.shrunk<-
  lapply(sim.data1,function(x){
  sm2vec(pcor.shrink(x,verbose=FALSE))
})

# optimal shrinkage value
shrunk.lambdas<-
  lapply(sim.data1,function(x){
  return(attr(pcor.shrink(x,verbose=FALSE),"lambda"))
})

# Un-shrunk partial correlation
unshrunk1<-
  lapply(sim.data1,function(x){
  unshrink_GGM(data = x, PLOT = F , l_0 = 0.01, corrected = T)
})
shrunk.lambdas[[1]]

#---------------------------
# Glasso partial correlations
glasso<-
  lapply(sim.data1,function(x){
  -sm2vec(cov2cor(
  huge.select(huge(x,method="glasso"),criterion="stars")$opt.icov
  ))
})


#--------------------------
# Compute ROC
simp_roc1<-simple_roc(c(sm2vec(TrueNET)!=0),abs(GGM1.shrunk[[1]]))
Usimp_roc1<-simple_roc(c(sm2vec(TrueNET)!=0),abs(unshrunk1[[1]]))
glasso_roc1<-simple_roc(c(sm2vec(TrueNET)!=0),abs(glasso[[1]]))

simp_pr1<-simple_PR(c(sm2vec(TrueNET)!=0),abs(GGM1.shrunk[[1]]))
Usimp_pr1<-simple_PR(c(sm2vec(TrueNET)!=0),abs(unshrunk1[[1]]))
glasso_pr1<-simple_PR(c(sm2vec(TrueNET)!=0),abs(glasso[[1]]))

#--------------------------

#emf(file=paste0("ROC_p",p,"n_",n1,"etA_",etaA,".emf"),width=7,height=7)
plot(simp_roc1$FPR,simp_roc1$TPR, pch=20,
  col="black",type="b",cex=0.35,
  ylim=c(0,1),xlim=c(0,1),
  ylab="TPR",xlab="FPR")

  #abline(a=0,b=1)
   lines(glasso_roc1$FPR,glasso_roc1$TPR,col="red",pch=20,type="b",cex=0.35)
   lines(Usimp_roc1$FPR,Usimp_roc1$TPR,col="blue",pch=20,type="b",cex=0.35)
   
  text(x=c(0.75),y=0.25,cex=01,
  labels=paste0('shrunk AUC=',signif(simple_auc(simp_roc1$TPR,simp_roc1$FPR), digits = 4)))
  text(x=c(0.75),y=0.15,cex=1,
  labels=paste0('un-shrunk AUC=',signif(simple_auc(Usimp_roc1$TPR,Usimp_roc1$FPR), digits = 4)))
  text(x=c(0.75),y=0.05,cex=1,
  labels=paste0('glasso AUC=',signif(simple_auc(glasso_roc1$TPR,glasso_roc1$FPR), digits = 4)))
  # text(x=c(0.25),y=0.05,cex=0.75,
  # labels=paste0('p=',p,"n=",n1,"etA=",etaA))
  legend(x =0.25, y =  0.5,
  legend=c("shrunk","un-shrunk","glasso"),
    col=c(rgb(0,0,0,0.7),
    rgb(0,0,1,0.7),
    rgb(1,0,0,0.7)),
    pch=c(20,20,20),
    #bty="n",
    #pt.cex=1,
    cex=1,
    text.col="black",
    horiz=T,
    #inset=c(0.1,0.1),
    y.intersp = 0.2)

dev.off()
  
#-------------------------
# PR curves
#emf(file=paste0("PR_p",p,"n_",n1,"etA_",etaA,".emf"),width=4,height=4)
plot(simp_pr1$PPV,simp_pr1$TPR,col="black",type="b",
  ylim=c(0,1),xlim=c(0,1),cex=0.5,
  ylab="PPV",xlab="TPR")
  
  abline(a=0,b=1)
  lines(Usimp_pr1$PPV,Usimp_pr1$TPR,col="blue",pch=20,type="b",cex=0.5)
  lines(glasso_pr1$PPV,glasso_pr1$TPR,col="red",pch=20,type="b",cex=0.5)
  text(x=c(0.20),y=0.75,cex=0.75,
  labels=paste0('shrunk AUC=', signif(simple_auc(simp_pr1$PPV,simp_pr1$TPR), digits = 4)))
  text(x=c(0.2),y=0.5,cex=0.75,
  labels=paste0('un-shrunk AUC=',signif(simple_auc(Usimp_pr1$PPV,Usimp_pr1$TPR), digits = 4)))
  text(x=c(0.2),y=0.25,cex=0.75,
  labels=paste0('glasso AUC=',signif(simple_auc(glasso_pr1$PPV,glasso_pr1$TPR), digits = 4)))
  text(x=c(0.2),y=0.05,cex=0.75,
  labels=paste0('p=',p,"n=",n1,"etA=",etaA))
  legend('center',
  legend=c("shrunk","un-shrunk","glasso"),
    col=c(rgb(0,0,0,0.7),
    rgb(0,0,1,0.7),
    rgb(1,0,0,0.7)),
    pch=c(20,20,20),
    bty="n",
    pt.cex=0.75,
    cex=0.7,
    text.col="black",
    horiz=F,
    inset=c(0.1,0.1),
    y.intersp=0.2)

#--------------------------
# Table 1.- ROC and PR tables
#------------------------

# ROC curves show the trade-off between the TPR and FPR
# Precision Recall curves show the trade-off between the TPR (recall) and PPV (precision)
# ROC curves are appropriate when the observations are balanced between each class,
# PR curves are appropriate for imbalanced datasets.
  
rm(list=ls())
  
Packages <- c("GeneNet", 
              "ggplot2", 
              "huge", 
              "stats4", 
              "devEMF")
  
lapply(Packages, library, character.only = TRUE)

source(file = "functions_Unshrink.R")
  
set.seed(1)

p <- 100
n1 <- seq( from = 10, to = 100, by = 10 )
etaA <- 0.03

# Fixed Simulate true network        
TrueNET <- ggm.simulate.pcor( p , etaA )
positives.idx <- which( sm2vec(TrueNET) != 0 )
non.positives.idx <- which( sm2vec(TrueNET) == 0 )

list_result = lapply( X = n1 , FUN = function(X){ roc_fun(X) } )
new_list_roc = unlist( list_result )

# list_result = lapply( X = n1 , FUN = function(X){ truncated_roc_fun(X) } )
# new_list_roc = unlist( list_result )

df = data.frame(  'sample'= n1,
                  'Unshrunk'= new_list_roc[names(new_list_roc)=='Unshrunk'],
                  'shrunk'= new_list_roc[names(new_list_roc)=='shrunk'],
                  'glasso'= new_list_roc[names(new_list_roc)=='glasso'],
                  'U-S' = new_list_roc[names(new_list_roc)=='U-S'],
                  'U-g' = new_list_roc[names(new_list_roc)=='U-g'])

df

write.table(df$shrunk , file = 'roc.txt',quote = F, row.names = F,col.names = F)

# f1 score
# set.seed(1)
# list_resultf1 = lapply( X = n1 , FUN = function(X){ f1_fun(X) } )
# new_list_f1 = unlist( list_resultf1 )
# 
# df1 = data.frame(  'sample'= n1,
#                   'Unshrunk'= new_list_f1 [names(new_list_f1 )=='Unshrunk'],
#                   'shrunk'= new_list_f1 [names(new_list_f1 )=='shrunk'],
#                   'glasso'= new_list_f1 [names(new_list_f1 )=='glasso'],
#                   'U-S' = new_list_f1 [names(new_list_f1 )=='U-S'],
#                   'U-g' = new_list_f1 [names(new_list_f1 )=='U-g'])
# 
# df1
# 
# write.table(df1$glasso , file = 'f1.txt',quote = F, row.names = F,col.names = F)

ggplot(df) + 
  geom_point(aes(x = sample, y = Unshrunk - shrunk), color = "black", size = 1)  +
  geom_point(aes(x = sample, y = Unshrunk - glasso ), color = "grey", size = 1) +
  theme_classic(base_size = 10, base_family = "" ) +
  labs(x = "sample size", y = "ROC AUC diff") + scale_x_continuous( breaks = n1 )


#----------------------------------------------
# Figure 4. Bland Altman plots
#----------------------------------------------
rm(list=ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "huge", 
              "stats4", 
              "devEMF")

lapply(Packages, library, character.only = TRUE)

p<-100
n1<-c(40, 2*40)
etaA<-0.02

TrueNET<-ggm.simulate.pcor(p,etaA)#ggm.simulate.pcor.MODIF
TrueNET
positives.idx<-which(sm2vec(TrueNET)!=0)

source(file = "functions_Unshrink.R")
source(file = "functions_pvalues.R")

set.seed(1)

positives.idx<-which(sm2vec(TrueNET)==0)

# Simulate data
sim.data1<-
  lapply(n1,function(x){
    ggm.simulate.data(x,TrueNET)
  })

# Concatenate the data with itself. The new data has double the sample size.
sim.data1[[2]] = rbind(sim.data1[[1]], sim.data1[[1]] +  0.0005 * sim.data1[[2]])
sim.data1[[3]] = rbind(sim.data1[[1]], sim.data1[[1]])

# Shrunk partial correlation
GGM1.shrunk<-
  lapply(sim.data1,function(x){
    sm2vec(pcor.shrink(x,verbose=FALSE))
  })

# Optimal shrinkage value
shrunk.lambdas<-
  lapply(sim.data1,function(x){
    return(attr(pcor.shrink(x,verbose=FALSE),"lambda"))
  })


# Un-shrunk partial correlation
unshrunk1<-
  lapply(sim.data1,function(x){
    unshrink_GGM(d = x,  PLOT=F,
                         l_0 = 0.01, 
                         corrected = F )
  })

# Glasso partial correlation
glasso<-
  lapply(sim.data1,function(x){
    -sm2vec(cov2cor(
      huge.select(huge(x,method="glasso"),criterion="stars")$opt.icov
    ))
  })

#------------------------------------
# Compare the partial correlation coefficients and the p-values
#------------------------------------
pval.UNSHRUNK <- p.standard(unshrunk1[[1]], p = p, n = n1[1])
pval.shrunk <- p.shrunk(GGM1.shrunk[[1]], p = p, n = n1[1], shrunk.lambdas[[1]])

pval.UNSHRUNK2 <- p.standard(unshrunk1[[2]], p, n1[2])
pval.shrunk2 <- p.shrunk(GGM1.shrunk[[2]], p, n1[2], shrunk.lambdas[[2]])

pval.UNSHRUNK3 <- p.standard(unshrunk1[[3]], p, n1[2])
pval.shrunk3 <- p.shrunk(GGM1.shrunk[[3]], p, n1[2], shrunk.lambdas[[3]])


plot(pval.UNSHRUNK2 - pval.UNSHRUNK)
points(pval.UNSHRUNK3 -pval.UNSHRUNK, col=3)

#------------------------------------
# Create data frame
id <- sm.index(matrix(0, p , p), diag = FALSE)
df <- data.frame(id[, 1], id[, 2])
names = 1:p

edge <-
  data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))

scat <- 
  data.frame(edge, 
             GGM1.shrunk[[1]], GGM1.shrunk[[2]],GGM1.shrunk[[3]],
             unshrunk1[[1]], unshrunk1[[2]], unshrunk1[[3]],
             pval.shrunk, pval.shrunk2, pval.shrunk3,
             pval.UNSHRUNK, pval.UNSHRUNK2 , pval.UNSHRUNK3)

colnames(scat) <- c("Edges", "shrunk", "shrunk2","shrunk3",
                    "unshrunk", "unshrunk2", "unshrunk3",
                    "pval.shrunk", "pval.shrunk2", "pval.shrunk3", 
                    "pval.unshrunk", "pval.unshrunk2", "pval.unshrunk3")

scat$Edges <- as.character(scat$Edges)

#------------------------------------
# Cohen criterium | pcor | >  0.1 and p values < 0.05

sum(scat$pval.shrunk<0.05 & scat$pval.shrunk2<0.05)
sum(scat$pval.unshrunk<0.05 & scat$pval.unshrunk2<0.05)

sum(abs(scat$shrunk)>0.1 & abs(scat$shrunk2)>0.1)
sum(abs(scat$unshrunk)>0.1 & abs(scat$unshrunk2)>0.1)

scat$Edges[which( abs(scat$shrunk) < 0.2 &
                    abs(scat$unshrunk) < 0.2)] <- NA

scat$Edges[which( abs(scat$pval.shrunk) > 0.05 &
                    abs(scat$pval.unshrunk) > 0.05)] <- NA

#--------------------
# Figure 4. Shrunk and Unshrunk - Bland Altman plot with different sample sizes

# With noise
fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(pval.shrunk + pval.shrunk2), y = (pval.shrunk - pval.shrunk2)), color = rgb(0.5,0.5,0.5,0.25), size = 1) +
  geom_point(aes(x = 0.5*(pval.unshrunk + pval.unshrunk2), y = (pval.unshrunk -  pval.unshrunk2)), color = rgb(1,0,0,0.25), size = 1) +
  labs(x = "ave p-value", y = "diff p-value") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black") +
  ylim(-1,1)
fig

fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(shrunk + shrunk2), y = (shrunk - shrunk2)), color = rgb(0.5,0.5,0.5,0.25), size = 1) +
  geom_point(aes(x = 0.5*(unshrunk + unshrunk2), y = (unshrunk - unshrunk2)), color = rgb(1,0,0,0.25), size = 1) +
  labs(x = "ave pcor", y = "diff pcor") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black") 
fig

# With-out noise
emf("Figure_4_DIFF_lambdas_c.emf", width = 4, height = 4)
fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(pval.shrunk + pval.shrunk3), y = (pval.shrunk - pval.shrunk3)), color = rgb(0.5,0.5,0.5,0.25), size = 3) +
  labs(x = "ave p-value", y = "diff p-value") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black")+ ylim(-1,1) 
  
fig
dev.off()

emf("Figure_4_DIFF_lambdas_d.emf", width = 4, height = 4)
fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(pval.unshrunk + pval.unshrunk3), y = (pval.unshrunk -  pval.unshrunk3)), color = rgb(0.75,0,0,0.25), size = 3) +
  labs(x = "ave p-value", y = "diff p-value") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black")+ ylim(-1,1) 
fig
dev.off()

emf("Figure_4_DIFF_lambdas_a.emf", width = 4, height = 4)
fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(shrunk + shrunk3), y = (shrunk - shrunk3)), color = rgb(0.5,0.5,0.5,0.25), size = 3) +
  labs(x = "ave pcor", y = "diff pcor") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black")+ ylim(-1,1) 
fig
dev.off()

emf("Figure_4_DIFF_lambdas_b.emf", width = 4, height = 4)
fig<-ggplot(data = scat) + 
  geom_point(aes(x = 0.5*(unshrunk + unshrunk3), y = (unshrunk - unshrunk3)), color = rgb(0.75,0,0,0.25), size = 3) +
  labs(x = "ave pcor", y = "diff pcor") +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= 0, linetype="dashed", color = "black") +ylim(-1,1)
fig
dev.off()


# How large are the differences?
max(abs(scat$unshrunk-scat$unshrunk3))
max(abs(scat$pval.unshrunk-scat$pval.unshrunk3))

#---------------------------------------
# Reviewer 2
#---------------------------------------

#---------------------------------------
# Figure 1a.
# Log transformation and partial correlations
#---------------------------------------

# Simulate expression
x=rpois(n=100,lambda=10)+5
y=x+0.6*rpois(n=100,lambda=10)+1

# Correlation
cor(x,y)
cor(log2(x),log2(y))

# log transformation keeps the trend,but change the correlation
par(mfrow=c(1,2),pty="s")

plot(x,y,col='blue',pch=20,bg='black')
text(x=20,y=15,labels=round(cor(x,y),4))

plot(log2(x),log2(y),col='orange',pch=20,bg='black')
text(x=log2(20),y=log2(15),labels=round(cor(log2(x),log2(y)),4))


# Influencial points can behave differents impact in the log scale
y[which.max(x)]=100

par(mfrow=c(1,2),pty="s")

plot(x,y,col='blue',pch=20,bg='black')
text(x=15,y=y[which.max(x)],labels='influencepoint->')
text(x=20,y=15,labels=round(cor(x,y),4))

plot(log2(x),log2(y),col='orange',pch=20,bg='black')
text(x=log2(15),y=log2(y[which.max(x)]),labels='influencepoint->')
text(x=log2(20),y=log2(15),labels=round(cor(log2(x),log2(y)),4))



#........................................
# 5) Real data: E. coli example
#........................................

rm(list = ls())

Packages <- c("GeneNet", 
              "ggplot2", 
              "igraph", 
              "stats4", 
              "limma", 
              "STRINGdb",
              "devEMF")

lapply(Packages, library, character.only = TRUE)

source("functions_Unshrink.R")
source("functions_pvalues.R")


#Load data
data(ecoli)

#
p <- ncol(ecoli)
n <- nrow(ecoli)
names <- colnames(ecoli)
ecoli <- scale(ecoli)

# Estimate GGM (i.e. the partial correlation coefficients)
shrunk.pcor <- sm2vec(pcor.shrink(ecoli))
lambda <- attr(pcor.shrink(ecoli), "lambda")
unshrunk.pcor <- unshrink_GGM(d = ecoli, PLOT = F, l_0 = 0.01, corrected = F)

# Compare the partial correlation coefficients, and the p-values
p.UNSHRUNK <- p.standard(unshrunk.pcor, ncol(ecoli), nrow(ecoli))
p.shrunk <- p.shrunk(shrunk.pcor, ncol(ecoli), nrow(ecoli), lambda)

# Cerate data frame
id <- sm.index(matrix(0, p , p), diag = FALSE)
df <- data.frame(id[, 1], id[, 2])

edge <-
  data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))

scat <- 
  data.frame(edge, 
             shrunk.pcor, unshrunk.pcor,
             p.shrunk, p.UNSHRUNK)

colnames(scat) <- c("Edges", "shrunk", "unshrunk", "pval.shrunk", "pval.unshrunk")

scat$Edges <- as.character(scat$Edges)

# Volcano plot
# filter using Cohen 0.1 and p values 0.05
scat$Edges[which( abs(scat$unshrunk) < 0.4)] <- NA

scat$Edges[which( abs(scat$pval.unshrunk) > 0.05)] <- NA

# plot the partial correlation coefficients. 
fig<-ggplot( scat, aes(x = unshrunk, y = pval.unshrunk)) + 
  geom_point(aes(x = unshrunk, y = -log10(pval.unshrunk)), size = 2) +
  geom_point(aes(x = shrunk, y = -log10(pval.shrunk)), color="grey", size = 2) +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "black")+
  geom_vline(xintercept= c(-0.1,0.1), linetype="dashed", color = "black")+
  # geom_text(aes(x = unshrunk, y = -log10(pval.unshrunk), label = Edges), 
  #            size = 2.5, nudge_x = -.05, nudge_y = 0)+
  # labs(x = "pcor", y = "-log_10(p-values)") + 
  scale_x_continuous(name="pcor", limits=c(-.6, .6), breaks = seq(-1,1,0.1)) +
  scale_y_continuous(name="-log_10(p-values)")

# Uncomment to save as .emf
emf(file = paste0("Ecoli_scatter.emf"), width = 12, height = 5 )
fig
dev.off()

fig+scale_x_log10()+ scale_y_log10()

#-------------------------
# How many edges are (i) statistically significant  
sum(p.shrunk < 0.05)
sum(p.UNSHRUNK < 0.05)

# How many edges are (i) statistically significant and (ii) strong effects 
sum( (abs(shrunk.pcor) > .1) & (p.shrunk < 0.05))
sum( (abs(unshrunk.pcor) > .1) & (p.UNSHRUNK < 0.05))

# better at small sample size
# BH
GGM <- vec2sm(  p.shrunk < .05 ) +
  2 * vec2sm( (abs(unshrunk.pcor) > .1) & p.UNSHRUNK < .05)
diag(GGM) <- 0

new_names = names
unconnected = which(rowMeans(GGM)==0)

if(length(unconnected)>0){
  GGM = GGM[-unconnected, -unconnected]
  rowMeans(GGM)
  new_names = names[-unconnected]
}

id.conn.s = which(GGM ==1 | GGM ==3 , arr.ind = T)
id.conn.us = which(GGM ==2 | GGM ==3 , arr.ind = T)
write.table(x = names[unique(c(id.conn.s[,1],id.conn.s[,2]))], file = 'ecoli_conn_s.txt',row.names = F,col.names = F, quote = F)
write.table(x = names[unique(c(id.conn.us[,1],id.conn.us[,2]))], file = 'ecoli_conn_us.txt',row.names = F,col.names = F, quote = F)
write.table(x = names, file = 'ecoli_all.txt',row.names = F,col.names = F, quote = F)

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

emf(file = paste0("FigureS5_ecoli_net.emf"), width = 11, height = 11 )
plot(g1,
  layout = layout_nicely(g1),
  edge.width = line.thickness[E(g1)$weight], 
  edge.color = colEdge[E(g1)$weight],
  vertex.label.color = "black",
  vertex.label.dist = -0.55,
  vertex.label.cex = 1 ,
  vertex.color = c("black"),
  vertex.size = 2 ,
  edge.label.family = "Times"
)
dev.off()


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

lapply(Packages, library, character.only = TRUE)

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

#--------------------------------------
## Network analysis
pdata_bot$strain
data.norm = t(data.norm)

# Network
shrunk.pcor <- sm2vec(pcor.shrink( data.norm , verbose = FALSE ))
lambda <-attr(pcor.shrink( data.norm, verbose = FALSE ),"lambda")
unshrunk.pcor <- unshrink_GGM( d = data.norm, PLOT = FALSE, l_0 = 0.01,corrected = F )
     
# P-values
pval.shrunk <- p.shrunk(shrunk.pcor, ncol(data.norm), nrow(data.norm), lambda)
pval.UNSHRUNK <- p.standard(unshrunk.pcor, ncol(data.norm), nrow(data.norm))
              
# Create data frame
id <- sm.index(matrix(0, p , p), diag = FALSE)
df <- data.frame(id[, 1], id[, 2])

edge <-
  data.frame(paste (names[df[, 1]], names[df[, 2]], sep = "-", collapse = NULL))

scat <- 
  data.frame(edge, 
             shrunk.pcor, unshrunk.pcor,
             pval.shrunk, pval.UNSHRUNK)

colnames(scat) <- c("Edges", "shrunk", "unshrunk", "pval.shrunk", "pval.unshrunk")

scat$Edges <- as.character(scat$Edges)

# Volcano plot
# filter using Cohen 0.1 and p values 0.05
scat$Edges[which( abs(scat$shrunk) < 0.3 &
                    abs(scat$unshrunk) < 0.3)] <- NA

scat$Edges[which( abs(scat$pval.shrunk) > 0.05 &
                    abs(scat$pval.unshrunk) > 0.05)] <- NA

# plot the partial correlation coefficients. 
fig<-ggplot(scat, aes(x = unshrunk, y = pval.unshrunk)) + 
  geom_point(aes(x = unshrunk, y = -log10(pval.unshrunk)), size = 1) +
  geom_point(aes(x = shrunk, y = -log10(pval.shrunk)), color="grey", size = 2) +
  theme_classic(base_size = 14, base_family = "") +
  geom_hline(yintercept= -log10(0.01), linetype="dashed", color = "black")+
  geom_vline(xintercept= c(-0.1,0.1), linetype="dashed", color = "black")+ 
  scale_x_continuous(name="pcor", limits=c(-.6, .6), breaks = seq(-1,1,0.1)) +
  scale_y_continuous(name="-log_10(p-values)")

# Uncomment to save as .emf
emf(file = paste0("MMUS_scatter.emf"), width = 12, height = 5 )
fig
dev.off()

fig+scale_x_log10()+ scale_y_log10()

#-------------------------
# How many edges are (i) statistically significant  
sum(pval.shrunk < 0.01)
sum(pval.UNSHRUNK < 0.01)

sum( p.adjust(pval.shrunk, method = 'BH') < 0.05)
sum(p.adjust(pval.UNSHRUNK, method = 'BH') < 0.05)

# How many edges are (i) statistically significant and (ii) strong effects 
sum( pval.shrunk < 0.01)
sum( (abs(unshrunk.pcor) > .1) & pval.UNSHRUNK < 0.01)

# Ensembl are mapped to external gene names
# library(biomaRt)
# ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# mouse_gene_ids = names
# all_new_gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
#                       filters = 'ensembl_gene_id', values = mouse_gene_ids, mart = ensembl)
# sum(is.na(all_new_gene))
# idd=match(names , all_new_gene$ensembl_gene_id)
# temp.names = all_new_gene$external_gene_name[]
# temp.names[is.na(temp.names)] = all_new_gene$ensembl_gene_id[is.na(temp.names)]
# new_names = temp.names



# better at small sample size
GGM <- vec2sm( pval.shrunk < 0.01 ) +
  2 * vec2sm( (abs(unshrunk.pcor) > .1) & pval.UNSHRUNK < 0.01)
diag(GGM) <- 0

unconnected = which(rowMeans(GGM)==0)

if(length(unconnected)>0){
  GGM = GGM[-unconnected, -unconnected]
  rowMeans(GGM)
  new_names = names[-unconnected]
}

# connected shrunk
id.cs = which( (GGM == 1) | (GGM == 3), arr.ind = T)
id.cu = which((GGM == 2) | (GGM == 3), arr.ind = T)

write.table(x = new_names[unique(c(id.cs[,1],id.cs[,2]))], file = 'mmus_conn_s.txt',row.names = F,col.names = F, quote = F)
write.table(x = new_names[unique(c(id.cu[,1],id.cu[,2]))], file = 'mmus_conn_us.txt',row.names = F,col.names = F, quote = F)
write.table(x = new_names, file = 'mmus_conn_all.txt',row.names = F,col.names = F, quote = F)


g1 <-
  graph_from_adjacency_matrix(
    abs(GGM),
    mode = c("undirected"),
    diag = FALSE,
    weighted = TRUE
  )

V(g1)$label <- as.matrix(names)
colEdge <- c("blue", "red", "magenta")
line.thickness <- c(1, 1, 3)

emf(file = paste0("FigureS6_mmus_net.emf"), width = 15, height = 15 )
plot(
  g1,
  layout = layout_nicely(g1),#layout.sphere(g1),
  edge.width = line.thickness[E(g1)$weight], 
  edge.color = colEdge[E(g1)$weight],
  vertex.label.color = "black",
  vertex.label.dist = -0.55,
  vertex.label.cex = 0.5 ,
  vertex.color = c("black"),
  vertex.size = 2 ,
  edge.label.family = "Times"
)
dev.off()


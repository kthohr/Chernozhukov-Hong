#
# Example 1 in Chernozhukov-Hong (2003)
# Censored Median Regression
# 
# Keith O'Hara
# 11/28/14
#
rm(list=ls())
#
library(ggplot2); library(grid)
# Save Graphs:
saveG <- FALSE
#
set.seed(8981)
#
checkfn <- function(pars,y,z){
  beta <- matrix(pars)
  #
  check <- y - pmax(0, z%*%beta)
  ret <- sum(abs(check))
  #
  return(ret)
}
# Objective Function:
ObjFn <- function(pars,y,z){
  #
  Crit <- -checkfn(pars,y,z)
  #
  return(Crit)
}
# Log Posterior Kernel
PKernel <- function(pars,y,z){
  prior <- sum(dunif(pars,-10,10,log=TRUE))
  #
  laplace <- ObjFn(pars,y,z)
  # log kernel:
  ret <- laplace + prior
  #
  return(-ret)
}
# True values of the parameters
beta0 <- matrix(c(-6,3,3,3))
#
n <- 200
X <- matrix(rnorm(3*n),ncol=3)
Z <- cbind(1,X)
epsilon <- rnorm(n)
#
u <- (X[,1]^2)*epsilon
#
Ys <- Z%*%beta0 + matrix(u)
Y <- pmax(0,Ys)
#
# Initialize MCMC run
#
initialvals <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
ehat <- Y - Z%*%initialvals
sigsq <- mean(ehat^2)
CovM <- solve(t(Z)%*%Z)*sigsq^2
#
postmode <- optim(par=initialvals,fn=PKernel,y=Y,z=Z,method="Nelder-Mead",hessian=TRUE)
# Set a scaling parameter for the covariance matrix (to adjust the acceptance rate), 
# the number of draws to keep, and the number of draws to use as burn-in.
scalepar <- 1
keep <- 10000
burnin <- 10000
Draws <- matrix(NA,nrow=(keep+1),ncol=length(postmode$par))
#
InitialDraw <- postmode$par
PrevLP <- (-1)*postmode$value
#
PickMeInstead <- matrix(c(InitialDraw))
#
CovM <- scalepar*CovM
CovMChol <- t(chol(CovM))
#
Acceptances <- 0
# For illustrative purposes, break this up into burn-in and keep draws.
for (i in 1:burnin){
  #
  proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(postmode$par)))
  #
  PropLP <- (-1)*PKernel(c(proposal),Y,Z)
  if(is.nan(PropLP)){
    PropLP <- -1000000
  }
  #
  if(runif(1) < exp(PropLP-PrevLP)){
    PickMeInstead <- proposal
    PrevLP <- PropLP
  }
}
#
Draws[1,] <- t(PickMeInstead)
#
for (i in 1:keep){
  #
  proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(postmode$par)))
  #
  PropLP <- (-1)*PKernel(c(proposal),Y,Z)
  if(is.nan(PropLP)){
    PropLP <- -1000000
  }
  #
  if(runif(1) < exp(PropLP-PrevLP)){
    Draws[i+1,] <- t(proposal)
    Acceptances <- Acceptances + 1
    #
    PickMeInstead <- proposal
    PrevLP <- PropLP
  }else{
    Draws[i+1,] <- Draws[i,]
  }
}
#
Draws <- Draws[-1,]
accept <- Acceptances/keep
#
colMeans(Draws)
#
df1 <- data.frame(1:keep,Draws)
colnames(df1) <- c("Draw","Beta0","Beta1","Beta2","Beta3")
BinDenom <- 40
#
ParamBin <- (max(df1[,2]) - min(df1[,2]))/BinDenom
gg1.1 <- ggplot(df1,aes(Beta0)) + xlab(expression(beta[0])) + ylab("")
gg1.2 <- gg1.1 + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(''),sep="") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
#
ParamBin <- (max(df1[,3]) - min(df1[,3]))/BinDenom
gg2.1 <- ggplot(df1,aes(Beta1)) + xlab(expression(beta[1])) + ylab("")
gg2.2 <- gg2.1 + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(''),sep="") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
#
ParamBin <- (max(df1[,4]) - min(df1[,4]))/BinDenom
gg3.1 <- ggplot(df1,aes(Beta2)) + xlab(expression(beta[2])) + ylab("")
gg3.2 <- gg3.1 + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(''),sep="") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
#
ParamBin <- (max(df1[,5]) - min(df1[,5]))/BinDenom
gg4.1 <- ggplot(df1,aes(Beta3)) + xlab(expression(beta[3])) + ylab("")
gg4.2 <- gg4.1 + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(''),sep="") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
#
vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
#
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_ps(file="CH1.eps",height=9,width=13)}
pushViewport(viewport(layout=grid.layout(2,2)))
print(gg1.2,vp = vplayout(1,1))
print(gg2.2,vp = vplayout(1,2))
print(gg3.2,vp = vplayout(2,1))
print(gg4.2,vp = vplayout(2,2))
if(saveG==TRUE){dev.off()}
#
#
#END
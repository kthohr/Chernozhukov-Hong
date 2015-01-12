#
# Example 2 in Chernozhukov-Hong (2003)
# Quantile IV Regression
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
# Moment Conditions:
gfn <- function(pars,y,z){
  #
  alpha <- pars[1]
  beta <- matrix(pars[2:4])
  #
  D <- z[,2:4]
  #
  check <- as.numeric(1 - 1*(y > alpha + D%*%beta))
  temp <- (0.5 - check)*z
  #
  g <- matrix(colMeans(temp))
  return(g)
}
# Weighting Matrix:
Wfn <- function(pars,y,z){
  #
  alpha <- pars[1]
  beta <- matrix(pars[2:4])
  #
  D <- z[,2:4]
  #
  check <- as.numeric(1 - 1*(y > alpha + D%*%beta))
  W <- t((0.5 - check)*z)%*%((0.5 - check)*z)/nrow(z)
  W <- solve(W)
  #
  return(W)
}
# Sample Objective Function:
ObjFn <- function(pars,y,z){
  #
  gbar <- gfn(pars,y,z)
  W <- Wfn(pars,y,z)
  #
  Crit <- -t(gbar)%*%W%*%gbar*(nrow(z)/2)
  #
  return(Crit)
}
# Log Posterior Kernel:
PKernel <- function(pars,y,z){
  #
  prior <- sum(dunif(pars,-10,10,log=TRUE))
  #
  laplace <- ObjFn(pars,y,z)
  # log kernel:
  ret <- laplace + prior
  #
  return(-ret)
}
#
# First, set the true values of the parameters
#
alpha0 <- 0
beta0 <- matrix(rep(0,3))
# Sample size
n <- 200
# Now create a sample
D <- matrix(rnorm(3*n),ncol=3)
Z <- cbind(1,D)
epsilon <- rnorm(n)
#
sigmaD <- numeric(n)
for(i in 1:n){
  sigmaD[i] <- (1 + sum(D[i,]))/5
}
# Heterskedastic errors
u <- sigmaD*epsilon
#
Y <- alpha0 + D%*%beta0 + matrix(u)
#
# Plot the objective function over an interval of intercept values
#
alphaseq <- seq(-3,3,length.out=1000)
critseq <- alphaseq
#
for(jj in 1:length(alphaseq)){
  critseq[jj] <- ObjFn(c(alphaseq[jj],0,0,0),Y,Z)
}
#
df1 <- data.frame(alphaseq,-critseq)
colnames(df1) <- c("Alpha","Qn")
#
gg1.1 <- ggplot(df1,aes(Alpha)) + xlab(expression(alpha)) + ylab("Objective Function")
gg1.2 <- gg1.1 + geom_line(aes(y = Qn), colour="darkblue") + labs(title="",sep="") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
#
if(class(dev.list()) != "NULL"){dev.off()}
if(saveG==TRUE){cairo_ps(file="CH2.1.eps",height=9,width=13)}
print(gg1.2)
if(saveG==TRUE){dev.off()}
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
# Acceptance rate
accept <- Acceptances/keep
Draws <- Draws[-1,]
#
colMeans(Draws)
#
df1 <- data.frame(1:keep,Draws)
colnames(df1) <- c("Draw","Alpha","Beta1","Beta2","Beta3")
BinDenom <- 40
#
ParamBin <- (max(df1[,2]) - min(df1[,2]))/BinDenom
gg1.1 <- ggplot(df1,aes(Alpha)) + xlab(expression(alpha)) + ylab("")
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
grid.newpage()
if(saveG==TRUE){cairo_ps(file="CH2.eps",height=9,width=13)}
pushViewport(viewport(layout=grid.layout(2,2)))
print(gg1.2,vp = vplayout(1,1))
print(gg2.2,vp = vplayout(1,2))
print(gg3.2,vp = vplayout(2,1))
print(gg4.2,vp = vplayout(2,2))
if(saveG==TRUE){dev.off()}
#
#
#END
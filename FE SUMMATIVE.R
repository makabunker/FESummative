##############CANDIDATE NUMBER: 71052,61303,72039,69522##############



rm(list=ls())
#load data
load("C:/Users/kornm/Desktop/FE/fm408_exam_data.RData")
#Question 3
#BS price
BSprice<-function(pc, S, k, vol, d, r, t)
{
  
  d1 <- (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  d2 <- d1 - vol * sqrt(t)
  
  BSprice <- pc * exp(-d * t) * S * 
    pnorm(pc * d1) - pc * k * exp(-r * t) * pnorm(pc * d2)
  return(BSprice)
}

########################## QUESTION 3 ##########################
#The interest and dividend rates are in annualized terms, so we need to make them continuous.
raw["Cinterest"] <- log(1+raw$InterestRate)
raw["Cdiv"] <- log(1+raw$DividendYield)
#find stock price
raw["St"] <- raw$Forward/exp((raw$Cinterest-raw$Cdiv)*(raw$Term/12))
#extract observed IV by using BS and uniroot
IV <- function(x)
{ 
  IV <- BSprice(1,St,K,x,div,int,t/12)-C
}
k <- nrow(raw)
for (i in 1:k)
{
  St <- raw$St[i]
  K <- raw$Strike[i]
  div <- raw$Cdiv[i]
  int<- raw$Cinterest[i]
  t<- raw$Term[i]
  C <- raw$CallPrice[i]
  Y <- uniroot(IV, interval = c(0,2))$root
  raw[i,11] <- Y
}
colnames(raw)[11] <- "IV"
raw["moneyness"] <- log(raw$Strike/raw$Forward)

#3d Scatter plot of the implied volatility 
library(plot3Drgl)
par(mfrow=c(1,1))
scatter3Drgl(x=raw$moneyness,y=raw$Term,z=raw$IV,colvar=raw$IV, type = c("shade", "wire", "dots"), xlab="Log moneyness", ylab="Months to maturity", zlab="Implied volatility")



########################## QUESTION 4 ##########################
#define matrix to store parameters
pmatrix <- matrix(nrow=480,ncol=5)

#svi function
vfunc=function(p){
  p=pmatrix[1,]
  a=p[1]
  b=p[2]
  rho=p[3]
  m=p[4]
  sigma2=p[5]^2
  
  kk=(k-m)^2+sigma2
  totalIV=a+b*(rho * (k-m)+sqrt(kk))
  
  if(!is.numeric(totalIV)) return(1000)
  if(min(totalIV)<0) return(1000)
  if(abs(rho)>1) return(1000)
  if(sqrt(sigma2)<=0) return(1000)
  if(a+(b*sqrt(sigma2))*(sqrt(1-rho^2))<0) return(1000)
  #code to solve butterfly arbitrage
  if (!all(butterfly(p,x) > 0)) return(1000)
  # sum of squares
  res=sum((totalIV-maturity*((v)^2))^2)
  
  return(res)
}

butterfly <- function(p,k){
  a= p[1]
  b=p[2]
  rho=p[3]
  m=p[4]
  sigma2=p[5]^2
  
  discr=sqrt((k-m)*(k-m)+sigma2)
  w <- a + b*(rho*(k-m)+discr)
  dw <- b*rho + b*(k-m)/discr
  d2w <- b*sigma2/(discr*discr*discr)
  gk <- 1-k*dw/w + dw*dw/4*(-1/w+k*k/(w*w)-4) + d2w/2
  return(gk)
}

library("DEoptim")


##typical criteria would be 
itermax <- 1e4
VTR <- 1e-4
l <- c(-20,0,-0.99,-20 ,1e-5)
u <- c(20,50,0.99,5,20)

#Using DEOptim to find the best parameters of p vector and store in pmatrix

x <- seq(from=-2, to =2, length.out = 100) #generate a sequence of log return to check butterfly abitrage
for (j in 1:48){ ###Calendar date here
  for (i in 1:10){ #monts to maturity
    
    i=1
    j=1
    maturity<-raw[((i*7)-6+(j-1)*70):((i*7)+(j-1)*70),2]/12
    v<-raw[((i*7)-6+(j-1)*70):((i*7)+(j-1)*70),11]
    k<- raw[((i*7)-6+(j-1)*70):((i*7)+(j-1)*70),12]
    outDEoptim<-DEoptim(vfunc,l,u,control=DEoptim.control(VTR=VTR,itermax=itermax))
    pmatrix[i+10*(j-1),]<-outDEoptim$optim$bestmem
  }}

#After get pmatrix, feed back to svi to get the estimated vol
svi<-function(a,b,rho,m,sigma2,maturity,k){
  kk=(k-m)^2+sigma2;
  totalIV= a+b*(rho * (k-m)+sqrt(kk));
  result <- sqrt(totalIV/(maturity))
  return(result)
}
svivol <- matrix(nrow=70,ncol=1)
cd <- 3
for (i in 1:10){
  for (j in 1:7){
    svivol[(i*7)-7+j] <- svi(pmatrix[i+(cd-1)*10,1],pmatrix[i+(cd-1)*10,2],pmatrix[i+(cd-1)*10,3],pmatrix[i+(cd-1)*10,4],pmatrix[i+(cd-1)*10,5]^2,
                             raw$Term[(i*7)-7+j+(cd-1)*70]/12,raw$moneyness[(i*7)-7+j+(cd-1)*70])
  }}

#Call price that has Vol from SVI
Cpfit<-vector(length=70)
for (i in 1:70){
  Cpfit[i] <- BSprice(1,raw$St[i],raw$Strike[i],svivol[i],raw$Cdiv[i],raw$Cinterest[i],
                      raw$Term[i]/12)
}
#plot Svi curve against actual observed IV
par(mfrow=c(2,5))
for (i in 1:10){
  plot(raw[((i*7)-6):(i*7),12],raw[((i*7)-6):(i*7),11],main=paste(raw$Term[(i*7 + - 6)], 'months to maturity', sep=' '),xlab="ln (K/F)",ylab="Implied volatility")
  lines(raw[((i*7)-6):(i*7),12],svivol[((i*7)-6):(i*7)])
  
  title("Jan 2006", outer=TRUE)
}

for (i in 1:10){
  plot(exp(raw[((i*7)-6):(i*7),12]),raw[((i*7)-6):(i*7),11],main=paste(raw$Term[(i*7 + - 6)], 'months to maturity', sep=' '),xlab="K/F",ylab="Implied volatility")
  lines(exp(raw[((i*7)-6):(i*7),12]),svivol[((i*7)-6):(i*7)])
  
  title("Jan 2006", outer=TRUE)
}
########################## QUESTION 5 ##########################
#plot call price surface to see if violate any NA
library(plot3Drgl)
par(mfrow=c(1,1))
scatter3Drgl((x=raw$Strike[1:70]),y=raw$Term[1:70],z=Cpfit,colvar=Cpfit,ticktype ="detailed", colkey= FALSE, facesets=FALSE)
#Check surface and then adjust restrictions up in vfunc

########################## QUESTION 6 ##########################
#plot actual vol in to rnd
#We have plotted this graph by using excel. The method is explained in the report.

########################## QUESTION 7 ##########################

#extract rnd using hessian function.


library(numDeriv)
#BSprice with SVI vol function inside, so the vol will change when stirke changes
Fitprice <- function(x){
  Mn <- log(x/Fd)
  kk=(Mn-m)^2+sigma2
  totalIV= a+b*(rho * (Mn-m)+sqrt(kk))
  vol <- sqrt(totalIV/(t))
  
  pc <- 1
  d1 <- (log(S / x) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  d2 <- d1 - vol * sqrt(t)
  
  fit <- pc * exp(-d * t) * S * 
    pnorm(pc * d1) - pc * x * exp(-r * t) * pnorm(pc * d2)
  return(fit)
}


nseg <- 1000 #number of segments
lret <- seq(from= -2, to=2, length.out = nseg) #sequences of log moneyness
dlret <- (lret[2] - lret[1])
q2 <- matrix(nrow=48*nseg, ncol=10)

for (k in 1:48){ ##date
  for (i in 1:10){
    Fd <- raw$Forward[(i*7)+(k-1)*70]
    x <- exp(lret)*Fd   #covernt log moneyness to strike price
    S <- raw$St[(i*7)+(k-1)*70]
    d <- raw$Cdiv[(i*7)+(k-1)*70]
    r <- raw$Cinterest[(i*7)+(k-1)*70]
    t <- raw$Term[(i*7)+(k-1)*70]/12
    a <- pmatrix[i+(k-1)*10,1]
    b <-pmatrix[i+(k-1)*10,2]
    rho <-pmatrix[i+(k-1)*10,3]
    m<-pmatrix[i+(k-1)*10,4]
    sigma2 <-pmatrix[i+(k-1)*10,5]^2
    
    for (j in 1:nseg){
      
      q2[j+(k-1)*nseg,i] <- hessian(Fitprice,x[j],method = 'Richardson')*exp((r-d)*t)
      if(q2[j+(k-1)*nseg,i]<0){q2[j+(k-1)*nseg,i]=0}  
      #even though we apply butterfly NA restriction(RND shouldn't be zero) there might be some extreme value because our log return in
      #DEoptim is more coarse than the simulated RND.
    }}}


#RND graph for each maturity
par(mfrow=c(2,5))
cd <- 2
for (j in 2:cd)
  for (i in 1:10){ 
    plot(lret,q2[(1+(j-1)*nseg):(nseg*j),i], type = "l", ylab="q", xlab ="log return")
  }
title("Risk neutral density", outer=TRUE)

Term <- c(1,3,6,12,24,36,48,60,84,120)

q2singles <- matrix(nrow = nseg, ncol = 10 )
q2singles <- q2[101:200,]

#RND surface

persp3Drgl(lret, Term, q2singles, xlab = "Ln Return", ylab = "Term", zlab = "RNP")


########################## QUESTION 9 ##########################
#extract return and variance
mu <- vector(length = 480)
varr <-vector(length = 480)
nu <- vector(length = 480)
for (j in 1:48){ # calendar dates
  for (i in 1:10){ # different term
    mu[i+(j-1)*10] <- (dlret*q2[(((j-1)*nseg)+1):(j*nseg),i])%*%(lret)
    varr[i+(j-1)*10] <- ((q2[(((j-1)*nseg)+1):(j*nseg),i])*dlret)%*%(lret^2)-mu[i+(j-1)*10]^2
    nu [i+(j-1)*10] <- varr[i+(j-1)*10]*12/(Term[i])
  }}


#draw graph

cl <- rainbow(11)

par(mfrow=c(1,2))
plot(Term,mu[1:10], type = "l", ylim = c(min(mu[1:120]),max(mu[1:120])), ylab = "mu" )
for (i in 1:11){
  lines(Term,mu[(1+10*i):((i+1)*10)],col = cl[i]) 
}
title("mu")
plot(Term,nu[1:10], type = "l", ylab = "nu", ylim = c(min(nu[1:120]),max(nu[1:120])))
for (i in 1:11){
  lines(Term,nu[(1+10*i):((i+1)*10)], col = cl[i]) 
}
title("nu")
###Q11 

library('quantmod')

t <- matrix(nrow=48,ncol=1)
for (i in 1:48){
  t[i] <- raw$ValuationDate[i*70]
}
t<- as.character(as.Date(t))
VIX<-get(getSymbols("^VIX",src="yahoo",from="2006-01-31",to="2009-12-31"))
VIX<-VIX[,4] #getting closing prices
# to get the values of VIX for the 48 dates from raw_data, use this
VIX<-as.numeric(VIX[which(as.character(index(VIX)) %in% t)]/1000)
# where t is the vector of the 48 dates
muplot <- vector(length = 48)
nuplot <- vector(length = 48)
for (i in 1:48){
  muplot[i] <- mu[i*10-9]*100
  nuplot[i] <- sqrt(nu[i*10-9]) #sqrt nu so nu and VIX are in similar scale
}

t <- as.Date(t)
par(mfrow=c(1,1))
plot(t,VIX, type ="l", ylim = c(min(muplot),max(VIX)), col="blue")
lines(t,muplot,type ="l",col="red")
lines(t,nuplot )
title("VIX(blue), SD (Black), and Mu (red)")

########################## QUESTION 13 ##########################

## here you will need to construct calculated variance from risk neutral density in the following format
## column : time to maturity; row: valuationdate
Q13 <- matrix(nrow=48, ncol = 10) 
for (i in 1:48){
  for (j in 1: 10){
    Q13 [i,j] <- nu[((i-1)*10+j)]
  }}
## if Q12 is the resulting matrix, to do the PCA, use
pz <- prcomp(Q13, tol = 0.1,na.rm=T)
summary(pz)

##################### QUESTION 15 & 16 ########################

###codify the Var equation from Q15 since it contains all the unknown params

q15varoptim <- function(p){
  
  t<-0
  T<-c(1,3,6,12,24,36,48,60,84,120)
  
  kappa <- p[1]
  eta <- p[2]
  thetabar <- p[3]
  sigma1 <- p[4]
  sigma2 <- p[5]
  rho <- p[6]
  
  T2 = (kappa - rho*sigma1 + exp(kappa*(t - T))*(rho*sigma1 + kappa*(-1 + rho*sigma1*(-t + T))))/kappa^2
  T3 = -(((-1 + exp(kappa*(t - T)))*rho*sigma1 + kappa*(kappa - rho*sigma1)*(t - T))/kappa^2)
  T4 = (exp(kappa*T)*(eta - kappa)^2*(kappa - rho*sigma1) + exp(eta*t - eta*T + kappa*T)*kappa^2*(eta - kappa + rho*sigma1) + exp(kappa*t)*eta*(kappa*(-eta + kappa) + rho*sigma1*(eta + kappa^2*(t - T) + kappa*(-2 - eta*t + eta*T))))/(exp(kappa*T)*eta*(eta - kappa)^2*kappa)
  
  Ttilde2 = (sigma1^2*(-exp(2*kappa*t) + exp(2*kappa*T) + 2*exp(kappa*(t + T))*kappa*(t - T)))/(2*exp(2*kappa*T)*kappa^3)
  Ttilde3 = -((3 - 4*exp(eta*(t - T)) + exp(2*eta*(t - T)))*kappa^6*sigma2^2 + eta*kappa^5*sigma2^2*(-1 + exp(2*eta*(t - T)) + 2*kappa*t - 2*kappa*T) + eta^6*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^5*kappa*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^4*kappa^2*(sigma1 - sigma2)*(sigma1 + sigma2)*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - 2*eta^2*kappa^4*sigma2^2*(2 - 2*(exp(eta*(t - T)) + exp(kappa*(t - T)) - exp((eta + kappa)*(t - T))) + kappa*t - kappa*T) + eta^3*kappa^3*(sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) + sigma2^2*(-1 + exp(2*kappa*(t - T)) - 2*kappa*t + 2*kappa*T)))/(4*eta^3*(eta - kappa)^2*kappa^3*(eta + kappa))
  Ttilde4 = (exp(-2*eta*t - 3*kappa*t - 7*eta*T - 6*kappa*T)*(2*exp(3*eta*t + 4*kappa*t + 6*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa^2*sigma2^2 - exp(4*eta*t + 3*kappa*t + 5*eta*T + 6*kappa*T)*(eta - 2*kappa)*kappa^4*sigma2^2 - exp(2*eta*t + 5*kappa*t + 7*eta*T + 4*kappa*T)*eta^3*((eta - kappa)^2*sigma1^2 - kappa^2*sigma2^2) + exp(2*eta*t + 3*kappa*t + 7*eta*T + 6*kappa*T)*(eta - 2*kappa)*(eta - kappa)^2*(eta^2*sigma1^2 + kappa^2*sigma2^2) + 2*exp(2*eta*t + 4*kappa*t + 7*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa*(-(kappa*sigma2^2) + eta*sigma1^2*(1 + (eta - kappa)*(t - T))) - 2*exp(3*(eta + kappa)*(t + 2*T))*eta*kappa^2*(-(eta*kappa*sigma1^2) + sigma2^2*(eta^2 + 2*kappa^3*(t - T) + eta*kappa*(-2 + eta*t - eta*T) + kappa^2*(2 - 3*eta*t + 3*eta*T)))))/(2*eta^3*(eta - 2*kappa)*(eta - kappa)^2*kappa^2)
  
  Tsigma2 <- T2 + Ttilde2 / 2
  Tthetabar <- T3 + Ttilde3 / 2
  Ttheta <- T4 + Ttilde4 / 2
  
  variances <- ((varassumed- thetabar) * Tsigma2) + (thetabar * Tthetabar) + ((thetaassumed - thetabar) * Ttheta)
  
  ###restrictions###
  if(kappa == 0) return(1000)
  if(eta == 0) return(1000)
  if(kappa == eta) return(1000)
  if(kappa == -eta) return(1000)
  if(eta == 2 * kappa) return(1000)
  if(!all(variances>0)) return(1000)
  
  
  ###sum squared (value to be minimized)###
  result = sum(((first10nu - variances))^2)
  
  return(result)
}

###map back the parameters to make sure they fit the variance equation###
q15varfit <- function(pvar,vt,thetat){
  
  t<-0
  T<-c(1,3,6,12,24,36,48,60,84,120)
  
  kappa <- pvar[1]
  eta <- pvar[2]
  thetabar <- pvar[3]
  sigma1 <- pvar[4]
  sigma2 <- pvar[5]
  rho <- pvar[6]
  
  
  T2 = (kappa - rho*sigma1 + exp(kappa*(t - T))*(rho*sigma1 + kappa*(-1 + rho*sigma1*(-t + T))))/kappa^2
  T3 = -(((-1 + exp(kappa*(t - T)))*rho*sigma1 + kappa*(kappa - rho*sigma1)*(t - T))/kappa^2)
  T4 = (exp(kappa*T)*(eta - kappa)^2*(kappa - rho*sigma1) + exp(eta*t - eta*T + kappa*T)*kappa^2*(eta - kappa + rho*sigma1) + exp(kappa*t)*eta*(kappa*(-eta + kappa) + rho*sigma1*(eta + kappa^2*(t - T) + kappa*(-2 - eta*t + eta*T))))/(exp(kappa*T)*eta*(eta - kappa)^2*kappa)
  
  Ttilde2 = (sigma1^2*(-exp(2*kappa*t) + exp(2*kappa*T) + 2*exp(kappa*(t + T))*kappa*(t - T)))/(2*exp(2*kappa*T)*kappa^3)
  Ttilde3 = -((3 - 4*exp(eta*(t - T)) + exp(2*eta*(t - T)))*kappa^6*sigma2^2 + eta*kappa^5*sigma2^2*(-1 + exp(2*eta*(t - T)) + 2*kappa*t - 2*kappa*T) + eta^6*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^5*kappa*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^4*kappa^2*(sigma1 - sigma2)*(sigma1 + sigma2)*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - 2*eta^2*kappa^4*sigma2^2*(2 - 2*(exp(eta*(t - T)) + exp(kappa*(t - T)) - exp((eta + kappa)*(t - T))) + kappa*t - kappa*T) + eta^3*kappa^3*(sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) + sigma2^2*(-1 + exp(2*kappa*(t - T)) - 2*kappa*t + 2*kappa*T)))/(4*eta^3*(eta - kappa)^2*kappa^3*(eta + kappa))
  Ttilde4 = (exp(-2*eta*t - 3*kappa*t - 7*eta*T - 6*kappa*T)*(2*exp(3*eta*t + 4*kappa*t + 6*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa^2*sigma2^2 - exp(4*eta*t + 3*kappa*t + 5*eta*T + 6*kappa*T)*(eta - 2*kappa)*kappa^4*sigma2^2 - exp(2*eta*t + 5*kappa*t + 7*eta*T + 4*kappa*T)*eta^3*((eta - kappa)^2*sigma1^2 - kappa^2*sigma2^2) + exp(2*eta*t + 3*kappa*t + 7*eta*T + 6*kappa*T)*(eta - 2*kappa)*(eta - kappa)^2*(eta^2*sigma1^2 + kappa^2*sigma2^2) + 2*exp(2*eta*t + 4*kappa*t + 7*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa*(-(kappa*sigma2^2) + eta*sigma1^2*(1 + (eta - kappa)*(t - T))) - 2*exp(3*(eta + kappa)*(t + 2*T))*eta*kappa^2*(-(eta*kappa*sigma1^2) + sigma2^2*(eta^2 + 2*kappa^3*(t - T) + eta*kappa*(-2 + eta*t - eta*T) + kappa^2*(2 - 3*eta*t + 3*eta*T)))))/(2*eta^3*(eta - 2*kappa)*(eta - kappa)^2*kappa^2)
  
  Tsigma2 <- T2 + Ttilde2 / 2
  Tthetabar <- T3 + Ttilde3 / 2
  Ttheta <- T4 + Ttilde4 / 2
  
  variances <- ((vt - thetabar) * Tsigma2) + (thetabar * Tthetabar) + ((thetat - thetabar) * Ttheta)
  
  return(variances)
}

###map back the parameters to make sure they fit the expected return equation###
q15returnsfit <- function(preturns,vt,thetat){
  
  t<-0
  T<-c(1,3,6,12,24,36,48,60,84,120)
  
  kappa <- preturns[1]
  eta <- preturns[2]
  thetabar <- preturns[3]
  sigma1 <- preturns[4]
  sigma2 <- preturns[5]
  rho <- preturns[6]
  
  T2 = (kappa - rho*sigma1 + exp(kappa*(t - T))*(rho*sigma1 + kappa*(-1 + rho*sigma1*(-t + T))))/kappa^2
  T3 = -(((-1 + exp(kappa*(t - T)))*rho*sigma1 + kappa*(kappa - rho*sigma1)*(t - T))/kappa^2)
  T4 = (exp(kappa*T)*(eta - kappa)^2*(kappa - rho*sigma1) + exp(eta*t - eta*T + kappa*T)*kappa^2*(eta - kappa + rho*sigma1) + exp(kappa*t)*eta*(kappa*(-eta + kappa) + rho*sigma1*(eta + kappa^2*(t - T) + kappa*(-2 - eta*t + eta*T))))/(exp(kappa*T)*eta*(eta - kappa)^2*kappa)
  
  Tsigma2tilde <- -T2 / 2
  Tthetabartilde <- -T3 / 2
  Tthetatilde <- -T4 / 2
  
  returns_expected <- ((vt - thetabar) * Tsigma2tilde) + (thetabar * Tthetabartilde) + ((thetat - thetabar) * Tthetatilde)
  
  return(returns_expected)
}

first10nu <- varr[1:10] #optimize variance curve by using first calendar date

library(DEoptim)

###DEoptim criteria###
itermax <- 1e4
VTR <- 1e-10
l <- c(1e-10,1e-10,0,0,0,-0.99)
u <- c(.01,.01,10,50,50,0.99)

###DEoptim###
t=0
pvar <- vector(length = 6)
varassumed <- 0.000005 #attempted different assumed values until we found a fit
thetaassumed <- 0.1
T <- c(1,3,6,12,24,36,48,60,84,120)

outDEoptim <- DEoptim(q15varoptim, lower = l, upper = u, control = DEoptim.control(VTR = VTR, itermax = itermax))
pvar <- outDEoptim$optim$bestmem

Variance_Structure <- q15varfit(pvar,varassumed,thetaassumed) #plot the fitted variance term structure against the observed nu
plot(Variance_Structure, type= "l")
lines(first10nu)

Returns_Structure <- q15returnsfit(pvar,varassumed,thetaassumed) #plot the returns term structure against the observed mu
plot(Returns_Structure, type = "l")
lines(first10mu)

##################### QUESTION 17 ########################


###Monte Carlo###
kappa <- pvar[1]
eta <- pvar[2]
thetabar <- pvar[3] 
sigma1 <- pvar[4]
sigma2 <- pvar[5]
rho <- pvar[6]

###Start Monte Carlo###
paths <- 101
times <- 48
Delta <- 1/12
theta <- matrix(nrow = paths, ncol = times) #initialize theta, v, and x for Monte Carlo
v <- matrix(nrow = paths, ncol = times)
x <- matrix(nrow = paths, ncol = times)
theta[,1] <- thetaassumed
v[,1] <- varassumed
x[,1] <- 0
Z1 <- matrix(nrow = paths, ncol = times) #initialize Z1-Z3 to store random distributions
Z2 <- matrix(nrow = paths, ncol = times)
Z3 <- matrix(nrow = paths, ncol = times)
for (i in 1:paths){
  for (j in 1:times){
    Z1[i,j] <- rnorm(1, mean = 0) * Delta #generate a random normal distribution
    Z2[i,j] <- rnorm(1, mean = 0) * Delta
    Z3[i,j] <- rnorm(1, mean = 0) * Delta
  }
}

for (i in 2:times){                       #use Milstein scheme to discretize SDE for theta
  theta[,i] = (sqrt(theta[,i-1]) + sigma2/2 * sqrt(Delta) * Z3[,i-1])^2 
  - eta * (theta[,i-1] - thetabar) * Delta - sigma2^2/4 * Delta
}


for (i in 2:times){                       #use Milstein scheme to discretize SDE for v
  v[,i] = (sqrt(v[,i-1]) + sigma1/2 * sqrt(Delta) * Z2[,i-1])^2 
  - kappa * (v[,i-1] - theta[,i-1]) * Delta - sigma1^2/4 * Delta
}


for (i in 2:times){  #create return time series using simple discretization of SDE
  x[,i] = x[,i-1] - v[,i-1]/2 * Delta + sqrt(Delta * v[,i-1]) * Z1[,i-1]
}


t <- vector(length = 48)  #define t as calendar dates for graph label  
for (i in 1:48){
  t[i] <- raw$ValuationDate[i*70]
}
t<- as.Date(t)

par(mfrow=c(2,2))
plot(t,x[1,], type = "l", ylim=c((min(x)),max(x))) #plot path of x across calendar dates
for (i in 2:paths){
  lines(t,x[i,])
}


mumatrices <- list()  #Create lists to store the mu and nu for each calendar date
numatrices <- list()
matrixxx <- list()
tau <- Term / 12

for (j in 1:48){      #Create a maxtrix stored with nu values for each calendar date and term
  for (i in 1:paths){
    
    thetasim <- theta[i,j]
    vsim <- v[i,j]
    matrixxx[[i]] <- q15varfit(pvar[1:6], thetasim, vsim) / tau
  }
  
  do.call(rbind, matrixxx)   #convert list to matrix
  numatrices[[j]] <- matrixxx    #store matrix in list of matrices
}

for (j in 1:48){     #Perform same as above for mu
  for (i in 1:paths){
    
    thetasim <- theta[i,j]
    vsim <- v[i,j]
    matrixxx[[i]] <- q15returnsfit(pvar[1:6], thetasim, vsim)
    
  }
  do.call(rbind, matrixxx)
  mumatrices[[j]] <- matrixxx
}

par(mfrow = c(4,6))
for (j in 1:48){        #plot simulated nu for all calendar dates
  plot(Term, numatrices[[j]][[1]], 'l', ylim=c(0,4))
  for (i in 2:paths){
    lines(Term, numatrices[[j]][[i]])
  }
}

for (j in 1:48){        #plot simulated mu for all calendar dates
  plot(Term, mumatrices[[j]][[1]], 'l',ylim=c(-5,0))
  for (i in 2:paths){
    lines(Term, mumatrices[[j]][[i]])
  }
}

for (i in 1:48){        #plot actual mu for all calendar dates
  plot(Term, mu[(i*10-9):(i*10)],type = "l")
}

for (i in 1:48){        #plot actual nu for all calendar dates
  plot(Term, nu[(i*10-9):(i*10)],type = "l")
}


##################### QUESTION 18 ########################

z01 <- vector(length = 48)
z05 <- vector(length = 48)
match01 <- vector(length = 48)
match05 <- vector(length = 48)
for (i in 1:48){      #find the 0.01 and 0.05 quantile value
  z01[i] <- quantile(x[,i],.01)
  z05[i] <- quantile(x[,i],.05)
  match01[i] <- match(z01[i],x[,i]) #identify the index associated with each quantile value
  match05[i] <- match(z05[i],x[,i])
}

theta01 <- vector(length = 48) 
v01 <-vector(length = 48)
theta05 <-vector(length = 48)
v05 <-vector(length = 48)

for (i in 1:48){           #extract theta and v
  theta01[i] <- theta[match01[i],i]
  v01[i] <- v[match01[i],i]
  theta05[i] <- theta[match05[i],i]
  v05[i] <- v[match05[i],i]
}

vars01 <- matrix(nrow = 48, ncol = 10)
vars05 <- matrix(nrow = 48, ncol = 10)
returns01 <- matrix(nrow = 48, ncol = 10)
returns05 <- matrix(nrow = 48, ncol = 10)

for (i in 1:48){        #calculate term structure of variances and returns from Q15 model
  vars01[i,] <- q15varfit(pvar, v01[i], theta01[i])/tau
  vars05[i,] <- q15varfit(pvar, v05[i], theta05[i])/tau
  returns01[i,] <- q15returnsfit(pvar, v01[i], theta01[i])
  returns05[i,] <- q15returnsfit(pvar, v05[i], theta05[i])
}


##################### QUESTION 19 ########################

par(mfrow=c(2,2))   
plot(Term,vars01[1,],type="l", ylim=(c(min(vars01),max(vars01)))) #plot term structure of 0.01 quantile variances
for (i in 2:48){
  lines(Term, vars01[i,])
}

plot(Term,vars05[1,],type="l", ylim=(c(min(vars05),max(vars05)))) #plot term structure of 0.05 quantile variances
for (i in 2:48){
  lines(Term, vars05[i,])
}

plot(Term,returns01[1,],type="l", ylim=(c(min(returns01),max(returns01)))) #plot term structure of 0.01 quantile returns
for (i in 2:48){
  lines(Term, returns01[i,])
}

plot(Term,returns05[1,],type="l", ylim=(c(min(returns05),max(returns05)))) #plot term structure of 0.05 quantile returns
for (i in 2:48){
  lines(Term, returns05[i,])
}

##################### QUESTION 20 ########################

for (i in 1:48){
  muplot[i] <- mu[i*10-9]*100
  nuplot[i] <- sqrt(nu[i*10-9])
}

t <- as.Date(t)

par(mfrow=c(1,1))

plot(t,VIX, type ="l", ylim = c(min(muplot),max(VIX)), col="blue") #plot the VIX for the period

for (i in 1:10){    #plot the sqrt(var) for each term against the VIX (for both 0.01 and 0.05 quantiles)
  lines(t,sqrt(vars01[,i]),col="green")
  lines(t, sqrt(vars05[,i]))
}  
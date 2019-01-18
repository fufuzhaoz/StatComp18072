## ---- echo=TRUE----------------------------------------------------------
cdf <- function(x,p,q) {
  n<-length(x)
  res<-numeric(n)
  for(i in 1:n){
    res[i]<-integrate(f, lower=-Inf, upper=x[i],
                      rel.tol=.Machine$double.eps^0.25,
                      p=p,q=q)$value
  }
  return(res)
}



## ---- echo=TRUE----------------------------------------------------------
p <- 0
q <- 1
x<-seq(-10,10,0.1)
f<-function(x,p,q){      #density function
   1/(q*pi*(1+((x-p)/q)^2))
 }
plot(x,cdf(x,p,q))

## ---- echo=TRUE----------------------------------------------------------
mytable <-function(x,y){
  x=as.factor(x)
  y=as.factor(y)
  xi=as.integer(x)
  xm=as.integer(max(xi))
  yi=(as.integer(y)-1L)*xm
  ym=as.integer(max(yi)+xm)
  matrix(.Internal(tabulate(xi+yi,ym)),xm)
}

## ---- echo=TRUE----------------------------------------------------------
library(microbenchmark)
x<-c(1,2,2,3,4)
y<-c(1,2,3,3,4)
mytable(x,y)
microbenchmark(table(x,y),mytable(x,y))

## ---- echo=TRUE----------------------------------------------------------
cramer <-function(x,y){ #compute the Cramer-von Mises statistic
  n <- length(x)
  m <- length(y)
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  W2 <- ((m*n)/(m+n)^2)*
    (sum((Fn(x)-Gm(x))^2)+sum((Fn(y)-Gm(y))^2))
  return(W2)
}

## ---- echo=TRUE----------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
w2.0 <- cramer(x,y)
print(w2.0)

## ---- echo=TRUE----------------------------------------------------------
library(datasets)
library(xtable)
data(tli)
tli.table <- xtable(tli[1:20,])
digits(tli.table)[c(2,6)] <- 0
print(tli[1:20,])

## ----ehco=TRUE-----------------------------------------------------------
    drnd <- function(x, p, n) {
    z <- NULL
    ps <- cumsum(p)
    r <- runif(n)
    for (i in 1:n) z <- c(z, x[which(r[i] <= ps)[1]])
    return(z)
}
x <- c(0,1, 2, 3, 4);
p <- c(0.1, 0.2, 0.2, 0.2, 0.3);
out <- drnd(x, p, 1000);
out
set.seed(1)
table(out)/10000
mean(out)
var(out)

## ---- echo=TRUE----------------------------------------------------------
n <- 1000; 
k <- 0 ;
j<- 0 ;
y <- numeric(n);
while (k < n) { 
  u <- runif(1) ;
  j <- j + 1 ;
  x <- runif(1)  ;
  if (x^2 * (1-x) > u) {
    k <- k + 1;
  y[k] <- x ;} }
j
set.seed(1)

## ---- echo=TRUE----------------------------------------------------------
p <- seq(.1, .9, .1);
Qhat <- quantile(y, p) ;
Q <- qbeta(p, 3, 2) 
se <- sqrt(p * (1-p) / (n * dbeta(Q, 2, 2)))
 round(rbind(Qhat, Q, se), 3) 
 par(mfcol=c(1,2)) 
 hist(Qhat, prob=TRUE) 
 hist(Q, prob=TRUE) 
 par(mfcol=c(1,1)) 

## ----echo=TRUE-----------------------------------------------------------
n <- 1000;
r <- 4 ;
beta <- 2 ;
lambda <- rgamma(n, r, beta) ;
x <- rexp(rgamma(n,lambda));
x
set.seed(1)

## ---- echo=TRUE----------------------------------------------------------
x <- seq(.1, 0.9, length = 9)
m <- 10000 
u <- runif(m) 
cdf <- numeric(length(x))
for (i in 1:length(x)) { 
  g <- 30*(x[i]^3)*(u^2)*((1-x[i]*u)^2)
  cdf[i] <- mean(g)
}
cdf

## ----echo=TRUE-----------------------------------------------------------
Phi <- pbeta(x,3,3)
print(round(rbind(x, cdf, Phi), 3)) 

## ---- echo=TRUE----------------------------------------------------------

MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
  u <- runif(R/2) 
  if (!antithetic) v <- runif(R/2) else 
    v <- 1 - u 
  u <- c(u, v) 
  cdf <- numeric(length(x)) 
  for (i in 1:length(x)) {
    g <- x[i]*u*exp(-(x[i]*u)^2 / 2) 
    cdf[i] <- mean(g)
  }
  cdf
}


## ---- echo=TRUE----------------------------------------------------------

MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
  u <- runif(R/2) 
  m <- 10000
  MC1 <- MC2 <- numeric(m)
  if (!antithetic) v <- runif(R/2) else 
    v <- 1 - u 
  u <- c(u, v) 
  cdf <- numeric(length(x)) 
  for (i in 1:length(x)) {
    g <- (x[i]^2)*u*exp(-(x[i]*u)^2 / 2) 
    cdf[i] <- mean(g)
  }
  cdf
}
x <- seq(1,3, length=10) 
set.seed(123)
MC1 <- MC.Phi(x, anti = TRUE) 
MC2 <- MC.Phi(x,anti=FALSE)
c(var(MC1),var(MC2))

## ---- echo=TRUE----------------------------------------------------------
print(c((var(MC2)-var(MC1))/var(MC2)))

## ---- echo=TRUE----------------------------------------------------------
m <- 10000 
theta.hat <- se <- numeric(5)
g <- function(x) { exp(-x^2 + log((x^2)/sqrt(2*pi)))}
x <- rexp(m) 
fg <- sqrt(2*pi)* g(x)/exp(x-(x^2)/2) 
theta.hat <- mean(fg) 
se[1] <- sd(fg)
se[1]

## ---- echo=TRUE----------------------------------------------------------
m <- 10000 
theta.hat <- se <- numeric(5)
g <- function(x) { exp(-x^2 + log((x^2)/sqrt(2*pi)))}
x <- runif(m) 
fg<- g(x)/(sqrt(2*pi)*(x^2))
theta.hat[2] <- mean(fg) 
se[2] <- sd(fg)
se[2]

## ---- echo=TRUE----------------------------------------------------------
f1 <- function(x) { exp(-x^2 + log((x^2)/sqrt(2*pi)))}
curve(f1,1,5)


## ---- echo=TRUE----------------------------------------------------------
f2 <- function(x){1/(sqrt(2*pi)*(x^2))}
curve(f2,1,5)


## ---- echo=TRUE----------------------------------------------------------
m <- 10000 
theta.hat <- se <- numeric(5)
g <- function(x) { exp(-x^2 + log((x^2)/sqrt(2*pi)))}
x <- runif(m) 
fg<- g(x)/(sqrt(2*pi)*(x^2))
theta.hat[2] <- mean(fg) 
se[2] <- sd(fg)
se[2]

## ---- echo=TRUE----------------------------------------------------------
library(bootstrap)
data(law,package = "bootstrap")

 n <- nrow(law) 
 y <- law$LSAT
 z <- law$GPA
 theta.hat <- cor(y,z)
 
 theta.jack <-numeric(n)
for (i in 1:n) {
  theta.jack[i] <- cor(y[-i],z[-i])
}
bias <- (n - 1)*(mean(theta.jack) - theta.hat)
se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2)) 
round(c(original=theta.hat,bias, se),3)


## ---- echo=TRUE----------------------------------------------------------
library(boot) 
theta.boot <- function(x,i) { 
  mean(x[i]) 
   }
x <- c(3,5,7,18,43,85,91,98,100,130,230,487)
boot.obj <- boot(x, statistic = theta.boot, R = 2000)
print(boot.obj)

## ---- echo=TRUE----------------------------------------------------------
print(boot.ci(boot.obj, type = c("basic", "norm", "perc")))

## ---- echo=TRUE----------------------------------------------------------
library(boot) 
n <- length(x)
B <- 2000

x <- c(3,5,7,18,43,85,91,98,100,130,230,487)
theta.hat <- mean(x) 
theta.b <- numeric(B)

for (b in 1:B) {
  i <- sample(1:n, size = n, replace = TRUE) 
  y <- x[i] 
  theta.b[b] <- mean(y) 
}

boot.BCa <- function(x, th0, th, stat, conf = .95) {
  x <- as.matrix(x) 
  n <- nrow(x) 
  N <- 1:n 
  alpha <- (1 + c(-conf, conf))/2 
  zalpha <- qnorm(alpha)
  
 z0 <- qnorm(sum(th < th0) / length(th))
 th.jack <- numeric(n)
 for (i in 1:n) { 
   J <- N[1:(n-1)] 
   th.jack[i] <- stat(x[-i, ], J)
   } 
 L <- mean(th.jack) - th.jack 
 a <- sum(L^3)/(6 * sum(L^2)^1.5)
 
 adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
 limits <- quantile(th, adj.alpha, type=6)
 return(list("est"=th0, "BCa"=limits))
}

  stat <- function(x,i) {
  mean(x[i])
} 
 
boot.BCa(x, th0 = theta.hat, th = theta.b, stat = stat)

## ---- echo=TRUE----------------------------------------------------------
library(bootstrap)
data(scor,package = "bootstrap")
n <- nrow(law)
v <- cov(scor)
E <- eigen(v)$value
theta.hat <- E[1]/sum(E)


theta.jack <- numeric(n)
for (i in 1:n) {
  scorstar <- scor[-i,]
  f <- eigen(cov(scorstar))$value
  theta.jack[i] <- f[1]/sum(f)
}
bias <- (n - 1)*(mean(theta.jack) - theta.hat)
se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2)) 
round(c(original=theta.hat,bias, se),3)

## ---- echo=TRUE----------------------------------------------------------
library(DAAG)
data(ironslag,package = "DAAG")
set.seed(1)
n <- length(ironslag$magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
e11 <- e21 <- e31 <- e41 <- numeric(n-1)
e12 <- e22 <- e32 <- e42 <- numeric(n-1)

for (k in 1:(n-1)) { 
  y <- ironslag$magnetic[c(-k,-(k+1))] 
  x <- ironslag$chemical[c(-k,-(k+1))]

  J1 <- lm(y ~ x)
  yhat11 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[k] 
  yhat12 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[k+1]
  e11[k] <- ironslag$magnetic[k] - yhat11
  e12[k+1] <- ironslag$magnetic[k+1] - yhat12
  e1[k] <- e11[k]+e12[k+1]
  
  J2 <- lm(y ~ x + I(x^2)) 
  yhat21 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[k] + 
    J2$coef[3] * ironslag$chemical[k]^2
  yhat22 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[k+1] + 
    J2$coef[3] * ironslag$chemical[k+1]^2 
  e21[k] <- ironslag$magnetic[k] - yhat21
  e22[k+1] <- ironslag$magnetic[k+1] - yhat22
  e2[k] <- e21[k]+e22[k+1]
  
  J3 <- lm(log(y) ~ x)
  logyhat31 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[k] 
  logyhat32 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[k+1]
  yhat31 <- exp(logyhat31) 
  yhat32 <- exp(logyhat32) 
  e31[k] <- ironslag$magnetic[k] - yhat31
  e32[k+1] <- ironslag$magnetic[k+1] - yhat32
  e3[k] <- e31[k]+e32[k+1]
  
  J4 <- lm(log(y) ~ log(x)) 
  logyhat41 <- J4$coef[1] + J4$coef[2] * log(ironslag$chemical[k])
  logyhat42 <- J4$coef[1] + J4$coef[2] * log(ironslag$chemical[k+1])
  yhat41 <- exp(logyhat41)
  yhat42 <- exp(logyhat42)
  e41[k] <- ironslag$magnetic[k] - yhat41
  e42[k+1] <- ironslag$magnetic[k+1] - yhat42
  e4[k] <- e41[k]+e42[k+1]
}
 c(mean(e1^2), mean(e2^2),mean(e3^2), mean(e4^2)) 

## ---- echo=TRUE----------------------------------------------------------
library(DAAG)
data(ironslag,package = "DAAG")
set.seed(1)
z <- sample(ironslag,length(ironslag)) #disrupt the sample
n <- length(z$magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
e11 <- e21 <- e31 <- e41 <- numeric(n-1)
e12 <- e22 <- e32 <- e42 <- numeric(n-1)

for (k in 1:(n-1)) { 
  y <- z$magnetic[c(-k,-(k+1))] 
  x <- z$chemical[c(-k,-(k+1))]

  J1 <- lm(y ~ x)
  yhat11 <- J1$coef[1] + J1$coef[2] * z$chemical[k] 
  yhat12 <- J1$coef[1] + J1$coef[2] * z$chemical[k+1]
  e11[k] <- z$magnetic[k] - yhat11
  e12[k+1] <- z$magnetic[k+1] - yhat12
  e1[k] <- e11[k]+e12[k+1]
  
  J2 <- lm(y ~ x + I(x^2)) 
  yhat21 <- J2$coef[1] + J2$coef[2] * z$chemical[k] + 
    J2$coef[3] * z$chemical[k]^2
  yhat22 <- J2$coef[1] + J2$coef[2] * z$chemical[k+1] + 
    J2$coef[3] * z$chemical[k+1]^2 
  e21[k] <- z$magnetic[k] - yhat21
  e22[k+1] <- z$magnetic[k+1] - yhat22
  e2[k] <- e21[k]+e22[k+1]
  
  J3 <- lm(log(y) ~ x)
  logyhat31 <- J3$coef[1] + J3$coef[2] * z$chemical[k] 
  logyhat32 <- J3$coef[1] + J3$coef[2] * z$chemical[k+1]
  yhat31 <- exp(logyhat31) 
  yhat32 <- exp(logyhat32) 
  e31[k] <- z$magnetic[k] - yhat31
  e32[k+1] <- z$magnetic[k+1] - yhat32
  e3[k] <- e31[k]+e32[k+1]
  
  J4 <- lm(log(y) ~ log(x)) 
  logyhat41 <- J4$coef[1] + J4$coef[2] * log(z$chemical[k])
  logyhat42 <- J4$coef[1] + J4$coef[2] * log(z$chemical[k+1])
  yhat41 <- exp(logyhat41)
  yhat42 <- exp(logyhat42)
  e41[k] <- z$magnetic[k] - yhat41
  e42[k+1] <- z$magnetic[k+1] - yhat42
  e4[k] <- e41[k]+e42[k+1]
}
 c(mean(e1^2), mean(e2^2),mean(e3^2), mean(e4^2)) 

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)

attach(chickwts) 
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts) 


cramer <-function(x,y){ #compute the Cramer-von Mises statistic
  n <- length(x)
  m <- length(y)
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  W2 <- ((m*n)/(m+n)^2)*
      (sum((Fn(x)-Gm(x))^2)+sum((Fn(y)-Gm(y))^2))
  W2
}

 
R <- 999 
z <- c(x, y)
n <- length(x)
m <- length(y)
K <- 1:(m+n) 
w2 <- numeric(R) 

w2.0 <- cramer(x,y) 


for (i in 1:R) { 
  k <- sample(K, size = n, replace = FALSE) 
  x1 <- z[k]
  y1 <- z[-k] 
  w2[i] <- cramer(x1,y1)
  }
p <- mean(c(w2.0, w2) >= w2.0) 

print(p)

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
f <- function(x, theta,eta) {
  stopifnot(theta > 0)
  return(1/(theta*pi*(1+((x-eta)/theta)^2))) 
} 
xt <- x[i-1] 
y <-  xt+runif(1,-1,1)
m <- 10000
theta <- 1
eta <- 0
x <- numeric(m) 
x[1] <- runif(1,-1,1) 
k <- 0
u <- runif(m)

for (i in 2:m) {
  xt <- x[i-1]
  y <- xt+runif(1,min=-1,max=1)
  num <- f(y, theta,eta)*dnorm(xt,mean=y,sd=1) 
  den <- f(xt, theta,eta)*dnorm(y,mean=xt,sd=1)
  if (u[i] <= num/den) x[i] <- y else {
    x[i] <- xt
    k <- k+1         #y is rejected 
  } 
}
print(k) 

## ---- echo=TRUE----------------------------------------------------------
index <- 1500:2000
y1 <- x[index] 
plot(index, y1, type="l", main="", ylab="x") 

## ---- echo=TRUE----------------------------------------------------------
gene<- quantile(x[-(1:1000)],seq(0.1,0.9,0.1))
empi<-qt(seq(0.1,0.9,0.1),1)
rbind(gene,empi)

## ---- echo=TRUE----------------------------------------------------------
b <- 1001     #discard the burnin sample 
y <- x[b:m] 
probs <- seq(0.1,0.9,0.01)
QR <- qt(probs,df=1)
Q <- quantile(x, probs)

qqplot(QR, Q, main="", 
       xlab="Cauchy Quantiles", ylab="Sample Quantiles")

hist(y, breaks="scott", main="", xlab="", freq=FALSE) 
lines(QR, f(QR, 1,0))

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)

w <- .25   #width of the uniform support set 
m <- 5000  #length of the chain 
burn <- 1000  #burn-in time 
total <- 197 
x <- numeric(m)  #the chain

## ---- echo=TRUE----------------------------------------------------------
group <- c(125,20,18,34)

prob <- function(y, group) { 
  # computes (without the constant) the target density 
  if (y < 0 || y >= 1) 
    return (0) 
  return(((2+y)/4)^group[1] * 
           ((1-y)/4)^group[2] * ((1-y)/4)^group[3] *
           ((y/4)^group[4]))
} 

u <- runif(m)          #for accept/reject step 
v <- runif(m, -w, w)   #proposal distribution
x[1] <- 0.5
for (i in 2:m) {
  y <- x[i-1] + v[i]
  if (u[i] <= prob(y, group) / prob(x[i-1], group)) 
    x[i] <- y else 
      x[i] <- x[i-1]
}
print(group) 
print(round(group/total, 3))

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
 xb <- x[(burn+1):m] 
 print(mean(xb))    

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)

Gelman.Rubin <- function(psi) { 
  # psi[i,j] is the statistic psi(X[i,1:j]) 
  # for chain in i-th row of X 
  
  psi <- as.matrix(psi) 
  n <- ncol(psi) 
  k <- nrow(psi)

  psi.means <- rowMeans(psi)     #row means 
  B <- n * var(psi.means)        #between variance est. 
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est. 
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est. 
  r.hat <- v.hat / W             #G-R statistic 
  return(r.hat) }

prob <- function(y,group) {
  # computes (without the constant) the target density)
  if (y < 0 || y >= 1) 
    return (0) 
  return(((2+y)/4)^group[1] * 
           ((1-y)/4)^group[2] * ((1-y)/4)^group[3] *
           ((y/4)^group[4]))
} 

unif.chain <- function(group,m,X1) { 
  x <- numeric(m)
  x[1] <- X1 
  u <- runif(m)
  v <- runif(m, -w, w)   #proposal distribution
  
  
for (i in 2:m) {
  y <- x[i-1] + v[i]
  if (u[i] <= prob(y, group) / prob(x[i-1], group)) 
    x[i] <- y else 
      x[i] <- x[i-1]
}
  return(x)
}

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
K <-c(4:25,100,500,1000)
solution <- numeric(length(K))
for (i in 1:length(K)) {
  k <- K[i]
  model <- function(a){
  pt(sqrt(a^2*(k-1)/(k-a^2)),df=k-1)-
  pt(sqrt(a^2*k/(k+1-a^2)),df=k)
     }
  solution[i] <- uniroot(model, c(1,2))$root
}
cbind(K,solution)

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
f<-function(x,eta,theta){      #density function
  1/(theta*pi*(1+((x-eta)/theta)^2))
}
cdf <- function(x,eta,theta) {
    n<-length(x)
    res<-numeric(n)
  for(i in 1:n){
     res[i]<-integrate(f, lower=-Inf, upper=x[i],
     rel.tol=.Machine$double.eps^0.25,
     theta=theta,eta=eta)$value
  }
    return(res)
}

## ---- echo=TRUE----------------------------------------------------------
#standard Cauchy ditribution
eta <- 0
theta <- 1
x<-seq(-10,10,0.1)
plot(x,cdf(x,eta,theta), lwd = 0.5,ylab="cdf",col="blue",main="cdf of Cauchy(eta=0, theta=1)")
curve(pcauchy(x,location = 0, scale = 1),add=TRUE,lwd = 2,col="green")

## ---- echo=TRUE----------------------------------------------------------
set.seed(123)
nA<-28 
nB<-24 
nOO<-41
nAB<-70  #observed data

tol <- .Machine$double.eps^0.5
N <- 10000 #max. number of iterations
p0<-0.2
q0<-0.2   #initial estimates
p<-p0
q<-q0
L<-numeric(N) #Record the maximum likelihood values in M-steps
k<-1  #Record the circle times

for (i in 1:N) {
a<-p/(2-p-2*q)
b<-q/(2-2*p-q)
Hp=nA*(1+a)+nAB
Hq=nB*(1+b)+nAB
Hpq=nA*(1-a)+nB*(1-b)+2*nOO

L[i]<-Hp*log(p)+Hpq*log(1-p-q)+Hq*log(q) #Record the    maximum likelihood values in M-steps

x<-p
y<-q                #storage previous estimates
p<-Hp/(Hp+Hq+Hpq)
q<-Hq/(Hp+Hq+Hpq)
k<-k+1

if (abs(p-x)/x<tol && abs(q-y)/y<tol) break
p.hat<-p
q.hat<-q
L<-L[1:k]
}
p.hat
q.hat
L

## ---- echo=TRUE----------------------------------------------------------
d<-L[c(2:length(L))]-L[c(1:(length(L)-1))]
plot(1:9, d, main = "max", xlab = "Order", ylab = "Height")

## ------------------------------------------------------------------------
library("MASS")
m <- length(Cars93)
# Create a data frame from the main data set.
car.data1 <- rbind(Cars93$AirBags, Cars93$Type)

# Create a table with the needed variables.
car.data2 <- table(Cars93$AirBags, Cars93$Type) 


# Perform the Chi-Square test.


## ------------------------------------------------------------------------
library(microbenchmark)

sample_rows <- function(df, i) sample.int(nrow(df), i, replace = TRUE)

boot_cor1 <- function(df, i) { sub <- df[sample_rows(df, i), , drop = FALSE] 
table(sub$x, sub$y) }

boot_cor2 <- function(df, i ) { idx <- sample_rows(df, i)
table(df$x[idx], df$y[idx]) }
df <- data.frame(x = runif(100), y = runif(100))
microbenchmark( 
  boot_cor1(df, 10),
  boot_cor2(df, 10)
) 


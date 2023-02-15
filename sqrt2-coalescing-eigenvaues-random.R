
###########################################
# Giovanni Motta
# Copyright © 2023
# random (deterministic) matrix-valued functions
# sqrt2-smoothing in the case of coalescing-eigenvalues 
###########################################

rm(list=ls())
#library(ggplot2)
library(grid)
#library(dplyr)
#library(knitr)
#library(data.table)
#library(pracma)
library(expm)
library(Matrix)
library(MASS)

#------some graphical parameters------

setpar=function(...) {
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(0.5,0.5,0.5),
  cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)}
# outer=oma(bottom, left,top,right)
# inner=mar(bottom, left,top,right)

N=4
T=1200
u=c(1:T)/T
r=2
L=array(0,c(N,r,T))

TT=T/2
v=seq(0,1,length.out=TT)
LL=array(0,c(N,r,TT))


tt=1:TT
c=15
c2=22
cc=(TT/2 -c+1):(TT/2 +c)
hh=.04
##############################
# loadings
##############################

fr=rep(0.75,N)
ph=rep(0.2,N)

nn=seq(1,N,2)
for (n in nn){
for (p in 1:r){ #for all t and all p, sum(L.sim[,p,t]^2) = 1 
  L[n,p,] =((sin(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
  LL[n,p,]=((sin(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
  L[n+1,p,] =((cos(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
  LL[n+1,p,]=((cos(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
}}

for (h in 1:3){
rk=0
for (tt in 1:T){
M=t(L[,,tt])%*%L[,,tt]
rk=c(rk,rankMatrix(M))
if (rankMatrix(M)==r){
M=sqrtm(M)
M=solve(M)
L[,,tt]=L[,,tt]%*%M}}
}
#L=L*sqrt(N)

NN=N
mi=min(L[1:NN,,])
ma=max(L[1:NN,,])
x11()
setpar(mfrow=c(r,NN),oma=c(1.2,3.4,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
for (j in 1:NN){
#if (j==1){
#plot(uu,L.sim[j,k,26:T],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none")}
#else {
plot(u,L[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt='n',xaxt="none")
text(0.5,0.8*ma,cex=2,bquote(Lambda[.(j)*','*.(k)]))
}}

dev.off()

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=1.6)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)
}
else{
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt="none",xaxt="none",cex.axis=1.6)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
text(0.5,0.7*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,cex.axis=1.6)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt="none",cex.axis=1.6)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
text(0.5,0.7*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}
dev.off()


M=array(0,c(r,r,T))
for (tt in 1:T){
M[,,tt]=t(L[,,tt])%*%L[,,tt]}

x11()
setpar(mfrow=c(r,r),oma=c(0.5,1.4,1,0.5),mar=c(.8,.3,1,0.8))
for (k in 1:r){
for (j in 1:r){
plot(u,M[j,k,],type="l",ylim=c(-2,2),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,1.3,cex=2,bquote(G[.(j)*','*.(k)]^P))
}}

dev.off()

##############################
# Eigenvalues = Var(Factors)
##############################

G.f=matrix(0,T,r)
G.f[,1]= .5*(1+cos(pi*u))
G.f[,2]= .5*(1-cos(pi*u))

x11()
plot(u,G.f[,1],type="l",lwd=2,ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1))
text(1.01,0,cex=2,bquote(lambda[.(1)]))
lines(u,G.f[,2],col="grey",lwd=2)
text(1.01,1,cex=2,bquote(lambda[.(2)]))

dev.off()

x11()
plot(u,apply(G.f,FUN=max,MARGIN=1),type="l",lwd=2,ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1))
text(1.01,0.98,cex=2,bquote(hat(lambda)[.(1)]))
lines(u,apply(G.f,FUN=min,MARGIN=1),col="grey",lwd=2)
text(1.01,0.02,cex=2,bquote(hat(lambda)[.(2)]))

dev.off()


##############################
# Factors
##############################

e=mvrnorm(n = T, mu=rep(0,r), Sigma=diag(r), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
e=t(e)
F=e

for (t in 1:T){
F[,t]=sqrt(diag(G.f[t,]))%*%matrix(e[,t],r,1)
}

# B L T R
x11()
setpar(mfrow=c(r,1),oma=c(1.5,2,0.5,0),mar=c(1.5,0.2,.3,0))
j=1
plot(1:T,F[j,],type="l",lwd=1,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(min(F),max(F)),xaxt='n')
text(T/2,max(F[j,]),cex=2,bquote(F[.(j)]))
abline(h=0)
j=2
plot(1:T,F[j,],type="l",lwd=1,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(min(F),max(F)))
text(T/2,max(F[j,]),cex=2,bquote(F[.(j)]))
abline(h=0)
dev.off()

##############################
# Common components
##############################

X=matrix(0,N,T)

for (t in 1:T){
X[,t]=L[,,t]%*%F[,t]
}

x11()
setpar(mfrow=c(N,1),oma=c(1.5,2,0.5,0),mar=c(1.5,.2,.3,0))
for (j in 1:(N-1)){
plot(1:T,X[j,],type="l",lwd=1,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5, ylim=c(min(X),max(X)),xaxt='n')
text(T/2,0.8*max(X[j,]),cex=2,bquote(X[.(j)]))
abline(h=0)
}
j=N
plot(1:T,X[j,],type="l",lwd=1,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(min(X),max(X)))
text(T/2,0.8*max(X[j,]),cex=2,bquote(X[.(j)]))
abline(h=0)

dev.off()
##############################
# Underlying time-varying Covariance
##############################

S=array(0,c(N,N,T))
for (t in 1:T){
S[,,t]=L[,,t]%*%diag(G.f[t,])%*%t(L[,,t])
}

#################
# Common components: MC
##############################
M=100
X=array(0,c(T,N,M))


for (m in 1:M){
F=mvrnorm(n = T, mu=rep(0,r), Sigma=diag(r), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
for (t in 1:T){
X[t,,m]=L[,,t]%*%sqrt(diag(G.f[t,]))%*%matrix(F[t,],r,1)}
}

########################################
# Estimated time-varying Covariance: MC
########################################

h=.055
G.h=array(0,c(N,N,TT,M))
for (s in 1:TT){ 	
Ws=diag(dnorm((v[s]- u)/h))/h
for (m in 1:M){ 
G.h[,,s,m]=crossprod(X[,,m],crossprod(Ws,X[,,m]))/sum(diag(Ws))
}}
######################################################################################
#setwd("/Users/GiovanniMotta/Documents/sqrt2-estimation/R codes")
#load("G.h.RData")
#h=.055
#N =dim(G.h)[1]
#TT=dim(G.h)[3]
#M =dim(G.h)[4]

G.M=apply(G.h,FUN=mean,MARGIN=c(1,2,3))
G.L=apply(G.h,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
G.U=apply(G.h,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

mi=min(S)
ma=max(S)
x11()
setpar(mfrow=c(N,N),oma=c(1.2,3.4,1.2,0.5),mar=c(.8,.3,.8,0.8))
for (j in 1:(N-1)){
for (k in 1:N){
G=S[j,k,]
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",cex.axis=1.5)}
else {
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt='n',xaxt="none")}
abline(h=0)
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)], "  vs  " ,tilde(Gamma)[.(j)*','*.(k)])))
}}
j=N
for (k in 1:N){
G=S[j,k,]
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,cex.axis=1.5)}
else {
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt='n',cex.axis=1.5)}
abline(h=0)
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)], "  vs  " ,tilde(Gamma)[.(j)*','*.(k)])))}

#title(paste("N=",N,",  T=",T,",  r=",r),outer=TRUE,cex.main=1.6)

# L are the smooth loadings
# P are the eigenvectors of the smooth matrix S
#dev.print(pdf,file="Gammatilde-est.pdf")

dev.off()

#########################################
# Estimated time-varying eigenvectors: MC
#########################################

P.h=array(0,c(N,r,TT,M))
V.h=array(0,c(TT,r,M))
	
for (s in 1:TT){
for (m in 1:M){
P.h[,,s,m]=		eigen(G.h[,,s,m],symmetric=TRUE)$vectors[,1:r]
V.h[s,,m]=eigen(G.h[,,s,m],symmetric=TRUE)$values[1:r]
}}


VhM=apply(V.h,FUN=mean,MARGIN=c(1,2))
VhL=apply(V.h,MARGIN=c(1,2),FUN=quantile,probs=0.025)
VhU=apply(V.h,MARGIN=c(1,2),FUN=quantile,probs=0.975)

#P.h=sqrt(N)*P.h
#x11()
plot(u,G.f[,1],type="l",col=2,lwd=2,lty=1,ylab="",xlab="", cex.axis=1.5,cex.lab=1.5,ylim=c(-0.035,1.03))
text(1.01,-0.04,cex=2,col=2,bquote(lambda[.(1)]))
lines(v, VhM[,1],lwd=2)
text(1.01,.93,cex=2,expression(tilde("\u2113")[1]))
lines(u,G.f[,2],col=2,lty=5,lwd=2)
text(1.01,1.03,cex=2,bquote(lambda[.(2)]),col=2)
lines(v, VhM[,2],lwd=2,lty=5,col=1)
text(1.01,0.04,cex=2, expression(tilde("\u2113")[2]),col=1)

#dev.print(pdf, file="evals-coal.pdf")
dev.off()

#############################################
c=15  # bandwidth for within-spans signing
#############################################
c2=22 # bandwidth for reconnecting the eventually discontinuous eigenvectors
#############################################

##############################
#swapping the eigenvalues
##############################

Vhat=array(0,c(TT,r,M))
for (m in 1:M){
Vhat[1:(TT/2),,m]      =V.h[1:(TT/2),,m]		# first half: identified
Vhat[(TT/2+1):TT,1,m]  =V.h[(TT/2+1):TT,2,m]	# second half: not-matched 
Vhat[(TT/2+1):(TT),2,m]=V.h[(TT/2+1):TT,1,m]	# we swap
}

VhatM=apply(Vhat,FUN=mean,MARGIN=c(1,2))
VhatL=apply(Vhat,MARGIN=c(1,2),FUN=quantile,probs=0.025)
VhatU=apply(Vhat,MARGIN=c(1,2),FUN=quantile,probs=0.975)


tt=1:TT
cc=(TT/2 -c+1):(TT/2 +c)
vv=v[-cc]

Vtilde=array(0,c(TT-2*c,r,M))
for (m in 1:M){
Vtilde[(1:(TT/2 -c)),,m]  =V.h[(1:(TT/2 -c)),,m]	# first half: identified
Vtilde[(TT/2-c+1):(TT-2*c),1,m]=V.h[(TT/2 +c+1):TT,2,m]	# second half: not-matched 
Vtilde[(TT/2-c+1):(TT-2*c),2,m]=V.h[(TT/2 +c+1):TT,1,m]	# we swap
}

VM=apply(Vtilde,FUN=mean,MARGIN=c(1,2))
VL=apply(Vtilde,MARGIN=c(1,2),FUN=quantile,probs=0.025)
VU=apply(Vtilde,MARGIN=c(1,2),FUN=quantile,probs=0.975)

V.tilde=array(0,c(TT,r,M))
for (m in 1:M){
for (k in 1:r){
y=Vtilde[,k,m]
fit=smooth.spline(vv, y,df=10)
V.tilde[,k,m]=predict(fit,v)$y
}}

V.M=apply(V.tilde,FUN=mean,MARGIN=c(1,2))
V.L=apply(V.tilde,MARGIN=c(1,2),FUN=quantile,probs=0.025)
V.U=apply(V.tilde,MARGIN=c(1,2),FUN=quantile,probs=0.975)

#font_add_google("Dancing Script", "Dancing")
#font_add_google("Parisienne","Parisienne")
#font_add_google("Sacramento","Sacramento")
#font_add_google("Helvetica","Helvetica")
#library(showtext)
#showtext_auto()

mi=min(G.f) 
ma=max(G.f) 

hh=.04
#x11()
#bottom left top right
setpar(mfrow=c(3,r),oma=c(2,3,2,0.5),mar=c(.8,.3,0,1))
#oma=c(1.2,3,1,0.5),mar=c(.8,.3,0,1)
k=1
plot(u,G.f[,k],type="l",ylim=c(-0.05,1.1),col=2,lty=1,lwd=2, cex.axis=1.5,xlab="",xaxt='n',yaxt='n')
axis(2,pos=-0.05,cex.axis=1.5)
text(0.99,0.07,cex=2,expression(lambda[1]),col=2)
text(0.5,0.85*ma,cex=2,expression(tilde("\u2113")[1]))
lines(v,VhM[,k],lty=1,lwd=1.5)
lines(v,VhU[,k],lty=2)
lines(v,VhL[,k],lty=2)
k=2
plot(u,G.f[,k],type="l",col=2,lty=1,lwd=2,cex.axis=1.5,xlab="",xaxt='n',yaxt='n')
text(0.5,0.85*ma,cex=2, bquote(tilde("\u2113")[2]))
lines(v,VhM[,k],lty=1,lwd=1.5)
lines(v,VhU[,k],lty=2)
lines(v,VhL[,k],lty=2)
text(1+0.011,G.f[T,k]-0.02,cex=2,bquote(lambda[.(k)]),col=2)
for (k in 1:r){
if (k==1){
plot(u,G.f[,k],type="l",col=2,lty=1,lwd=2,cex.axis=1.5,xlab="", xaxt='n', yaxt='n')
axis(2,pos=-0.05,cex.axis=1.5)
}
else {
plot(u,G.f[,k],type="l",col=2,lty=1,lwd=2,cex.axis=1.5,xlab="",xaxt='n',yaxt='n')}
lines(vv,VM[,k],lty=1,lwd=1.5)
lines(vv,VU[,k],lty=2)
lines(vv,VL[,k],lty=2)
rect(v[cc[1]],0,v[cc[2*c]],0.85,border=NA,col="white")
rect(v[cc[1]],mi-2*hh,v[cc[2*c]],mi + hh/4,col = rgb(0.5,0.5,0.5,1/4))
lines(u[cc[1:(2*c)]*2],G.f[cc[1:(2*c)]*2,k],col=2,lty=5,lwd=2)
segments(x0=v[cc[1]],   y0=0, x1 = v[cc[1]],   y1=VU[tt[cc[1]-1],k])
segments(x0=v[cc[2*c]], y0=0, x1 = v[cc[2*c]], y1=VU[tt[cc[1]],k])
text(0.5,0.85*ma,cex=2,bquote(hat(lambda)[.(k)]))}
for (k in 1:r){
if (k==1){
plot(u,G.f[,k],type="l",col=2,lty=1,lwd=2,cex.axis=1.5,yaxt='n', xaxt='n', xlab="")
axis(2,pos=-0.05,cex.axis=1.5)}
else {
plot(u,G.f[,k],type="l",col=2,lty=1,lwd=2,cex.axis=1.5,xlab="", yaxt='n',xaxt='n' )}
axis(1,pos=-0.05,cex.axis=1.5)
lines(v,V.M[,k],lty=1,lwd=1.5)
lines(v,V.U[,k],lty=2)
lines(v,V.L[,k],lty=2)
text(0.5,0.9*ma,cex=2,bquote(tilde(lambda)[.(k)]))}

#dev.print(pdf,file="evals-coal-est2.pdf")
dev.off()

###########
#matching the eigenvectors
P.bar=array(0,c(N,r,TT,M))
P.tilde=array(0,c(N,r,TT,M))
###########

for (m in 1:M){
P=P.h[,,,m]
Q1=P[,,(1:(TT/2))]			# first half: identified up-to-sign
Q2=array(0,c(N,r,TT/2))		 
Q2[,1,]=P[,2,((TT/2 + 1):TT)] 	# second half: not-matched and identified up-to-sign
Q2[,2,]=P[,1,((TT/2 + 1):TT)]		# we swap

R1=Q1
R2=Q2

#############################################
# within-spans signing
#############################################

for (k in 1:r){
Q=matrix(Q1[,k,1:(TT/2 - c)],N,TT/2 -c) # "sign-adjusted" eigenvector: N x TT/2
for (t in 1:(TT/2 - c -1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
b=matrix(Q[,t+1]+Q[,t],N,1)
if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}}
R1[,k,1:(TT/2 - c)]=Q
}

for (k in 1:r){
Q=matrix(Q2[,k,((c+1):(TT/2))],N,TT/2 -c) # "signed-adjusted" eigenvector: N x TT/2
for (t in 1:(TT/2 - c -1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
b=matrix(Q[,t+1]+Q[,t],N,1)
if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}}
R2[,k,((c+1):(TT/2))]=Q
}

R=array(0,c(N,r,TT))	
for (k in 1:r){
    R[,k,1:(TT/2)]=R1[,k,]
R[,k,((TT/2)+1):TT]=R2[,k,]
}

P.bar[,,,m]=R

for (k in 1:r){	# picking the right sign, just for simulation purposes
a=matrix(R[,k,1:(TT/2)],N,TT)
b=matrix(LL[,k,1:(TT/2)],N,TT)
if (norm(a-b,"F")> norm(a+b,"F")){P.bar[,k,1:(TT/2),m]=-R1[,k,]}}

for (k in 1:r){	# picking the right sign, just for simulation purposes
a =matrix(R[,k,((TT/2)+1):TT],N,TT)
b=matrix(LL[,k,((TT/2)+1):TT],N,TT)
if (norm(a-b,"F")> norm(a+b,"F")){P.bar[,k,((TT/2)+1):TT,m]=-R2[,k,]}}
#############################################
# anchoring the discontinuous eigenvectors
#############################################
Q=R
for (k in 1:r){
a=matrix(R1[,k,TT/2-c2+1],N,1)
b=matrix(R2[,k,c2],N,1)
if (norm(a-b,"F")> norm(a+b,"F")){Q[,k,(TT/2+1):TT]=-R2[,k,]}
}

P.tilde[,,,m]=Q

for (k in 1:r){	# picking the right sign, just for simulation purposes
a=matrix(Q[,k,],N,TT)
b=matrix(LL[,k,],N,TT)
if (norm(a-b,"F")> norm(a+b,"F")){P.tilde[,k,,m]=-Q[,k,]}
}} #MC loop

Pbar=P.bar
for (m in 1:M){
Pbar[,1,(TT/2 +1):TT,m]=P.bar[,2,(TT/2 +1):TT,m]
Pbar[,2,(TT/2 +1):TT,m]=P.bar[,1,(TT/2 +1):TT,m]
}

PhM=apply(Pbar,MARGIN=c(1,2,3),FUN=mean)
PhL=apply(Pbar,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
PhU=apply(Pbar,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

P.barM=apply(P.bar,MARGIN=c(1,2,3),FUN=mean)
P.barL=apply(P.bar,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
P.barU=apply(P.bar,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

P.tildeM=apply(P.tilde,MARGIN=c(1,2,3),FUN=mean)
P.tildeL=apply(P.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
P.tildeU=apply(P.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

mi=min(L) 
ma=max(L) 

#############
mix=1.5*mi
max=1.5*ma
#############
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2,yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)
axis(2,pos=-0.05,cex.axis=1.6)
}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,-1.1*ma,cex=2,col=1,bquote(paste(Pi[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xaxt='n',yaxt='n',xlab="",ylab="",lwd=2,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt='n',cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,-1.1*ma,cex=2,col=1,bquote(paste(Pi[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}
#dev.print(pdf,file="loads-coal-sim.pdf")
dev.off()

#############

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2,yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)
axis(2,pos=-0.05,cex.axis=1.6)
}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,P.h[j,k,,m],lwd=2)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xaxt='n',yaxt='n',xlab="",ylab="",lwd=2,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt='n',cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,P.h[j,k,,m],lwd=2)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(P)[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}

#dev.print(pdf,file="ignore0.pdf")
dev.off()
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2,yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)
axis(2,pos=-0.05,cex.axis=1.6)
}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,PhM[j,k,],lwd=2)
lines(v,PhL[j,k,],lwd=1.5,lty=2)
lines(v,PhU[j,k,],lwd=1.5,lty=2)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xaxt='n',yaxt='n',xlab="",ylab="",lwd=2,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt='n',cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,PhM[j,k,],lwd=2)
lines(v,PhL[j,k,],lwd=1.5,lty=2)
lines(v,PhU[j,k,],lwd=1.5,lty=2)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}

#dev.print(pdf,file="ignore1.pdf")

dev.off()

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2,yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)
axis(2,pos=-0.05,cex.axis=1.6)
}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1,ylim=c(mix,max), xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,P.barM[j,k,],lwd=2)
lines(v,P.barL[j,k,],lwd=1.5,lty=2)
lines(v,P.barU[j,k,],lwd=1.5,lty=2)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,hat(Pi)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xaxt='n',yaxt='n',xlab="",ylab="",lwd=2,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt='n',cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
lines(v,P.barM[j,k,],lwd=2)
lines(v,P.barL[j,k,],lwd=1.5,lty=2)
lines(v,P.barU[j,k,],lwd=1.5,lty=2)
abline(h=0)
abline(v=0.5)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,hat(Pi)[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}
#dev.print(pdf,file="ignore2.pdf")
dev.off()

hh=.08

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=1.6)
lines(vv,P.barM[j,k,tt[-cc]],lwd=2)
lines(vv,P.barL[j,k,tt[-cc]],lwd=1.5,lty=2)
lines(vv,P.barU[j,k,tt[-cc]],lwd=1.5,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)
}
else{
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt="none",xaxt="none",cex.axis=1.6)
lines(vv,P.barM[j,k,tt[-cc]],lwd=2)
lines(vv,P.barL[j,k,tt[-cc]],lwd=1.5,lty=2)
lines(vv,P.barU[j,k,tt[-cc]],lwd=1.5,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
text(0.5,0.7*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,cex.axis=1.6)
lines(vv,P.barM[j,k,tt[-cc]],lwd=2)
lines(vv,P.barL[j,k,tt[-cc]],lwd=1.5,lty=2)
lines(vv,P.barU[j,k,tt[-cc]],lwd=1.5,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt="none",cex.axis=1.6)
lines(vv,P.barM[j,k,tt[-cc]],lwd=2)
lines(vv,P.barL[j,k,tt[-cc]],lwd=1.5,lty=2)
lines(vv,P.barU[j,k,tt[-cc]],lwd=1.5,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
abline(v=0.5)}
text(0.5,0.7*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}
dev.off()

cc2=(TT/2 -c2+1):(TT/2 +c2)
vv2=v[-cc2]

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",ylab="",lwd=2, xaxt="none",yaxt="none",cex.axis=1.6)
if (k==1){
axis(2,pos=-0.05,cex.axis=1.6)}
lines(vv2,P.tildeM[j,k,tt[-cc2]],lwd=1.5)
lines(vv2,P.tildeL[j,k,tt[-cc2]],lty=2)
lines(vv2,P.tildeU[j,k,tt[-cc2]],lty=2)
rect(v[cc2[1]],mix-0.08,v[cc2[2*c2]],mix+0.1,col = "lightgrey")
rect(v[cc2[1]],mi+1.3*hh,v[cc2[2*c2]],0.95*ma,border=NA,col="white")
segments(x0=v[cc2[1]],   y0=mix-0.08, x1 = v[cc2[1]],   y1=P.tildeU[j,k,tt[cc2[1]-1]])
segments(x0=v[cc2[2*c2]],y0=mix-0.08, x1 = v[cc2[2*c2]], y1=P.tildeU[j,k,tt[cc2[2*c2]]])
lines(u[cc2[1:(2*c2)]*2],L[j,k,cc2[1:(2*c2)]*2],col=2,lty=1,lwd=2)
points(v[cc2[1]],P.tildeM[j,k,tt[cc2[1]]],pch=16)
points(v[cc2[2*c2]],P.tildeM[j,k,tt[cc2[2*c2]]],pch=16)
abline(h=0)
text(0.2,-1.1*ma,cex=2,
bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(Pi)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2,lty=1, ylim=c(mix,max),xlab="",xaxt='n',yaxt='n',ylab="",lwd=2,cex.axis=1.6)
if (k==1){
axis(2,pos=-0.05,cex.axis=1.6)}
lines(vv2,P.tildeM[j,k,tt[-cc2]],lwd=1.5)
lines(vv2,P.tildeL[j,k,tt[-cc2]],lwd=1,lty=2)
lines(vv2,P.tildeU[j,k,tt[-cc2]],lwd=1,lty=2)
rect(v[cc2[1]],mix-0.08,v[cc2[2*c2]],mix+0.1,col = "lightgrey")
rect(v[cc2[1]],mi+1.3*hh,v[cc2[2*c2]],0.95*ma,border=NA,col="white")
segments(x0=v[cc2[1]],   y0=mix-0.08, x1 = v[cc2[1]],   y1=P.tildeU[j,k,tt[cc2[1]-1]])
segments(x0=v[cc2[2*c2]],y0=mix-0.08, x1 = v[cc2[2*c2]], y1=P.tildeU[j,k,tt[cc2[2*c2]]])
lines(u[cc2[1:(2*c2)]*2],L[j,k,cc2[1:(2*c2)]*2],col=2,lty=1,lwd=2)
points(v[cc2[1]],P.tildeM[j,k,tt[cc2[1]]],pch=16)
points(v[cc2[2*c2]],P.tildeM[j,k,tt[cc2[2*c2]]],pch=16)
abline(h=0)
text(0.2,-1.1*ma,cex=2,bquote(
paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(Pi)[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}

#dev.print(pdf,file="Phat-est.pdf")
dev.off()
###################################

Q.tilde=P.tilde

for (m in 1:M){
for (j in 1:N){
for (k in 1:r){
y=P.tilde[j,k,tt[-cc],m]
fit=smooth.spline(vv, y,df=10)
Q.tilde[j,k,,m]=predict(fit,v)$y
}}}

Q.tildeM=apply(Q.tilde,MARGIN=c(1,2,3),FUN=mean)
Q.tildeL=apply(Q.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
Q.tildeU=apply(Q.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

mi=min(L) 
ma=max(L) 
hh=.08

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2, ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.6)
if (k==1){
axis(2,cex.axis=1.6,pos=-0.05)}
lines(v,Q.tildeM[j,k,],lwd=1.5)
lines(v,Q.tildeL[j,k,],lwd=1,lty=2)
lines(v,Q.tildeU[j,k,],lwd=1,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(Pi)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2,ylim=c(mix,max),xlab="",ylab="", yaxt='n',xaxt='n',lwd=2,cex.axis=1.5)
if (k==1){
axis(2,cex.axis=1.6,pos=-0.05)
}	
lines(v,Q.tildeM[j,k,],lwd=1.5)
lines(v,Q.tildeL[j,k,],lwd=1,lty=2)
lines(v,Q.tildeU[j,k,],lwd=1,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
text(0.2,-1.1*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(Pi)[.(j)*','*.(k)])))
axis(1,pos=mi-0.5,cex.axis=1.6)
}
dev.off()
#dev.print(pdf,file="Ptilde-est.pdf")
############################################################

R.tilde=Q.tilde

for (m in 1:M){
Q=Q.tilde[,,,m]
for (h in 1:3){
for (tt in 1:TT){
QQ=t(Q[,,tt])%*%Q[,,tt]
if (rankMatrix(QQ)==r){
QQ=sqrtm(QQ)
QQ=solve(QQ)
R.tilde[,,tt,m]=Q.tilde[,,tt,m]%*%QQ}}
}}

QQ=array(0,c(r,r,TT))
for (tt in 1:TT){
QQ[,,tt]=t(R.tilde[,,tt,m])%*%R.tilde[,,tt,m]}


x11()
setpar(mfrow=c(r,r),oma=c(1.2,1,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
for (j in 1:r){
plot(v,QQ[j,k,],type="l",ylim=c(-2,2),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,1.3,cex=2,bquote(G[.(j)*','*.(k)]^P))
}}
dev.off()


R.tildeM=apply(R.tilde,MARGIN=c(1,2,3),FUN=mean)
R.tildeL=apply(R.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
R.tildeU=apply(R.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2,lty=5,ylim=c(mix,max),xlab="",ylab="",lwd=2, yaxt="none",xaxt="none",cex.axis=1.5)
if (k==1){
	axis(2,at=c(-1.5,-0.5,0.5,1.5),c(-1.5,-0.5,0.5,1.5),cex.axis=1.5)}
lines(v,R.tildeM[j,k,],lwd=1.5)
lines(v,R.tildeL[j,k,],lwd=1,lty=2)
lines(v,R.tildeU[j,k,],lwd=1,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
text(0.8,0.6*ma,cex=1.5,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(Pi)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
plot(u,(L[j,k,]),type="l",col=2,lty=5,ylim=c(mix,max),xlab="",ylab="", 
yaxt="none",lwd=2, cex.axis=1.5)
if (k==1){
	axis(2,at=c(-1.5,-0.5,0.5,1.5),c(-1.5,-0.5,0.5,1.5),cex.axis=1.5)}	
lines(v,R.tildeM[j,k,],lwd=1.5)
lines(v,R.tildeL[j,k,],lwd=1,lty=2)
lines(v,R.tildeU[j,k,],lwd=1,lty=2)
rect(v[cc[1]],mi-2*hh,v[cc[2*c2]],mi+hh,col = rgb(0.5,0.5,0.5,1/4))
abline(h=0)
text(0.8,0.6*ma,cex=1.5,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,tilde(Pi)[.(j)*','*.(k)])))}

#dev.print(pdf,file="Pitilde-est.pdf")
dev.off()
####################################

G.tilde=array(0,c(N,N,TT,M))
for (m in 1:M){ 
for (s in 1:TT){ 	
G.tilde[,,s,m]=R.tilde[,,s,m]%*%diag(V.tilde[s,,m])%*%t(R.tilde[,,s,m])}}

G.M=apply(G.tilde,FUN=mean,MARGIN=c(1,2,3))
G.L=apply(G.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
G.U=apply(G.tilde,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

mi=min(S)
ma=max(S)

x11()
setpar(mfrow=c(N,N),oma=c(1.2,3.4,1.2,0.5),mar=c(.8,.3,.8,0.8))
for (j in 1:(N-1)){
for (k in 1:N){
G=S[j,k,]
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G,type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=1.5)}
else {
plot(u,G,type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt="none")}
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)]^X, "  vs  " ,G[.(j)*','*.(k)]^X)))
}}
j=N
for (k in 1:N){
G=S[j,k,]
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G,type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,cex.axis=1.5)}
else {
plot(u,G,type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',cex.axis=1.5)}
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)]^X, "  vs  " ,G[.(j)*','*.(k)]^X)))}
#title(paste("N=",N,",  T=",T,",  r=",r),outer=TRUE,cex.main=1.6)

##################################
dev.off()


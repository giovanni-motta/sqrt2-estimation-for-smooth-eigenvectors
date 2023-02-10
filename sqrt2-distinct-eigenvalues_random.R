###########################################
# Giovanni Motta
# Copyright © 2023
# estimating eigenvectors of random matrix-valued functions
# sqrt2-smoothing in the case of distinct-eigenvalues 
############################################

##############################
# random loadings: MC simulations
##############################

rm(list=ls())
library(grid)
library(expm)
library(MASS)
setpar=function(...) {
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(0.5,0.5,0.5),
  cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)}
# outer=oma(bottom, left,top,right)
# inner=mar(bottom, left,top,right)
T=400
N=10
r=2
L=array(0,c(N,r,T))
TT=200
u=c(1:T)/T
v=c(1:TT)/TT
LL=array(0,c(N,r,TT))
fr=runif(N,0.6,0.9)
ph=runif(N,0.1,0.3)

G.f=matrix(0,T,r)
G.f[,1]=0.2+2*(u^2-u^3)
G.f[,2]=u^2-u^3

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
MM=t(L[,,tt])%*%L[,,tt]
rk=c(rk,rankMatrix(MM))
if (rankMatrix(MM)==r){
MM=sqrtm(MM)
MM=solve(MM)
L[,,tt]=L[,,tt]%*%MM}}
}
#L=L*sqrt(N)

#####
M=100
#####
F  =array(0,c(T,r,M))
Y  =array(0,c(T,N,M))
X  =array(0,c(T,N,M))
Z  =array(0,c(T,N,M))

x11()
plot(u,G.f[,1],type="l",lwd=2,ylim=c(0,0.5),ylab="",cex.axis=1.5,cex.lab=1.5)
text(1,G.f[T,1],cex=2,bquote(lambda[1]))
lines(u,G.f[,2],col="grey",lwd=2,ylab="")
text(1,G.f[T,2],cex=2,bquote(lambda[2]))

dev.off()

for (m in 1:M){ 
for (t in 1:T){ 
F[t,,m]=mvrnorm(n=1,mu=rep(0,r),Sigma=diag(G.f[t,]), empirical = FALSE)}}

# Common components
for (m in 1:M){ 
for (t in 1:T){
X[t,,m]=L[,,t]%*%matrix(F[t,,m],r,1)
}}

# Idiosyncratic components
g=.8
G.z=matrix(0,N,N)
for (i in 1:N){ 
for (j in 1:N){ 
G.z[i,j]=g^abs(i-j)}}

for (m in 1:M){ 
Z[,,m]=mvrnorm(n = T, mu=rep(0,N), Sigma=G.z, tol = 1e-6, empirical = FALSE)}

# Observations

for (m in 1:M){ 
Y[,,m]=X[,,m]+Z[,,m]}

# Covariance matrix

G.x=array(0,c(N,N,T))
G.y=array(0,c(N,N,T))

for (t in 1:T){
G.x[,,t]=L[,,t]%*%diag(G.f[t,])%*%t(L[,,t])
G.y[,,t]=G.x[,,t]+G.z}

# Smooth covariance matrix  

h=.055
W=array(0,c(T,T,TT))
for (s in 1:TT){
W[,,s]=diag(dnorm((v[s]- u)/h))/h
}

G.h=array(0,c(N,N,TT,M))
for (m in 1:M){ 
for (s in 1:TT){ 	
G.h[,,s,m]=crossprod(sqrt(W[,,s])%*%X[,,m],sqrt(W[,,s])%*%X[,,m])/sum(diag(W[,,s]))
#t(X[,,m])%*%W[,,s]%*%X[,,m]/sum(diag(W[,,s]))
}}
NN=3
mi=min(G.x[1:NN,1:NN,])
ma=max(G.x[1:NN,1:NN,])
# b l t r

x11()
setpar(mfrow=c(NN,NN),oma=c(5.5,6.5,1.2,0.5),mar=c(.8,.4,.8,0.8))
for (j in 1:(NN-1)){
for (k in 1:NN){
G=G.x[j,k,]
if (k==1){
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2, xaxt="none",cex.axis=3.5)
}
else {
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt='n',
xaxt="none")}
abline(h=0)
text(0.3,0.6*ma,cex=3.5,col=2,bquote(Gamma[.(j)*','*.(k)]))
text(0.7,0.6*ma,cex=3.5,col=1,bquote(tilde(Gamma)[.(j)*','*.(k)]))
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1),col=1,lwd=2)
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1)+1.96*apply(G.h[j,k,,],FUN=sd,MARGIN=1), col=1,lty=2)
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1)-1.96*apply(G.h[j,k,,],FUN=sd,MARGIN=1), col=1,lty=2)
}}
j=NN
for (k in 1:NN){
G=G.x[j,k,]
if (k==1){
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt='n',cex.axis=3.5)
axis(1,padj=0.5,pos=mi-0.02,cex.axis=3.5)}
#axis(2,pos=-0.01,cex.axis=2.5)
else {
plot(u,G,type="l",col=2,lty=1,ylim=c(mi,ma),xlab="",ylab="",lwd=2,
yaxt='n',xaxt='n',cex.axis=3.5)
axis(1,c(0,0.25,0.5,0.75,1), c(0,0.25,0.5,0.75,1),padj=0.5,pos=mi-0.02,cex.axis=3.5)
}
abline(h=0)
text(0.3,0.6*ma,cex=3.5,col=2,bquote(Gamma[.(j)*','*.(k)]))
text(0.7,0.6*ma,cex=3.5,col=1,bquote(tilde(Gamma)[.(j)*','*.(k)]))
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1),col=1,lwd=2)
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1)+1.96*apply(G.h[j,k,,],FUN=sd,MARGIN=1),col=1,lty=2)
lines(v,apply(G.h[j,k,,],FUN=mean,MARGIN=1)-1.96*apply(G.h[j,k,,],FUN=sd,MARGIN=1),col=1,lty=2)
}

#dev.print(pdf,file="gamma-distinct-new.pdf")
dev.off()

##############################
  V=matrix(0,T,r)
P.h=array(0,c(N,r,TT,M))
V.h=array(0,c(TT,r,M))

for (t in 1:T){
	V[t,]=eigen(G.x[,,t],symmetric=TRUE)$values[1:r]}
	
for (s in 1:TT){
for (m in 1:M){
P.h[,,s,m]=		eigen(G.h[,,s,m],symmetric=TRUE)$vectors[,1:r]
V.h[s,,m]=eigen(G.h[,,s,m],symmetric=TRUE)$values[1:r]
}}

##############################
# eigenvalues
##############################

mi=min(V) 
ma=max(V) 
# b l t r


setpar(oma=c(2.5,2.5,1.5,0.5), mar=c(2.5,.5,1.5,0.5))
plot(u,V[,1],type="l",col=1,ylim=c(mi,ma), ylab="",lwd=2,xlab="",cex.axis=2.5,cex.lab=2.5,axes=FALSE)
axis(1,padj=0.3,pos=-0.02,cex.axis=2.5)
axis(2,pos=-0.01,cex.axis=2.5)
text(1.02,V[T,1],cex=3,expression("\U2113"[1]))
Vh=apply(V.h[,1,],FUN=mean,MARGIN=1)
lines(v,Vh,col=1,lwd=2,lty=2)
text(1.01,Vh[TT],cex=3, expression(tilde("\U2113")[1]))
lines(u,V[,2],col="darkgrey",lwd=2)
text(1.01,0.0008+V[T,2],cex=3,expression("\U2113"[2]))
Vh=apply(V.h[,2,],FUN=mean,MARGIN=1)
lines(v,Vh,col="darkgrey",lty=2,lwd=2)
text(1.01,Vh[TT],cex=3, expression(tilde("\U2113")[2]))
#title(paste("N=",N,",  T=",T,",  r=",r,",  h=",h),outer=TRUE,cex.main=2.5)

# dev.print(pdf,file="evals-distinct.pdf")
dev.off()


mi=min(V)/N
ma=max(V)/N

x11()
setpar(mfrow=c(1,r),oma=c(1.2,3.4,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
plot(u,V[,k]/N,type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none")
Vh=V.h[,k,]/N
lines(v,apply(Vh,FUN=mean,MARGIN=1),col=1)
lines(v,apply(Vh,FUN=mean,MARGIN=1)+1.96*apply(Vh,FUN=sd,MARGIN=1),col=1,lty=2)
lines(v,apply(Vh,FUN=mean,MARGIN=1)-1.96*apply(Vh,FUN=sd,MARGIN=1),col=1,lty=2)
text(1.02,V[T,k]/N,cex=1,bquote(V[.(k)]))}
abline(h=0)
dev.off()

R.h=P.h
for (k in 1:r){
for (m in 1:M){	
Q=matrix(P.h[,k,,m],N,TT)		# N x 1 "adjusted (the right sign)" eigenvector 
for (t in 1:(TT-1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
b=matrix(Q[,t+1]+Q[,t],N,1)
if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}}
R.h[,k,,m]=Q
}}

L.h=R.h
for (i in 1:NN){
for (k in 1:r){
for (m in 1:M){	
a=matrix(R.h[i,k,,m],TT,1)
b=matrix(  LL[i,k,],TT,1)
if (norm(a-b,"F")> norm(a+b,"F"))
{L.h[i,k,,m] = - R.h[i,k,,m]}		# TT x 1 "adjusted (the right sign)" eigenvector for the MC replications 
}}}

mm=M/2
mi=min(L) 
ma=max(L) 
# b l t r
x11()

setpar(mfrow=c(r,NN),oma=c(2,3,1.3,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:(r-1)){
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=2)
}
else{
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt="none",cex.axis=2)}
abline(h=0)
text(0.5,0.9*ma,cex=2,col=2,bquote(P[.(j)*','*.(k)]))}}
k=r
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",xaxt='n', ylab="",lwd=3,cex.axis=2)
abline(h=0)}
else {
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n', ,xaxt='n',cex.axis=1.6)}
axis(1,pos=mi-0.1,cex.axis=2)
abline(h=0)
text(0.5,0.9*ma,cex=2,col=2,bquote(P[.(j)*','*.(k)]))}

#title(paste("N=",N,",  T=",T,",  r=",r,",  h=",h),outer=TRUE,cex.main=1.6)
#dev.print(pdf,file="loads-distinct-sim.pdf")
dev.off()

mm=M/2
mi=min(L) 
ma=max(L) 
# b l t r
x11()
setpar(mfrow=c(r,NN),oma=c(5,4.5,1.3,0.5),mar=c(.8,.5,.8,1))
for (k in 1:(r-1)){
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,lty=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=3)
}
else{
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3, lty=2, yaxt='n',xaxt="none",cex.axis=3)}
lines(v,P.h[j,k,,mm],col=1,lwd=2)
text(0.3,0.9*ma,cex=3,col=2,bquote(paste(Pi[.(j)*','*.(k)])))
text(0.7,0.9*ma,cex=3,col=1,bquote(paste(tilde(P)[.(j)*','*.(k)])))
abline(h=0)
}}
k=r
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,lty=2,ylim=c(mi,ma),xlab="",xaxt='n', ylab="",lwd=3,cex.axis=3)
abline(h=0)}
else {
plot(u,L[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lty=2, lwd=3,yaxt='n', ,xaxt='n',cex.axis=3)}
axis(1,c(0,0.25,0.5,0.75,1), c(0,0.25,0.5,0.75,1),pos=mi-0.08,padj=0.4,cex.axis=3)
abline(h=0)
lines(v,P.h[j,k,,mm],col=1,lwd=2)
text(0.3,0.9*ma,cex=3,col=2,bquote(paste(Pi[.(j)*','*.(k)])))
text(0.7,0.9*ma,cex=3,col=1,bquote(paste(tilde(P)[.(j)*','*.(k)])))
}

#title(paste("N=",N,",  T=",T,",  r=",r,",  h=",h),outer=TRUE,cex.main=2)
#dev.print(pdf,file="loads-distinct-new.pdf")
dev.off()

L.M=apply(L.h,FUN=mean,MARGIN=c(1,2,3))
L.lbd=array(0,c(NN,r,TT))
L.ubd=L.lbd
for (i in 1:NN){
for (j in 1:r){
for (tt in 1:TT){
L.lbd[i,j,tt]=quantile(L.h[i,j,tt,],probs=0.025)
L.ubd[i,j,tt]=quantile(L.h[i,j,tt,],probs=0.975)}}}

L.lbd=apply(L.h,FUN=quantile,probs=0.025,MARGIN=c(1,2,3))
L.ubd=apply(L.h,FUN=quantile,probs=0.975,MARGIN=c(1,2,3))

x11()
setpar(mfrow=c(r,NN),oma=c(5,4.5,1.3,0.5),mar=c(.8,.5,.8,1))
for (k in 1:(r-1)){
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,lty=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=3)
}
else{
plot(u,L[j,k,],type="l",col=2,lty=2, ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt="none",cex.axis=3)}
lines(v,L.M[j,k,],col=1,lwd=2)
lines(v,L.lbd[j,k,],col=1,lty=2)
lines(v,L.ubd[j,k,],col=1,lty=2)
text(0.3,0.9*ma,cex=3,col=2,bquote(paste(Pi[.(j)*','*.(k)])))
text(0.7,0.9*ma,cex=3,col=1,bquote(paste(tilde(Pi)[.(j)*','*.(k)])))
abline(h=0)
}}
k=r
for (j in 1:NN){
if (j==1){
plot(u,L[j,k,],type="l",col=2,lty=2,ylim=c(mi,ma),xlab="",xaxt='n', ylab="",lwd=3,cex.axis=3)
abline(h=0)}
else {
plot(u,L[j,k,],type="l",col=2,lty=2,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n', ,xaxt='n',cex.axis=3)}
axis(1,c(0,0.25,0.5,0.75,1), c(0,0.25,0.5,0.75,1),pos=mi-0.08,padj=0.4,cex.axis=3)
abline(h=0)
lines(v,L.M[j,k,],col=1,lwd=2)
lines(v,L.lbd[j,k,],col=1,lty=2)
lines(v,L.ubd[j,k,],col=1,lty=2)
text(0.3,0.9*ma,cex=3,col=2,bquote(paste(Pi[.(j)*','*.(k)])))
text(0.7,0.9*ma,cex=3,col=1,bquote(paste(tilde(Pi)[.(j)*','*.(k)])))
}
#dev.print(pdf,file="loads-distinct-adjusted-new.pdf")
dev.off()




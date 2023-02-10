
#############################################
# Wei Biao Wu and Giovanni Motta
# Registration Function 
# Sample (random) matrix-valued functions
# Copyright © 2023 
#############################################

rm(list=ls())
library(expm)
library(MASS)
library(splines)
setpar=function(...) {
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(0.5,0.5,0.5),
  cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)}
# outer=oma(bottom, left,top,right)
# inner=mar(bottom, left,top,right)
# library(Matrix)

T=3000

u=seq(0,1,1/T)
u=u[-1]


##############################
# eigenvalues
##############################

l3=sin(u)+0.02
l2=.2+cos(2*pi*u)/6
l1=.7-u^3/4

u1=u[which(abs(l2-l3)<0.00034)]
u1 #0.21
u2=u[which(abs(l1-l3)<0.0002)]
u2 #0.65

u1=0.21
u2=0.65

d=0.05

s1=c(rep(1,u2*T),rep(3,(1-u2)*T))
s2=c(rep(2,u1*T),rep(3,(u2-u1)*T),rep(1,(1-u2)*T))
s3=c(rep(3,u1*T),rep(2,(1-u1)*T))

length(s1)
length(s2)
length(s3)


x11()
par(mfrow=c(1,2),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(0,u1-d,u1+d,u2-d,u2+d,1),expression(0,u[1]-delta,u[1]+delta,u[2]-delta,u[2]+delta,1),cex.axis=1.5)
lines(u,l1,lty=1,lwd=3)
lines(u,l2,lty=1,lwd=1)
segments(u1,0,u1,l3[T*u1],col="grey") 
segments(u2,0,u2,l1[T*u2],col="grey") 
segments(u1-d,0,u1-d,l3[T*(u1-d)],col="grey") 
segments(u1+d,0,u1+d,l3[T*(u1+d)],col="grey")
segments(u2-d,0,u2-d,l3[T*(u2-d)],col="grey") 
segments(u2+d,0,u2+d,l3[T*(u2+d)],col="grey")

text(1+0.03,l3[T],expression(lambda[3]),cex=2)
text(0,l3[1],expression(lambda[3]),cex=2)
text(0,l2[1]-0.01,expression(lambda[2]),cex=2)
text(0,l1[1]-0.01,expression(lambda[1]),cex=2)
text(1+0.03,l2[T],expression(lambda[2]),cex=2)
text(1+0.03,l1[T],expression(lambda[1]),cex=2)

plot(u,s1,xlim=c(0,1.03),ylim=c(1,3),type="l",lty=1,lwd=2.5,ylab="",xaxt='n',yaxt='n',cex.lab=1.5)
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=1.5)
lines(u,s2,lty=1,lwd=1)
lines(u,s3,lty=2,lwd=1)
text(1+0.03,s3[T],expression(psi[3]),cex=2)
text(1+0.03,s2[T],expression(psi[2]),cex=2)
text(1+0.03,s1[T],expression(psi[1]),cex=2)
text(-0.01,s3[1]-0.02,expression(psi[3]),cex=2)
text(-0.01,s2[1]-0.02,expression(psi[2]),cex=2)
text(-0.01,s1[1]-0.03,expression(psi[1]),cex=2)

dev.off()


N=3
Lambda=matrix(cbind(l1,l2,l3),T,N)
r=N

##############################
# eigenvectors
##############################
P=array(0,c(N,r,T))

TT=T/2
v=seq(0,1,length.out=TT)
PP=array(0,c(N,r,TT))


fr=rep(0.75,N)
ph=rep(0.2,N)

n=1
for (p in 1:r){ #for all t and all p, sum(L.sim[,p,t]^2) = 1 
  P[n,p,]   =((sin(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
 PP[n,p,]   =((sin(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
  P[n+1,p,] =((cos(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
 PP[n+1,p,] =((cos(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
  P[n+2,p,] =((exp(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
 PP[n+2,p,] =((exp(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
}

for (i in 1:6){
rk=0
for (tt in 1:T){
A=t(P[,,tt])%*%P[,,tt]
rk=c(rk,rankMatrix(A))
if (rankMatrix(A)==r){
A=sqrtm(A)
A=solve(A)
P[,,tt]=P[,,tt]%*%A}}
}


A=array(0,c(r,r,T))
for (tt in 1:T){
A[,,tt]=t(P[,,tt])%*%P[,,tt]}

x11()
setpar(mfrow=c(r,r),oma=c(1.2,3.4,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
for (j in 1:r){
plot(u,A[j,k,],type="l",ylim=c(-2,2),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,0.8*n,cex=2,bquote(G[.(j)*','*.(k)]^P))
}}

NN=N
mi=min(P[1:NN,,])
ma=max(P[1:NN,,])
x11()
setpar(mfrow=c(r,NN),oma=c(1.2,3.4,0.5,0.5),mar=c(1.8,.3,.8,0.8))
for (j in 1:(NN-1)){
	for (k in 1:NN){
if (k==1){
plot(u,P[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",cex.axis=2)}
else {
plot(u,P[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",yaxt="none",cex.axis=2)}
text(0.5,0.8*ma,cex=3,bquote(Pi[.(j)*','*.(k)]))}}
j=NN
	for (k in 1:NN){
if (k==1){
plot(u,P[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,cex.axis=2)}
else {
plot(u,P[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt="none",cex.axis=2)}
text(0.5,0.8*ma,cex=3,bquote(Pi[.(j)*','*.(k)]))}

dev.off()


G=P
for (tt in 1:T){
G[,,tt]=P[,,tt]%*%diag(Lambda[tt,])%*%t(P[,,tt])}

M=200
X=array(0,c(T,N,M))

for (m in 1:M){
e=mvrnorm(n=T, mu=rep(0,r), Sigma=diag(r), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
for (t in 1:T){
X[t,,m]=P[,,t]%*%sqrt(diag(Lambda[t,]))%*%matrix(e[t,],r,1)
}}

########################################
# Estimated time-varying Covariance: MC
########################################

h=.04
G.h=array(0,c(N,N,TT,M))
for (s in 1:TT){ 	
Ws=diag(dnorm((v[s]- u)/h))/h
for (m in 1:M){ 
G.h[,,s,m]=crossprod(X[,,m],crossprod(Ws,X[,,m]))/sum(diag(Ws))
}}


N =dim(G.h)[1]
TT=dim(G.h)[3]
M =dim(G.h)[4]

G.M=apply(G.h,FUN=mean,MARGIN=c(1,2,3))
G.L=apply(G.h,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
G.U=apply(G.h,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

mi=min(G)
ma=max(G)
x11()
setpar(mfrow=c(N,N),oma=c(1.2,3.4,1.2,0.5),mar=c(.8,.3,.8,0.8))
for (j in 1:(N-1)){
for (k in 1:N){
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G[j,k,],type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",cex.axis=1.5)}
else {
plot(u,G[j,k,],type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt="none")}
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)]^X, "  vs  " ,tilde(G)[.(j)*','*.(k)]^X)))
}}
j=N
for (k in 1:N){
Gh=G.M[j,k,]
Gu=G.U[j,k,]
Gl=G.L[j,k,]
if (k==1){
plot(u,G[j,k,],type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,cex.axis=1.5)}
else {
plot(u,G[j,k,],type="l",col="grey",ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',cex.axis=1.5)}
lines(v,Gh,lwd=2)
lines(v,Gl,lwd=1.5,lty=2)
lines(v,Gu,lwd=1.5,lty=2)
text(0.5,0.6*ma,cex=2,bquote(paste(Gamma[.(j)*','*.(k)]^X, "  vs  " ,tilde(G)[.(j)*','*.(k)]^X)))}

#title(paste("N=",N,",  T=",T,",  r=",r),outer=TRUE,cex.main=1.6)

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

##########
# order statistics (ordered eigenvalues)
##########

x11()
par(mfrow=c(1,2),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
axis(1,at=c(0,1),expression(0,1),cex.axis=1.5)
lines(u,l1,lty=1,lwd=3)
lines(u,l2,lty=1,lwd=1)
plot(v,VhM[,3],xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",xlab="u",
lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5)
lines(v,VhM[,2],lty=1,lwd=1)
lines(v,VhM[,1],lty=1,lwd=3)
axis(1,at=c(0,u1-d,u1+d,u2-d,u2+d,1),expression(0,u[1]-delta,u[1]+delta,u[2]-delta,u[2]+delta,1),cex.axis=1.5)
axis(1,at=c(0,1),expression(0,1),cex.axis=1.5)

#####
#for (d1 in seq(0.04,0.08,1/100)){
#for (d2 in seq(0.04,0.08,1/100)){

d1=0.04
d2=0.08

t1=round(TT*u1)+2
t1md=t1-TT*d1
t1pd=t1+TT*d1-1
t2=round(TT*u2)
t2md=t2-TT*d2
t2pd=t2+TT*d2
#####

S=cbind(s1,s2,s3)
S=S[seq(2,T,2),]
PP=P[,,seq(2,T,2)]

Pi.hat=P.h
Lambda.hat=V.h
for (m in 1:M){
for (t in 1:TT){
for (i in 1:r){
Pi.hat[,S[t,i],t,m]  =P.h[,i,t,m]
Lambda.hat[t,S[t,i],m]=V.h[t,i,m]
#Pi.hat[,sigma[t,i,m],t,m]=P.h[,i,t,m]
#Lambda.hat[t,sigma[t,i,m],m]=V.h[t,i,m]
}}}


############
t0=c(1:t1md,(t1pd+1):(t2md),(t2pd+1):TT)
mu=12
Lambda.tilde=Lambda.hat
for (m in 1:M){
for (i in 1:r){
#	fit=smooth.spline(x=v[t0],y=Lambda.hat[t0,i,m],df=8)
	fit=smooth.spline(x=v[t0],y=Lambda.hat[t0,i,m],df=12)
Lambda.tilde[,i,m]=predict(fit,v)$y}}

t0=c(1:t1md,(t1pd+1):(t2md),(t2pd+1):TT)
Lambda.tilde=Lambda.hat
x=v[t0]
for (m in 1:M){
for (i in 1:r){
	y=Lambda.hat[t0,i,m]
	df=data.frame(x,y)
fit=lm(y~ns(x,knots=seq(0.01,0.99,length.out=10)),data=df)
fit=lm( y ~ bs(x,degree=3,knots=c(u1,u2)))
Lambda.tilde[,i,m]=predict(fit,newdata=list(x=v))
}}
############

LhatM=apply(Lambda.hat,FUN=mean,MARGIN=c(1,2))
LhatL=apply(Lambda.hat,MARGIN=c(1,2),FUN=quantile,probs=0.025)
LhatU=apply(Lambda.hat,MARGIN=c(1,2),FUN=quantile,probs=0.975)

LtildeM=apply(Lambda.tilde,FUN=mean,MARGIN=c(1,2))
LtildeL=apply(Lambda.tilde,MARGIN=c(1,2),FUN=quantile,probs=0.025)
LtildeU=apply(Lambda.tilde,MARGIN=c(1,2),FUN=quantile,probs=0.975)

###################
# swapped eigenvalues
###################

x11()
par(mfrow=c(1,2),mar=c(5,3,2,3),oma=c(0,1,1,0)) # set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",lty=2,ylab="",xlab="",xaxt='n',
cex.lab=1.5,cex.axis=3)
axis(1,at=c(0,u1,u2,1),padj=0.4,expression(0,u[1],u[2],1),cex.axis=3)
lines(u,l1,lty=1,lwd=3)
lines(u,l2,lty=1,lwd=2)
text(1+0.03,l3[T],expression(lambda[3]),cex=3)
text(0,l3[1],expression(lambda[3]),cex=3)
text(0,l2[1]-0.01,expression(lambda[2]),cex=3)
text(0,l1[1]-0.01,expression(lambda[1]),cex=3)
text(1+0.03,l2[T],expression(lambda[2]),cex=3)
text(1+0.03,l1[T],expression(lambda[1]),cex=3)

plot(v, LhatM[,3],xlim=c(0,1.03),ylim=c(0.032,max(l3)),
type="l",lty=2,xlab="",ylab="",xaxt='n',cex.lab=3,cex.axis=3)
axis(1,padj=0.4,tick = TRUE, at=c(0,1),expression(0,1),cex.axis=3)
axis(1,padj=0.4,tick = FALSE, at=c(u1-d1*1.5,u1+d1*1.5,u2-d2,u2+d2), 
expression(u[1]-delta[1],u[1]+delta[1],u[2]-delta[2],u[2]+delta[2]),cex.axis=2)
lines(v,LhatM[,2],lwd=2,col=1)
lines(v,LhatM[,1],lwd=3)
lines(v[t1md:t1pd],LhatM[t1md:t1pd,3],col="white",lwd=5)
lines(v[t1md:t1pd],LhatM[t1md:t1pd,2],col="white",lwd=5)
lines(v[t2md:t2pd],LhatM[t2md:t2pd,3],col="white",lwd=5)
lines(v[t2md:t2pd],LhatM[t2md:t2pd,1],col="white",lwd=5)
text(1+0.03,l3[T],expression(tilde(lambda)[3]),cex=3)
text(0,l3[1]+0.03,expression(tilde(lambda)[3]),cex=3)
text(0,l2[1]-0.01,expression(tilde(lambda)[2]),cex=3)
text(0,l1[1]-0.01,expression(tilde(lambda)[1]),cex=3)
text(1+0.03,l2[T]-0.01,expression(tilde(lambda)[2]),cex=3)
text(1+0.04,l1[T]+0.02,expression(tilde(lambda)[1]),cex=3)
segments(v[t1md],0,v[t1md], LhatM[t1md,2],col="grey") 
segments(v[t1pd],0,v[t1pd], LhatM[t1pd,3],col="grey") 
segments(v[t2md],0,v[t2md], LhatM[t2md,1],col="grey") 
segments(v[t2pd],0,v[t2pd], LhatM[t2pd,3],col="grey") 
text(u1,0.8,expression(paste(delta[1],"=")),cex=3)
text(u1+0.12,0.804,d1,cex=3)
text(u2,0.8,expression(paste(delta[2],"=")),cex=3)
text(u2+0.12,0.804,d2,cex=3)


#title(main=paste("h=",h),cex=1.5,outer=F)

dev.print(pdf,"Lambda-bar.pdf")

###################
# interpolated eigenvalues
###################

L=cbind(l1,l2,l3)
x11()
par(mfrow=c(1,3),mar=c(5,3,1,0.5),oma=c(0,1,1,0)) # set the margins 
plot(v, LtildeM[,1],lwd=2,xlim=c(0,1.03),ylim=c(0.032,max(l3)),
type="l",lty=1,xlab="",ylab="",cex.lab=2,axes=FALSE)
axis(1,at=c(0,0.25,0.5,0.75,1),c(0,0.25,0.5,0.75,1),padj=0.5,cex.axis=4)
axis(2,at=c(0,0.2,0.4,0.6,0.8),c(0,0.2,0.4,0.6,0.8),padj=0.25,cex.axis=4)
lines(u,L[,1],col=2,lwd=2)
lines(v,LtildeU[,1],col=1,lty=2)
lines(v,LtildeL[,1],col=1,lty=2)
text(0.4,0.85,bquote(lambda[.(1)]),cex=4.5,col=2)
text(0.6,0.85,bquote(tilde(lambda)[.(1)]),cex=4.5,col=1)
for (j in 2:r){
plot(v, LtildeM[,j],lwd=2,xlim=c(0,1.03),ylim=c(0.032,max(l3)),
type="l",lty=1,xlab="",ylab="",yaxt='n',cex.lab=2,axes=FALSE)
axis(1,at=c(0,0.25,0.5,0.75,1),c(0,0.25,0.5,0.75,1),padj=0.5,cex.axis=4)
axis(2,at=c(0,0.2,0.4,0.6,0.8),c(0,0.2,0.4,0.6,0.8),padj=0.25,cex.axis=4)
lines(u,L[,j],col=2,lwd=2)
lines(v,LtildeU[,j],col=1,lty=2)
lines(v,LtildeL[,j],col=1,lty=2)
text(0.4,0.85,bquote(lambda[.(j)]),cex=4.5,col=2)
text(0.6,0.85,bquote(tilde(lambda)[.(j)]),cex=4.5,col=1)
}

#}}

dev.print(pdf,"Lambda-tilde.pdf")

dev.off()

######
# signing
######
Pi.bar= Pi.hat
for (m in 1:M){
for (j in 1:r){
	for (t in 2:TT){
a= Pi.bar[,j,t,m]-Pi.bar[,j,t-1,m]
#b= Pi.bar[,j,t,m]+Pi.bar[,j,t-1,m]
if (sum(a^2)>2){Pi.bar[,j,t,m]=-Pi.bar[,j,t,m]}}}}

######
# anchoring
######
P.bar=Pi.bar
for (m in 1:M){
for (j in 1:r){
	a=P.bar[,j,t1pd,m]-P.bar[,j,t1md,m]
if (sum(a^2)>2){P.bar[,j,t1pd:t2md,m]=-P.bar[,j,t1pd:t2md,m]}
}}
for (m in 1:M){
for (j in 1:r){
	a=P.bar[,j,t2pd,m]-P.bar[,j,t2md,m]
if (sum(a^2)>2){P.bar[,j,t2pd:TT,m]=-P.bar[,j,t2pd:TT,m]}
}}

############
t0=c(1:t1md,(t1pd+1):(t2md),(t2pd+1):TT)
P.tilde=P.bar
for (m in 1:M){
for (i in 1:r){
for (j in 1:r){
fit=smooth.spline(x=v[t0],y=P.bar[i,j,t0,m],df=10)
P.tilde[i,j,,m]=predict(fit,v)$y}}}


t0=c(1:t1md,(t1pd+1):(t2md),(t2pd+1):TT)
P.tilde=P.bar
x=v[t0]
for (m in 1:M){
for (i in 1:r){
for (j in 1:r){
	y=P.bar[i,j,t0,m]
	df=data.frame(x,y)
#fit=lm(y~ns(x,knots=seq(0.01,0.99,length.out=10)),data=df)
#fit=lm( y ~ bs(x,degree=3,knots=c(0,u1,0.43,u2,1)))
fit=lm( y ~ bs(x,degree=3,knots=c(u1,u2)))
P.tilde[i,j,,m]=predict(fit,newdata=list(x=v))}}}
############

P.tilde.plot=P.tilde
for (m in 1:M){
for (s in 1:TT){
for (j in 1:r){
a=PP[,j,s]-P.tilde[,j,s,m]
b=PP[,j,s]+P.tilde[,j,s,m]
#if (norm(a,"F")>norm(b,"F")){P.tilde[,j,t0,m]=-P.tilde[,j,t0,m]}
if (sum(a^2)>sum(b^2)){P.tilde.plot[,j,s,m]=-P.tilde[,j,s,m]}
}}}

PtildeM=apply(P.tilde.plot,FUN=mean,MARGIN=c(1,2,3))
PtildeL=apply(P.tilde.plot,MARGIN=c(1,2,3),FUN=quantile,probs=0.025)
PtildeU=apply(P.tilde.plot,MARGIN=c(1,2,3),FUN=quantile,probs=0.975)

NN=N
mi=min(PP[1:NN,,])
ma=max(PP[1:NN,,])
x11()
setpar(mfrow=c(r,NN),oma=c(1.2,3.4,0.5,0.5),mar=c(1.8,.3,.8,0.8))
for (j in 1:(NN-1)){
	for (k in 1:NN){
if (k==1){
plot(v,PP[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",cex.axis=2)
lines(v,PtildeM[j,k,],lwd=2,lty=1)
lines(v,PtildeL[j,k,],lwd=1,lty=2)
lines(v,PtildeU[j,k,],lwd=1,lty=2)
}
else {
plot(v,PP[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",yaxt="none",cex.axis=2)
lines(v,PtildeM[j,k,],lwd=2,lty=1)
lines(v,PtildeL[j,k,],lwd=1,lty=2)
lines(v,PtildeU[j,k,],lwd=1,lty=2)
}
text(0.4,0.7,cex=3,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.6,0.72,cex=3,bquote(tilde(Pi)[.(j)*','*.(k)]))
}}
j=NN
for (k in 1:NN){
if (k==1){
plot(v,PP[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=2,cex.axis=2)
lines(v,PtildeM[j,k,],lwd=2,lty=1)
lines(v,PtildeL[j,k,],lwd=1,lty=2)
lines(v,PtildeU[j,k,],lwd=1,lty=2)
}

else {
plot(v,PP[j,k,],type="l",col=2,ylim=c(mi,ma),xlab="",ylab="",lwd=2,yaxt="none",cex.axis=2)
lines(v,PtildeM[j,k,],lwd=2,lty=1)
lines(v,PtildeL[j,k,],lwd=1,lty=2)
lines(v,PtildeU[j,k,],lwd=1,lty=2)
}
text(0.4,0.7,cex=3,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.6,0.72,cex=3,bquote(tilde(Pi)[.(j)*','*.(k)]))
}

#
dev.print(pdf,"Pi-tilde.pdf")


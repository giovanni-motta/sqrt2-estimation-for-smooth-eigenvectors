
#################################################
# Wei Biao Wu and Giovanni Motta
# Registration Function 
# Deterministic (non-random) matrix-valued functions
# Copyright © 2023 
#################################################

rm(list=ls())
library(expm)
setpar=function(...) {
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(0.5,0.5,0.5),
  cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)}
# outer=oma(bottom, left,top,right)
# inner=mar(bottom, left,top,right)
# library(Matrix)

T=1000

u=seq(0,1,1/T)
u=u[-1]


##############################
# eigenvalues
##############################

l3=sin(u)+0.02
l2=.2+cos(2*pi*u)/6
l1=.7-u^3/4

u1=u[which(abs(l2-l3)<0.001)]
u1
u2=u[which(abs(l1-l3)<0.00055)]
u2

d=0.05

s1=c(rep(1,u2*T),rep(3,(1-u2)*T))
s2=c(rep(2,u1*T),rep(3,(u2-u1)*T),rep(1,(1-u2)*T))
s3=c(rep(3,u1*T),rep(2,(1-u1)*T))

x11()
par(mfrow=c(1,2),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",lty=1,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5,col="darkgrey",lwd=2)
#segments(u2,0,u2,l1[T*u2],col="darkgrey") 
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
lines(u,l1,lty=1,lwd=2)
lines(u,l2,lty=1,lwd=2,col=2)
#segments(u1,0,u1,l3[T*u1],col="grey") 
text(1+0.03,l3[T],expression(lambda[3]),cex=2)
text(0,l3[1],expression(lambda[3]),cex=2)
text(0,l2[1]-0.01,expression(lambda[2]),cex=2)
text(0,l1[1]-0.01,expression(lambda[1]),cex=2)
text(1+0.03,l2[T],expression(lambda[2]),cex=2)
text(1+0.03,l1[T],expression(lambda[1]),cex=2)

plot(u,s1,xlim=c(0,1.03),ylim=c(1,3),type="l",lty=1,lwd=2,ylab="",xaxt='n',yaxt='n',cex.lab=1.5)
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=1.5)
lines(u,s2,lty=1,lwd=2,col=2)
lines(u,s3,lty=1,lwd=2,col="darkgrey") 
rect(u1-0.005,1-0.005,u1+0.005,3+0.005,col="white",border="transparent")
rect(u2-0.005,1-0.005,u2+0.005,3+0.005,col="white",border="transparent")
text(1+0.03,s3[T],expression(psi[3]),cex=2)
text(1+0.03,s2[T],expression(psi[2]),cex=2)
text(1+0.03,s1[T],expression(psi[1]),cex=2)
text(-0.01,s3[1]-0.02,expression(psi[3]),cex=2)
text(-0.01,s2[1]-0.02,expression(psi[2]),cex=2)
text(-0.01,s1[1]-0.03,expression(psi[1]),cex=2)

#dev.print(pdf,file="registration.pdf")
dev.off()

S=cbind(s1,s2,s3)
N=dim(S)[2]
r=N

x11()
layout(matrix(c(1,1,1,2,3,4), 3, 2, byrow = F),)
par(mar=c(5,3,2,3),oma=c(0,1,0,0)) 
# set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l", ylab="",xlab="",xaxt='n',cex.lab=3,cex.axis=3,col=1,lwd=2,lty=3)
#segments(u2,0,u2,l1[T*u2],col="darkgrey") 
axis(1,padj=0.3,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=3)
lines(u,l1,lty=1,lwd=2)
lines(u,l2,lty=2,lwd=2)
#segments(u1,0,u1,l3[T*u1],col="grey") 
text(1+0.03,l3[T],expression(lambda[3]),cex=3)
text(0,l3[1],expression(lambda[3]),cex=3)
text(0,l2[1]-0.01,expression(lambda[2]),cex=3)
text(0,l1[1]-0.01,expression(lambda[1]),cex=3)
text(1+0.03,l2[T],expression(lambda[2]),cex=3)
text(1+0.03,l1[T],expression(lambda[1]),cex=3)

for (j in 1:r){ 
plot(u,S[,j],xlim=c(0,1.03),xlab="",ylim=c(1,3),type="l",lty=1,lwd=2, ylab="",xaxt='n',yaxt='n',cex.lab=3,cex.axis=3)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=3)
text(0.5,2.5,bquote(psi[.(j)]),cex=3)
}
axis(1,at=c(0,u1,u2,1),padj=0.3,expression(0,u[1],u[2],1),cex.axis=3)

#dev.print(pdf,file="registration-new.pdf")
dev.off()

Lambda=matrix(cbind(l1,l2,l3),T,N)
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
  PP[n,p,]  =((sin(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
  P[n+1,p,] =((cos(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
  PP[n+1,p,]=((cos(p*pi*fr[n]*v - ph[n]))*sqrt(2/N)) 
  P[n+2,p,] =((exp(p*pi*fr[n]*u- ph[n]))*sqrt(2/N)) 
  PP[n+2,p,]=((sin(p*pi*fr[n]*v- ph[n]))*sqrt(2/N)) 
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


dev.off()


G=P
for (tt in 1:T){
G[,,tt]=P[,,tt]%*%diag(Lambda[tt,])%*%t(P[,,tt])}

P.hat=P
L.hat=Lambda

for (t in 1:T){
GG=G[,,t]
L.hat[t,] =eigen(GG,symmetric=TRUE)$values
P.hat[,,t]=eigen(GG,symmetric=TRUE)$vectors
}


#######################
NN=N
mi=min(P.hat[1:NN,,])
ma=max(P.hat[1:NN,,])

x11()
setpar(mfrow=c(r,NN),oma=c(3,9,0.5,0.5),mar=c(2,.5,.8,0.8))
for (j in 1:(NN-1)){
	for (k in 1:NN){
if (k==1){
plot(u,P.hat[j,k,],cex.axis=4,type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2, xaxt="none",yaxt="none",axes=F)
axis(2,pos =-0.1, padj=0.15,at=c(-0.8,-0.4,0,0.4,0.8), expression(-0.8,-0.4,0,0.4,0.8),cex.axis=4)}
else {
plot(u,P.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2, xaxt="none",yaxt="none",axes=F)}
abline(h=0,lty=1,col="grey")
abline(v=c(u1,u2),lty=2,col="grey")
lines(u,P[j,k,],col=2,lwd=2,lty=5)
text(0.3,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4.5,bquote(P[.(j)*','*.(k)]))}}
j=NN
for (k in 1:NN){
if (k==1){
plot(u,P.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,
xaxt="none",yaxt="none",axes=FALSE)
axis(1,padj=0.7,at=c(0,u1,0.5,u2,1),expression(0,u[1],0.5,u[2],1),cex.axis=4,pos = -1.1)
axis(2,pos=-0.1, padj=0.15,at=c(-0.8,-0.4,0,0.4,0.8),expression(-0.8,-0.4,0,0.4,0.8),cex.axis=4)
}
else{
plot(u,P.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",yaxt="none",axes=FALSE)
axis(1,padj=0.7,at=c(0,u1,0.5,u2,1),expression(0,u[1],0.5,u[2],1),cex.axis=4,pos = -1.1)}
abline(h=0,lty=1,col="grey")
abline(v=c(u1,u2),lty=2,col="grey")
lines(u,P[j,k,],col=2,lwd=2,lty=5)
text(0.3,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4.5,bquote(P[.(j)*','*.(k)]))}

#dev.print(pdf,file="P.pdf")
dev.off()

G.hat=G
for (tt in 1:T){
G.hat[,,tt]=P.hat[,,tt]%*%diag(L.hat[tt,])%*%t(P.hat[,,tt])}

x11()
par(mfrow=c(r,r),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
for (i in 1:r){
	for (j in 1:r){
plot(u,G[i,j,],xlim=c(0,1.03),ylim=c(min(G.hat),max(G.hat)),
type="l",lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5,lwd=3)
lines(u,G.hat[i,j,],col=2)
}}


#####
d=0.02
t1=T*u1
t1md=T*(u1-d)
t1pd=T*(u1+d)
t2=T*u2
t2md=T*(u2-d)
t2pd=T*(u2+d)
#####

sigma=matrix(0,T,r)
for (t in 1:t1md){
for (i in 1:r){sigma[t,i]=i}}

for (t in t1pd:t2md){
for (i in 1:r){
dP=matrix(0,2,r)	
for (j in 1:r){
a=P.hat[,i,t1pd]-P.hat[,sigma[t1md,j],t1md]
b=P.hat[,i,t1pd]+P.hat[,sigma[t1md,j],t1md]
dP[1,j]=sqrt(sum(a^2))
dP[2,j]=sqrt(sum(b^2))}
dPmin=apply(dP,FUN=min,MARGIN=2)
sigma[t,i]=which(dPmin==min(dPmin))
}}

for (t in t2pd:T){
for (i in 1:r){
dP=matrix(0,2,r)	
for (j in 1:r){
a= P.hat[,i,t2pd]-P.hat[,sigma[t2md,j],t2md]
b=-P.hat[,i,t2pd]-P.hat[,sigma[t2md,j],t2md]
dP[1,j]=sum(a^2)
dP[2,j]=sum(b^2)}
dPmin=apply(dP,FUN=min,MARGIN=2)
sigma[t,i]=which(dPmin==min(dPmin))
}}

for (t in t2pd:T){
for (i in 1:r){
dP=matrix(0,2,r)	
for (j in 1:r){
a= P.hat[,i,t2pd]-P.hat[,sigma[t2md,j],t2md]
b=-P.hat[,i,t2pd]-P.hat[,sigma[t2md,j],t2md]
dP[1,j]=sum(a^2)
dP[2,j]=sum(b^2)}
dPmin=apply(dP,FUN=min,MARGIN=2)
sigma[t,i]=which(dPmin==min(dPmin))
}}

#t0=c((t1md+1):(t1pd-1),(t2md+1):(t2pd-1))
#u0=u[-t0]

x11()
plot(u,sigma[,1],xlim=c(0,1.03),ylim=c(1,3),type="l",lty=1,lwd=2, ylab="",xaxt='n',yaxt='n',cex.lab=1.5)
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=1.5)
lines(u,sigma[,2],col=2,lty=1,lwd=2)
lines(u,sigma[,3],lty=5,lwd=2,col="darkgrey")
rect(u1-d-0.01,1-0.075,u1+d+0.01,3+0.05, col="white",border="transparent")
rect(u2-d-0.01,1-0.075,u2+d+0.01,3+0.05, col="white",border="transparent")
rect(u1-d-0.01,1-0.075,u1+d+0.01,0.95, col="lightgrey")
rect(u2-d-0.01,1-0.075,u2+d+0.01,0.95, col="lightgrey")
text(1+0.03,sigma[T,3],expression(psi[3]),cex=2)
text(1+0.03,sigma[T,2],expression(psi[2]),cex=2)
text(1+0.03,sigma[T,1],expression(psi[1]),cex=2)
text(-0.01,sigma[1,3]-0.02,expression(psi[3]),cex=2)
text(-0.01,sigma[1,2]-0.02,expression(psi[2]),cex=2)
text(-0.01,sigma[1,1]-0.03,expression(psi[1]),cex=2)

#dev.print(pdf,file="reg-non_random.pdf")
dev.off()

Lambda.hat=L.hat
Pi.hat=P.hat
for (t in 1:T){
for (i in 1:r){
Pi.hat[,sigma[t,i],t]=P.hat[,i,t]
Lambda.hat[t,sigma[t,i]]=L.hat[t,i]}}

x11()
par(mfrow=c(1,2),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
plot(u,l3,xlim=c(0,1.03),ylim=c(0.032,max(l3)),type="l",lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(0,u1-d,u1+d,u2-d,u2+d,1),expression(0,u[1]-delta,u[1]+delta,u[2]-delta,u[2]+delta,1),cex.axis=1.5)
axis(1,at=c(0,1),expression(0,1),cex.axis=1.5)
lines(u,l1,lty=1,lwd=3)
lines(u,l2,lty=1,lwd=1)
plot(u, Lambda.hat[,3],col=1,xlim=c(0,1.03),ylim=c(0.032,max(l3)),
type="l",lty=2,ylab="",xaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(0,u1-d,u1+d,u2-d,u2+d,1),expression(0,u[1]-delta,u[1]+delta,u[2]-delta,u[2]+delta,1),cex.axis=1.5)
lines(u, Lambda.hat[,1],lty=1,lwd=3,col=1)
lines(u, Lambda.hat[,2],lty=1,lwd=1,col=1)
lines(u[t1md:t1pd],Lambda.hat[t1md:t1pd,2],lty=1,lwd=4,col="white")
lines(u[t1md:t1pd],Lambda.hat[t1md:t1pd,3],lty=1,lwd=4,col="white")
lines(u[t2md:t2pd],Lambda.hat[t2md:t2pd,3],lty=1,lwd=2,col="white")
lines(u[t2md:t2pd],Lambda.hat[t2md:t2pd,1],lty=1,lwd=4,col="white")

dev.off()

d.plot=2*d
x11()
par(mfrow=c(1,2),mar=c(4,2,2,3),oma=c(0,1,0,0)) # set the margins of the plots
plot(u,Lambda.hat[,3],col="darkgrey",xlim=c(0,1.03), ylim=c(0.032,max(l3)),type="l",lty=1,ylab="",xaxt='n', cex.lab=1.5,cex.axis=1.5,lwd=2)
#axis(1,at=c(u1,u2),expression(u[1],u[2]),cex.axis=1.5)
axis(1,at=c(u1-d,u1+d),expression("",""),cex.axis=1)
axis(1,at=c(u1-d.plot,u1+d.plot),expression(u[2]-delta,u[2]+delta),cex.axis=1.3,lwd.ticks= 0)
axis(1,at=c(u2-d,u2+d),expression("",""),cex.axis=1)
axis(1,at=c(u2-d.plot,u2+d.plot),expression(u[2]-delta,u[2]+delta),cex.axis=1.3,lwd.ticks= 0)
axis(1,at=c(0,1),expression(0,1),cex.axis=1.5)
lines(u,Lambda.hat[,1],lty=1,lwd=2)
lines(u,Lambda.hat[,2],lty=1,lwd=2,col=2)
#segments(u1,0,u1,l3[T*u1],col="grey") 
#segments(u2,0,u2,l1[T*u2],col="grey") 
# segments(u1-d,0,u1-d,l3[T*(u1-d)],col="grey") 
# segments(u1+d,0,u1+d,l3[T*(u1+d)],col="grey")
# segments(u2-d,0,u2-d,l3[T*(u2-d)],col="grey") 
# segments(u2+d,0,u2+d,l3[T*(u2+d)],col="grey")
lines(u[(t1md-1):(t1pd+1)],Lambda.hat[(t1md-1):(t1pd+1),2],lty=1,lwd=4,col="white")
lines(u[(t1md-1):(t1pd+1)],Lambda.hat[(t1md-1):(t1pd+1),3],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),1],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),2],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),3],lty=1,lwd=4,col="white")

abline(v=c(u1-d,u1+d,u2-d,u2+d),lty=2,col="grey")

text(0,l3[1],     expression(lambda[3]^o),cex=2)
text(0,l2[1],expression(lambda[2]^o),cex=2)
text(0,l1[1],expression(lambda[1]^o),cex=2)
text(1+0.03,l3[T],expression(lambda[3]^o),cex=2)
text(1+0.03,l2[T],expression(lambda[2]^o),cex=2)
text(1+0.03,l1[T],expression(lambda[1]^o),cex=2)

plot(u,sigma[,1],xlim=c(0,1.03),ylim=c(1,3),type="l",lty=1,lwd=2,ylab="",xaxt='n',yaxt='n',cex.lab=1.5)
axis(1,at=c(0,u1,u2,1),expression(0,u[1],u[2],1),cex.axis=1.5)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=1.5)
lines(u,sigma[,2],lty=1,lwd=2,col=2)
lines(u,sigma[,3],lty=1,lwd=2,col="darkgrey")
text(1+0.03,sigma[T,3],expression(psi[3]^o),cex=2)
text(1+0.03,sigma[T,2],expression(psi[2]^o),cex=2)
text(1+0.03,sigma[T,1],expression(psi[1]^o),cex=2)
text(-0.01,sigma[1,3],expression(psi[3]^o),cex=2)
text(-0.01,sigma[1,2],expression(psi[2]^o),cex=2)
text(-0.01,sigma[1,1],expression(psi[1]^o),cex=2)
# lines(u[(t1md-1):(t1pd+1)],rep(1,length((t1md-1):(t1pd+1))),lty=1,lwd=4,col="white")
# lines(u[(t1md-1):(t1pd+1)],Lambda.hat[(t1md-1):(t1pd+1),3],lty=1,lwd=4,col="white")
# lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),3],lty=1,lwd=2,col="white")
# lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),1],lty=1,lwd=4,col="white")
# lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),2],lty=1,lwd=2,col="white")
rect(u1-d-0.01,1-0.075,u1+d+0.01,3+0.05, col="white",border="transparent")
rect(u2-d-0.01,1-0.075,u2+d+0.01,3+0.05, col="white",border="transparent")
rect(u1-d-0.01,1-0.08,u1+d+0.01,0.95, col="lightgrey")
rect(u2-d-0.01,1-0.08,u2+d+0.01,0.95, col="lightgrey")

#dev.print(pdf,"reg-non_random.pdf")
dev.off()

S=cbind(s1,s2,s3)
N=dim(S)[2]
r=N

x11()
layout(matrix(c(1,1,1,2,3,4), 3, 2, byrow = F),)
par(mar=c(5,3,2,3),oma=c(0,1,0,0)) 
# set the margins of the plots
plot(u,Lambda.hat[,3],col=1,xlim=c(0,1.03), ylim=c(0.032,max(l3)), type="l",xlab="",ylab="",xaxt='n', cex.lab=3.5,cex.axis=3,lwd=2,lty=3)
axis(1,padj=0.3,at=c(u1-d,u1+d),expression("",""),cex.axis=2.5)
axis(1,padj=0.3,at=c(u1-d.plot*1.1,u1+d.plot*1.1), 
expression(u[1]-delta,u[1]+delta),cex.axis=2.5,lwd.ticks= 0)
axis(1,padj=0.3,at=c(u2-d,u2+d),expression("",""),cex.axis=2.5)
axis(1,padj=0.3,at=c(u2-d.plot*1.1,u2+d.plot*1.1),
expression(u[2]-delta,u[2]+delta),cex.axis=2.5,lwd.ticks= 0)
axis(1,at=c(0,1),padj=0.3,expression(0,1),cex.axis=3)
lines(u,Lambda.hat[,1],lty=1,lwd=2)
lines(u,Lambda.hat[,2],lty=2,lwd=2)
lines(u[(t1md-1):(t1pd+1)],Lambda.hat[(t1md-1):(t1pd+1),2],lty=1,lwd=4,col="white")
lines(u[(t1md-1):(t1pd+1)],Lambda.hat[(t1md-1):(t1pd+1),3],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),1],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),2],lty=1,lwd=4,col="white")
lines(u[(t2md-1):(t2pd+1)],Lambda.hat[(t2md-1):(t2pd+1),3],lty=1,lwd=4,col="white")
abline(v=c(u1-d,u1+d,u2-d,u2+d),lty=2,col="grey")
text(0,l3[1],     expression(lambda[3]^o),cex=3)
text(0,l2[1],expression(lambda[2]^o),cex=3)
text(0,l1[1],expression(lambda[1]^o),cex=3)
text(1+0.03,l3[T],expression(lambda[3]^o),cex=3)
text(1+0.03,l2[T],expression(lambda[2]^o),cex=3)
text(1+0.03,l1[T],expression(lambda[1]^o),cex=3)

for (j in 1:r){ 
plot(u,sigma[,j],xlim=c(0,1.03),ylim=c(1,3),type="l",lty=1,lwd=2,
ylab="",xlab="",xaxt='n',yaxt='n',cex.lab=2.5)
axis(2,at=c(1,2,3),c(1,2,3),cex.axis=3)
text(0.5,2.5,bquote(psi[.(j)]^o),cex=3)
rect(u1-d-0.01,1-0.075,u1+d+0.01,3+0.05, col="white",border="transparent")
rect(u2-d-0.01,1-0.075,u2+d+0.01,3+0.05, col="white",border="transparent")
}
rect(u1-d-0.01,1-0.08,u1+d+0.01,1, col="lightgrey")
rect(u2-d-0.01,1-0.08,u2+d+0.01,1, col="lightgrey")
axis(1,padj=0.3,at=c(u1-d*1.5,u1+d*1.5),expression("",""),cex.axis=2.5)
axis(1,padj=0.3,at=c(u1-d.plot*1.1,u1+d.plot*1.1),
expression(u[1]-delta,u[1]+delta),cex.axis=2.5,lwd.ticks= 0)
axis(1,padj=0.3,at=c(u2-d*1.5,u2+d*1.5),expression("",""),cex.axis=2.5)
axis(1,padj=0.3,at=c(u2-d.plot*1.1,u2+d.plot*1.1),expression(u[2]-delta,u[2]+delta),cex.axis=2.5,lwd.ticks= 0)
axis(1,padj=0.3,at=c(0,1),expression(0,1),cex.axis=3)

#dev.print(pdf,"reg-non_random-new.pdf")
dev.off()

#######################
NN=N
mi=min(P.hat[1:NN,,])
ma=max(P.hat[1:NN,,])

x11()
setpar(mfrow=c(r,NN),oma=c(3,9,0.5,0.5),mar=c(2,.5,.8,0.8))
for (j in 1:(NN-1)){
	for (k in 1:NN){
if (k==1){
plot(u,Pi.hat[j,k,],cex.axis=4,type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2, xaxt="none",yaxt="none",axes=F)
axis(2,pos =-0.1, padj=0.15,at=c(-0.8,-0.4,0,0.4,0.8), expression(-0.8,-0.4,0,0.4,0.8),cex.axis=4)}
else {
plot(u,Pi.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2, xaxt="none",yaxt="none",axes=F)}
abline(h=0,lty=1,col="grey")
abline(v=c(u1,u2),lty=2,col="grey")
lines(u,P[j,k,],col=2,lwd=2,lty=5)
text(0.3,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]^o))}}
j=NN
for (k in 1:NN){
if (k==1){
plot(u,Pi.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,
xaxt="none",yaxt="none",axes=FALSE)
axis(1,padj=0.7,at=c(0,u1,0.5,u2,1),expression(0,u[1],0.5,u[2],1),cex.axis=4,pos = -1.1)
axis(2,pos=-0.1, padj=0.15,at=c(-0.8,-0.4,0,0.4,0.8),expression(-0.8,-0.4,0,0.4,0.8),cex.axis=4)
}
else{
plot(u,Pi.hat[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none",yaxt="none",axes=FALSE)
axis(1,padj=0.7,at=c(0,u1,0.5,u2,1),expression(0,u[1],0.5,u[2],1),cex.axis=4,pos = -1.1)}
abline(h=0,lty=1,col="grey")
abline(v=c(u1,u2),lty=2,col="grey")
lines(u,P[j,k,],col=2,lwd=2,lty=5)
text(0.3,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4.5,bquote(Pi[.(j)*','*.(k)]^o))}

#dev.print(pdf,"Pi.pdf")
dev.off()

##########
#signing
###########

Pi.bar= Pi.hat
t0=c(2:t1md,(t1pd+1):(t2md),(t2pd+1):T)
for (j in 1:r){
	for (t in t0){
a= Pi.bar[,j,t]-Pi.bar[,j,t-1]
b= Pi.bar[,j,t]+Pi.bar[,j,t-1]
if (sum(a^2)>2){Pi.bar[,j,t]=-Pi.bar[,j,t]}}}


##########
#anchoring
###########

P.bar=Pi.bar

 for (j in 1:r){
 if (sign(Pi.bar[1,j,t1md])+sign(Pi.hat[1,j,t1md])==0){P.bar[,j,1:t1md]=-P.bar[,j,1:t1md]}
 if (sign(Pi.bar[1,j,t2md])+sign(Pi.hat[1,j,t2md])==0){P.bar[,j,1:t2md]=-P.bar[,j,1:t2md]}
 }

for (j in 1:r){
	a=P.bar[,j,t1pd]-P.bar[,j,t1md]
if (sum(a^2)>2){P.bar[,j,t1pd:t2md]=-P.bar[,j,t1pd:t2md]}
}

for (j in 1:r){
	a=P.bar[,j,t2pd]-P.bar[,j,t2md]
if (sum(a^2)>2){P.bar[,j,t2pd:T]=-P.bar[,j,t2pd:T]}
}

NN=N
mi=min(P[1:NN,,])
ma=max(P[1:NN,,])

x11()
setpar(mfrow=c(r,NN),oma=c(4,6,0.5,0.5),mar=c(1.8,.3,.8,0.8))
for (j in 1:(NN-1)){
	for (k in 1:NN){
if (k==1){
plot(u,P.bar[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2.5, yaxt="none",xaxt="none",cex.axis=3)
axis(2,pos=-0.07,cex.axis=3)
}
else {
plot(u,P.bar[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2.5,
xaxt="none",yaxt="none",cex.axis=3)}
lines(u,P[j,k,],col=2,lwd=1.5,lty=2)
text(0.3,0.8*ma,cex=4,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4,bquote(Pi[.(j)*','*.(k)]^o))
abline(h=0,lty=1,col="grey")
}}
j=NN
	for (k in 1:NN){
if (k==1){
plot(u,P.bar[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2.5, xaxt="none",yaxt="none",cex.axis=3)
axis(1,padj=0.5,pos=1.2*mi,cex.axis=3)
axis(2,pos=-0.07,cex.axis=3)
}
else {
plot(u,P.bar[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2.5, xaxt="none",yaxt="none",cex.axis=3)
axis(1,padj=0.5,pos=1.2*mi,cex.axis=3)
}
lines(u,P[j,k,],col=2,lwd=1.5,lty=2)
text(0.3,0.8*ma,cex=4,bquote(Pi[.(j)*','*.(k)]),col=2)
text(0.7,0.8*ma,cex=4,bquote(Pi[.(j)*','*.(k)]^o))
abline(h=0,lty=1,col="grey")
}

#dev.print(pdf,"Piobar-new.pdf")
dev.off()


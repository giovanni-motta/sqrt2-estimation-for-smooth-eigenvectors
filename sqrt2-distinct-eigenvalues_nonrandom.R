
###########################################
# Giovanni Motta
# Copyright © 2023
# non-random (deterministic) matrix-valued functions
# sqrt2-smoothing in the case of distinct-eigenvalues 
###########################################
rm(list=ls())
#library(ggplot2)
library(grid)
#library(dplyr)
#library(knitr)
#library(data.table)
#library(pracma)
#library(Matrix)
library(expm)
library(MASS)

#------some graphical parameters------

setpar=function(...) {
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(0.5,0.5,0.5),
  cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)}
# outer=oma(bottom, left,top,right)
# inner=mar(bottom, left,top,right)

N=4
T=400
u=c(1:T)/T
uu=u[-T]
r=2
L=array(0,c(N,r,T))

##############################
# eigenvectors
##############################

fr=runif(N,0.6,0.9)
ph=runif(N,0.1,0.3)

fr=seq(0.6,0.9,length.out=N)
ph=seq(0.1,0.3,length.out=N)

nn=seq(1,N,2)
for (n in nn){
for (p in 1:r){ #for all t and all p, sum(L.sim[,p,t]^2) = 1 
  L[n,p,] =((sin(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
  L[n+1,p,] =((cos(p*pi*fr[n]*u - ph[n]))*sqrt(2/N)) 
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
plot(u,L[j,k,],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,0.8*ma,cex=2,bquote(Lambda[.(j)*','*.(k)]))
}}

dev.off()

M=array(0,c(r,r,T))
for (tt in 1:T){
M[,,tt]=t(L[,,tt])%*%L[,,tt]}

x11()
setpar(mfrow=c(r,r),oma=c(1.2,3.4,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
for (j in 1:r){
#if (j==1){
#plot(uu,L.sim[j,k,26:T],type="l",ylim=c(mi,ma),xlab="",ylab="",lwd=2,xaxt="none")}
#else {
plot(u,M[j,k,],type="l",ylim=c(-2,2),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,0.75,cex=2,bquote(paste(L[.(j)*','*.(k)]^T, L[.(j)*','*.(k)])))
}}

dev.off()

G.f=matrix(0,T,r)
G.f[,1]=0.2+2*(u^2-u^3)
G.f[,2]=u^2-u^3

#library(showtext)
#showtext_auto()
# b l t r
mi=min(G.f) 
ma=max(G.f) 
setpar(oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
plot(u,G.f[,1],type="l",ylim=c(mi,ma),xlab="x",ylab="",lwd=2,cex.axis=1.6)
text(1.02,.2,cex=2, expression("\u2113"[1]))
lines(u,G.f[,2],lwd=2)
text(1.02,0,cex=2,expression("\u2113"[2]))
#dev.print(pdf, file="evals1.pdf")

dev.off()
S=array(0,c(N,N,T))
for (t in 1:T){
S[,,t]=L[,,t]%*%diag(G.f[t,])%*%t(L[,,t])
}

NN=N
mi=min(S)
ma=max(S)

x11()
setpar(mfrow=c(NN,NN),oma=c(1.2,3.4,1.2,0.5),mar=c(.8,.3,.8,0.8))
for (j in 1:(NN-1)){
for (k in 1:NN){
G=S[j,k,]
if (k==1){
plot(u,G,type="l",col=1,ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt="none",yaxt="none", cex.axis=1.5)
abline(h=0)
axis(2,pos=-0.07,cex.axis=1.5)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}
else {
plot(u,G,type="l",col=1,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt="none")
abline(h=0)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}
}}
j=NN
for (k in 1:NN){
G=S[j,k,]
if (k==1){
plot(u,G,type="l",col=1,ylim=c(mi,ma),xlab="",ylab="",lwd=3,xaxt='n',cex.axis=1.5,yaxt='n')
axis(1,pos=mi-0.04,cex.axis=1.5)
axis(2,pos=-0.07,cex.axis=1.5)
abline(h=0)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}
else {
plot(u,G,type="l",col=1,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt='n')
axis(1,pos=mi-0.04,cex.axis=1.5)
abline(h=0)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}}
#dev.print(pdf,file="Ax-dist-talk.pdf")
dev.off()

V=matrix(0,T,r)
P=array(0,c(N,r,T))	
for (t in 1:T){
	V[t,]=eigen(S[,,t],symmetric=TRUE)$values[1:r]
	P[,,t]=eigen(S[,,t],symmetric=TRUE)$vectors[,1:r]}

R=P
for (k in 1:r){
Q=matrix(P[,k,],N,T)		# N x 1 "adjusted (the right sign)" eigenvector 
for (t in 1:(T-1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
b=matrix(Q[,t+1]+Q[,t],N,1)
if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}}
R[,k,]=Q
}

NN=N
mi=min(L) 
ma=max(L) 
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.3,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),col=2,lwd=2, ylim=c(mi,ma),type="l", xlab="",ylab="",xaxt="none",yaxt="none",cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),col=2,lwd=2,type="l",
ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(u,(L[j,k,]),col=2,lwd=2,type="l",xaxt="none",yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(u,(L[j,k,]),col=2,lwd=2,type="l",yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)])))
}
#dev.print(pdf,file="P-talk.pdf")
dev.off()

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,P[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,P[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,R[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,P[j,k,1:(T-1)],type="l",col=1,lwd=1.5,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,P[j,k,1:(T-1)],type="l",col=1,lwd=1.5,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,R[.(j)*','*.(k)])))
}

#dev.print(pdf,file="loads-distinct-nonr-talk.pdf")
dev.off()

Q=abs(P)

uu=u[-T]
NN=N
mi=min(L) 
ma=max(L) 
# b l t r

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(abs(P)[.(j)*','*.(k)], "  vs  " ,abs(R)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=1.5,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=1.5,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(abs(P)[.(j)*','*.(k)], "  vs  " ,abs(R)[.(j)*','*.(k)])))
}

#dev.print(pdf,file="loads-distinct-nonr-abs.pdf")
dev.off()

R=-R

NN=N
mi=min(L) 
ma=max(L) 
# b l t r
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2.5, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2.5, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2.5,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2.5,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(P[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}

#dev.print(pdf,file="loads-distinct-nonr-adj.pdf")

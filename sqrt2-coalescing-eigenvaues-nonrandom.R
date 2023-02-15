
###########################################
# Giovanni Motta
# Copyright © 2023
# non-random (deterministic) matrix-valued functions
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
T=400
u=c(1:T)/T
r=2
L=array(0,c(N,r,T))

TT=400
v=seq(0,1,length.out=TT)
LL=array(0,c(N,r,TT))

##############################
# loadings
##############################

fr=runif(N,0.6,0.9)
ph=runif(N,0.1,0.3)

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

M=array(0,c(r,r,T))
for (tt in 1:T){
M[,,tt]=t(L[,,tt])%*%L[,,tt]}

x11()
setpar(mfrow=c(r,r),oma=c(1.2,3.4,0.5,0.5),mar=c(.8,.3,.8,0.8))
for (k in 1:r){
for (j in 1:r){
plot(u,M[j,k,],type="l",ylim=c(-2,2),xlab="",ylab="",lwd=2,xaxt="none")
text(0.5,0.8*n,cex=2,bquote(G[.(j)*','*.(k)]^P))
}}
dev.off()


 library(showtext)
 library(jsonlite)
 library(curl)
 
 font_add_google("Sacramento","Sacramento")
 showtext_auto()

##############################
# Factors
G.f=matrix(0,T,r)
G.f[,1]= .5*(1+cos(pi*u))
G.f[,2]= .5*(1-cos(pi*u))

x11()
plot(u,G.f[,1],type="l",lwd=2,ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1),xlim=c(0,1.05))
text(1.04,0,cex=2,bquote(lambda[.(1)]))
lines(u,G.f[,2],col=2,lwd=2,lty=5)
text(1.04,1,cex=2,bquote(lambda[.(2)]),col=2)
#dev.print(pdf,file="evals-coal-smooth.pdf")
dev.off()

x11()
plot(u,apply(G.f,FUN=max,MARGIN=1),type="l",lwd=2,ylab="",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1),xlim=c(0,1.05))
 text(1.02,1,cex=2,bquote(l),family="Sacramento")
 text(1.045,0.99,cex=2,bquote(""[.(1)]))
 lines(u,apply(G.f,FUN=min,MARGIN=1),col=2,lwd=2,lty=5)
 text(1.02,0,cex=2,bquote(l),family="Sacramento",col=2)
 text(1.045,-0.01,cex=2,bquote(""[.(2)]),col=2)

#dev.print(pdf,file="evals-coal-unsmooth.pdf")
dev.off()

S=array(0,c(N,N,T))
for (t in 1:T){
S[,,t]=L[,,t]%*%diag(G.f[t,])%*%t(L[,,t])
}

mi=min(S[1:NN,1:NN,])
ma=max(S[1:NN,1:NN,])

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
abline(h=0)
axis(1,pos=mi-0.05,cex.axis=1.5)
axis(2,pos=-0.07,cex.axis=1.5)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}
else {
plot(u,G,type="l",col=1,ylim=c(mi,ma),xlab="",ylab="",lwd=3,yaxt='n',xaxt='n')
abline(h=0)
axis(1,pos=mi-0.05,cex.axis=1.5)
text(0.5,0.6*ma,cex=2,bquote(A[.(j)*','*.(k)]))}}
title(paste("N=",N,",  T=",T,",  r=",r),outer=TRUE,cex.main=1.6)

#dev.print(pdf,file="Ax-coal-talk.pdf")
dev.off()

# L are the smooth loadings

V=matrix(0,T,r)
P=array(0,c(N,r,T))	
for (t in 1:T){
	V[t,]=eigen(S[,,t],symmetric=TRUE)$values[1:r]
	P[,,t]=eigen(S[,,t],symmetric=TRUE)$vectors[,1:r]}

#P=P*sqrt(N)

# P are the eigenvectors of the smooth matrix S

c=5

###########
#matching the eigenvectors
###########

Q1=P[,,(1:(T/2))]			# first half: identified up-to-sign
Q2=array(0,c(N,r,T/2))		 

Q2[,1,]=P[,2,((T/2 + 1):T)] 	# second half: not-aligned and identified up-to-sign
Q2[,2,]=P[,1,((T/2 + 1):T)]

R1=Q1
R2=Q2

for (k in 1:r){
Q=matrix(Q1[,k,1:(T/2 - c)],N,T/2 -c) # "sign-adjusted" eigenvector: N x T/2
for (t in 1:(T/2 - c -1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
#b=matrix(Q[,t+1]+Q[,t],N,1)
#if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}
if (norm(a,"F")> sqrt(2)){Q[,t+1] = - Q[,t+1]}
}
R1[,k,1:(T/2 - c)]=Q
}

for (k in 1:r){
Q=matrix(Q2[,k,((c+1):(T/2))],N,T/2 -c) # "signed-adjusted" eigenvector: N x T/2
for (t in 1:(T/2 - c -1)){
a=matrix(Q[,t+1]-Q[,t],N,1)
#b=matrix(Q[,t+1]+Q[,t],N,1)
#if (norm(a,"F")> norm(b,"F")){Q[,t+1] = - Q[,t+1]}
if (norm(a,"F")>sqrt(2)){Q[,t+1] = - Q[,t+1]}
}
R2[,k,((c+1):(T/2))]=Q
}

R=array(0,c(N,r,T))	
for (k in 1:r){
    R[,k,1:(T/2)]=R1[,k,]
R[,k,((T/2)+1):T]=R2[,k,]
}

Q=R
for (k in 1:r){
a=matrix(R1[,k,T/2-c],N,1)
b=matrix(R2[,k,c+1],N,1)
#if (norm(a-b,"F")> norm(a+b,"F")){Q[,k,(T/2+1):T]=-R2[,k,]}
if (norm(a-b,"F")> sqrt(2)){Q[,k,(T/2+1):T]=-R2[,k,]}
}

uu=u[-T]
mi=min(L) 
ma=max(L) 


NN=N
# b l t r
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
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)])))
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
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)])))
}

#dev.print(pdf,file="Pi-talk.pdf")
dev.off()

# b l t r
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
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,R[.(j)*','*.(k)])))
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
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,R[.(j)*','*.(k)])))
}

#dev.print(pdf,file="Phat-talk.pdf")
dev.off()

# b l t r
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,abs(P)[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,abs(P)[j,k,1:(T-1)],type="l",col=1,lwd=1.5, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(abs(Pi)[.(j)*','*.(k)], "  vs  " ,abs(R)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,abs(P)[j,k,1:(T-1)],type="l",col=1,lwd=1.5,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,abs(P)[j,k,1:(T-1)],type="l",col=1,lwd=1.5,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,abs(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(abs(Pi)[.(j)*','*.(k)], "  vs  " ,abs(R)[.(j)*','*.(k)])))}

#dev.print(pdf,file="Pabs-talk.pdf")
dev.off()


# b l t r
x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,R[j,k,1:(T-1)],type="l",col=1,lwd=2,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=2,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(P)[.(j)*','*.(k)])))
}


#dev.print(pdf,file="Pbar-talk.pdf")
dev.off()

Q=-Q

x11()
setpar(mfrow=c(NN,r),oma=c(1.2,3,1,0.5),mar=c(.8,.5,1,0.8))
for (j in 1:(NN-1)){
for (k in 1:r){
if (k==1){
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=2, ylim=c(mi,ma), xlab="",ylab="",xaxt="none",yaxt='n',cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=2, ylim=c(mi,ma),xlab="",ylab="",yaxt="none",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(Pi)[.(j)*','*.(k)])))
}}
j=NN
for (k in 1:r){
if (k==1){
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=2,xaxt="none", yaxt="none", ylim=c(mi,ma),xlab="",ylab="",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
axis(2,pos=-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
else{
plot(uu,Q[j,k,1:(T-1)],type="l",col=1,lwd=2,yaxt="none", ylim=c(mi,ma),xlab="",ylab="",xaxt="none",cex.axis=1.6)
lines(u,(L[j,k,]),col=2,lwd=1,lty=1)
axis(1,pos=mi-0.05,cex.axis=1.6)
abline(h=0)
abline(v=0.5)}
text(0.2,0.75*ma,cex=2,bquote(paste(Pi[.(j)*','*.(k)], "  vs  " ,bar(Pi)[.(j)*','*.(k)])))
}

#dev.print(pdf,file="Ptilde-talk.pdf")
dev.off()
#############################################################
#############################################################

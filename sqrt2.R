
############################
# the sqrt2-idea
# Giovanni Motta 
# Copyright © 2023
############################

rm(list=ls())
d2=seq(0,2,length.out=200)
d1=2*sqrt(1 - (d2^2)/4)
#library(fontcm)
library(showtext)

par(oma=c(0,0,1,1),mar=c(5,5,1.5,1))
plot(d2,d1,lwd=2,type="l",xlab=bquote(paste("||",P(t+1)-P(t),"||")), ylab=bquote(paste("||",P(t+1)+P(t),"||")),xaxt='n',yaxt='n',cex.lab=1.5,main=bquote(paste("the ",sqrt(2),"-idea")))
axis(side=1, at=c(0,sqrt(2),2), labels=expression(0,sqrt(2),2),cex.axis=1.5,tck=0.02)
axis(side=2, at=c(0,sqrt(2),2), labels=expression(0,sqrt(2),2),cex.axis=1.5,tck=0.02)
abline(h=sqrt(2),lty=2) #sqrt-2
abline(v=sqrt(2),lty=2) #is a uniform bound
text(x = 0.5*(2+sqrt(2)), y = 0, '}', srt = 90, cex = 10, family = 'Helvetica Neue UltraLight')
text(x = 0.5*(2+sqrt(2)), y = 0.1, bquote(paste(tilde(P)(t+1),"=", -P(t+1))),cex=1.1)

vectors <- function(X, origin=c(0,0),
                    lwd=2, angle=13, length=0.15,
                    labels=TRUE, cex.lab=1.5, pos.lab=4, frac.lab=1,  ...) {
  if (is.vector(X)) X <- matrix(X, ncol=2)
  arrows(origin[1], origin[2], X[,1], X[,2], lwd=lwd, length=length, angle=angle, ...)
  if (is.logical(labels) && labels) {
    labels <- rownames(X)
  }
  if (!is.null(labels)) {
    # DONE: allow for labels to be positioned some fraction of the way from origin to X
    # FIXME: it is dangerous to use ... for both arrows() and text(), e.g., for col=
    xl = origin[1] + frac.lab * (X[,1]-origin[1])
    yl = origin[2] + frac.lab * (X[,2]-origin[2])
    text(xl, yl, labels, cex=cex.lab, pos=pos.lab, ...)
  }
}


u <- c(3,1)
u=u/sqrt(sum(u^2))
v <- c(1,3)
v=v/sqrt(sum(v^2))

######################################
sum <- u+v
diff=v-u
xlim <- c(-1.25,1.25)
ylim <- c(-1.25,1.25)
# proper geometry requires asp=1
x11()
plot(xlim, ylim, type="n", xlab="", ylab="", asp=1)
abline(v=0, h=0, col="gray")
vectors(rbind(u,v,`u+v`=sum,`v-u`=diff), 
labels=c("P(t)","P(t+1)","P(t+1)+P(t)","P(t+1)-P(t)"),  col=c("black", "black", "grey","red"),lty=c(1,1,2,2),
cex.lab=c(2, 2, 2.2))
# show the opposing sides of the parallelogram
vectors(sum, origin=u, col="black", lty=2)
vectors(sum, origin=v, col="black", lty=2)
vectors(v, origin=u, col="red", lty=2) # do NOT change the sign
#dev.print(pdf,file="parall1.pdf")



######################################
z=-v
sum <- u+z
diff=z-u
x11()
plot(xlim, ylim, type="n", xlab="", ylab="", asp=1)
abline(v=0, h=0, col="gray")
vectors(rbind(u,z,`u+z`=sum,`z-u`=diff), 
labels=c("P(t)","P(t+1)","P(t+1)+P(t)","P(t+1)-P(t)"),  col=c("black", "black", "grey","red"),lty=c(1,1,2,2),
cex.lab=c(2, 2, 2.2))

# show the opposing sides of the parallelogram
vectors(sum, origin=u, col="black", lty=2)
vectors(sum, origin=z, col="black", lty=2)
vectors(z, origin=u, col="red", lty=2) # change the sign
#dev.print(pdf,file="parall2.pdf")
######################################


w=-z
sum <- u+w
diff=w-u
x11()
plot(xlim, ylim, type="n", xlab="", ylab="", asp=1)
abline(v=0, h=0, col="gray")
vectors(rbind(u,z,w,`w+u`=sum,`w-v`=diff), 
labels=c("P(t)","P(t+1)","-P(t+1)","-P(t+1)+P(t)","-P(t+1)-P(t)"),  col=c("black", "black","black","grey","red"),lty=c(1,1,1,2,2),
cex.lab=c(2, 2, 2.2))
# show the opposing sides of the parallelogram
vectors(sum, origin=u, col="black", lty=2)
vectors(sum, origin=w, col="black", lty=2)
vectors(w, origin=u, col="red", lty=2) # change-of-sign geometry 
#dev.print(pdf,file="parall3.pdf")






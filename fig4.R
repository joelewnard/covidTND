wu = c(0.85,0.7,0.55,0.4)
beta1 = c(2,4)
beta2 = c(0.25,0.5)

trueEst = 10^seq(-3,0,0.01)

or1 = or2 = array(NA,dim=c(length(wu),2,length(trueEst)))
for (i in 1:length(wu)){
  for (j in 1:2){
    wv1 = exp(-beta1[j]*log(wu[i]))    
    wv2 = exp(-beta2[j]*log(wu[i]))    
    or1[i,j,] = trueEst*beta1[j]*wv1/wu[i]
    or2[i,j,] = trueEst*beta2[j]*wv2/wu[i]
  }
}


plot(or1[1,2,],type='l',ylim=c(0,2))
lines(trueEst,col='red')



pdf(file='fig4.pdf',width=6.5,height=2.25)
par(mfrow=c(1,3))

par(lwd=0.5)
par(tck=-0.02)
par(mgp=c(3,0.35,0))
par(mar=c(3.5,3,1,1))
cols = c('darkorchid4','darkslategray4','darkorange','forestgreen')
ltys = c(1,5)

plot(1,type='n',axes=F,ann=F,ylim=c(-7.5,7.5),xlim=c(-7,0))
lines(y=log(trueEst),x=log(trueEst),col='red')
for (i in 1:length(wu)) for (j in 1:2){
  lines(log(or1[i,j,]),x=log(trueEst),col=cols[i],lty=ltys[j])
}
box(bty='l')
axis(side=2,at=log(c(0.001,0.01,0.1,1,10,100,1000)),
     labels=c(expression(10^-3),expression(10^-2),expression(10^-1),1,
              expression(10^1),expression(10^2),expression(10^3)),cex.axis=0.65,las=1,
     lwd=0,lwd.ticks=0.5)
axis(side=1,at=log(c(0.001,0.01,0.1,1)),c(expression(10^-3),expression(10^-2),expression(10^-1),1),
     lwd=0,lwd.ticks=0.5,labels=NA)
text(x=log(c(0.001,0.01,0.1,1)),y=-8.8,adj=1,xpd=T,srt=45,
     c(expression(10^-3),expression(10^-2),expression(10^-1),1),cex=0.65)
mtext(side=2,expression(hat(OR)),cex=0.6,line=1.75,las=1)
mtext(side=1,expression(1-theta[italic(S)]*theta[italic(P)]),cex=0.65,line=2)
text(x=c(-6.6,-5.6),y=log(10^3),
     c(expression(beta==2),expression(beta==4)),cex=0.75,xpd=T,srt=45,adj=0)
for (i in 1:length(wu)){
  lines(y=rep(log(10^(3-(i-1)*0.4)),2)-0.75,x=c(-7,-6.2),col=cols[i])
  lines(y=rep(log(10^(3-(i-1)*0.4)),2)-0.75,x=c(-6,-5.2),col=cols[i],lty=5)
}
text(x=-5,y=log(10^(3-c(0,1,2,3)*0.4))-0.75,
     c(expression(omega[italic(U)]==0.9),
       expression(omega[italic(U)]==0.7),
       expression(omega[italic(U)]==0.5),
       expression(omega[italic(U)]==0.3)),cex=0.75,adj=0)




plot(1,type='n',axes=F,ann=F,ylim=c(-7.5,7.5),xlim=c(-7,0))
lines(y=log(trueEst),x=log(trueEst),col='red')
for (i in 1:length(wu)) for (j in 1:2){
  lines(log(or2[i,j,]),x=log(trueEst),col=cols[i],lty=ltys[j])
}
box(bty='l')
axis(side=2,at=log(c(0.001,0.01,0.1,1,10,100,1000)),
     labels=c(expression(10^-3),expression(10^-2),expression(10^-1),1,
              expression(10^1),expression(10^2),expression(10^3)),cex.axis=0.65,las=1,
     lwd=0,lwd.ticks=0.5)
axis(side=1,at=log(c(0.001,0.01,0.1,1)),c(expression(10^-3),expression(10^-2),expression(10^-1),1),
     lwd=0,lwd.ticks=0.5,labels=NA)
text(x=log(c(0.001,0.01,0.1,1)),y=-8.8,adj=1,xpd=T,srt=45,
     c(expression(10^-3),expression(10^-2),expression(10^-1),1),cex=0.65)
mtext(side=2,expression(hat(OR)),cex=0.6,line=1.75,las=1)
mtext(side=1,expression(1-theta[italic(S)]*theta[italic(P)]),cex=0.65,line=2)

text(x=c(-6.6,-5.6),y=log(10^3),
     c(expression(beta==0.25),expression(beta==0.5)),cex=0.75,xpd=T,srt=45,adj=0)
for (i in 1:length(wu)){
  lines(y=rep(log(10^(3-(i-1)*0.4)),2)-0.75,x=c(-7,-6.2),col=cols[i])
  lines(y=rep(log(10^(3-(i-1)*0.4)),2)-0.75,x=c(-6,-5.2),col=cols[i],lty=5)
}
text(x=-5,y=log(10^(3-c(0,1,2,3)*0.4))-0.75,
     c(expression(omega[italic(U)]==0.9),
       expression(omega[italic(U)]==0.7),
       expression(omega[italic(U)]==0.5),
       expression(omega[italic(U)]==0.3)),cex=0.75,adj=0)



plot(1,xlim=c(-5,5),ylim=c(0,1),axes=F,ann=F,type='n')
xs = exp(seq(-5,5,0.01)); xs[xs==1] = NA
lines(y=xs^(1/(1-xs)),x=log(xs))
box(bty='l')
mtext(side=2,expression(omega[U]*'*'),cex=0.6,line=1.75,las=1)
mtext(side=1,expression(beta),cex=0.65,line=2)
axis(side=1,at=log(c(0.01,0.1,1,10,100)),
     lwd=0,lwd.ticks=0.5,labels=NA)
text(x=log(c(0.01,0.1,1,10,100)),y=-0.09,adj=1,xpd=T,srt=45,
     c(expression(10^-2),expression(10^-1),1,
       expression(10^1),expression(10^2)),cex=0.65)
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
dev.off()






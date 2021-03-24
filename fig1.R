theta1 = c(1,0.5,0.2)
theta2 = seq(0.02,1,0.01)



ks = c(0.184,0.655*897961/330e6,0.682*103900/330e6)
### prev any symptom among uninfected in fw study (https://www.medrxiv.org/content/10.1101/2020.12.27.20248894v1.full.pdf+html)
### hosp beds at normal 65.5% occupancy (https://www.cdc.gov/nchs/data/hus/2017/089.pdf)
### icu beds at normal 68.2% occupany (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5520980/)

ds = c(0.5,0.029,0.029*0.407)
### prob any symptoms (assume 50%)
### prob hospital (salje science 2.9%)
### prob ICU (lewnard bmj 40.7% given hosp)

as = c(0.01,0.1)


### index: shed/prog, as, k/ds, thetas
estsTND1 = estsTCC1 = estsTND2 = estsTCC2 = array(NA,dim=c(3,3,3,length(theta2)))

for (i in 1:3){
   for (j in 1:2){
     for (k in 1:3){
       for (l in 1:dim(estsTND1)[4]){
 
         thetaP1 = theta1[i]
         thetaS1 = theta2[l]
         
         thetaP2 = theta2[l]
         thetaS2 = theta1[i]
         
         a = as[j]; d = ds[k]; n = ks[k]
         
         estsTCC1[i,j,k,l] = 1 - thetaS1*(thetaP1*d + (1-thetaP1*d)*n)*(1-a*d)/((1-thetaS1*thetaP1*a*d)*(d+(1-d)*n))
         estsTCC2[i,j,k,l] = 1 - thetaS2*(thetaP2*d + (1-thetaP2*d)*n)*(1-a*d)/((1-thetaS2*thetaP2*a*d)*(d+(1-d)*n))
         
         estsTND1[i,j,k,l] = 1 - thetaS1*(thetaP1*d + (1-thetaP1*d)*n)*(1-a)/((1-thetaS1*a)*(d+(1-d)*n))
         estsTND2[i,j,k,l] = 1 - thetaS2*(thetaP2*d + (1-thetaP2*d)*n)*(1-a)/((1-thetaS2*a)*(d+(1-d)*n))
       }
     }
   }
}


trueVE = array(NA,dim=c(3,dim(estsTND1)[4]))
for (i in 1:3){
  trueVE[i,] = 1 - theta1[i]*theta2
}




setwd('~/Dropbox/covid tnd')

pdf(file='fig1.pdf',height=4,width=6.5)
layout(matrix(c(rep(14,3),18,rep(16,3),
              1,3,5,18,2,4,6,
              rep(15,3),18,rep(17,3),
              7,9,11,18,8,10,12,
              rep(13,7)),nrow=5,ncol=7,byrow=T),heights=c(0.1,1,0.1,1,0.3),widths=c(rep(1,3),0.25,rep(1,3)))
xs = 1-theta2
par(lwd=0.5)
par(tck=-0.04)
par(mgp=c(3,0.3,0))
par(mar=c(3.5,3,1,0.5))
cols = c('darkorchid4','darkslategray4','goldenrod3')

plot.fn = function(objTCC,objTND,lab,i,j,ymin,ymax){
  plot(y=trueVE[i,],x=log(1-xs),type='n',axes=F,ann=F,ylim=c(ymin,ymax))
  lines(y=log(1-trueVE[i,]),x=log(1-xs),col='red')
  for (k in 1:3){
    lines(y=objTCC[i,j,k,],x=log(1-xs),col=cols[k])
    lines(y=objTND[i,j,k,],x=log(1-xs),col=cols[k],lty=5)
  }
  if (lab=='S'){
    if (i==1){
      mtext(expression(theta[italic(S)]==1),adj=0,cex=0.6)    
    } else{
      if (i==2){
        mtext(expression(theta[italic(S)]==0.5),adj=0,cex=0.6)
      } else{
        mtext(expression(theta[italic(S)]==0.2),adj=0,cex=0.6)
      }
    }
    mtext(side=1,line=1.5,expression(theta[P]),cex=0.6)
  } else{
    if (i==1){
      mtext(expression(theta[italic(P)]==1),adj=0,cex=0.6)    
    } else{
      if (i==2){
        mtext(expression(theta[italic(P)]==0.5),adj=0,cex=0.6)
      } else{
        mtext(expression(theta[italic(P)]==0.2),adj=0,cex=0.6)
      }
    }
    mtext(side=1,line=1.5,expression(theta[S]),cex=0.6)
  }
  box(bty='l')
  axis(side=1,at=log(c(0.01,0.02,0.05,0.1,0.2,0.5,1)),labels=NA,lwd=0,lwd.ticks=0.5)
  axis(side=2,at=log(c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1)),
       labels=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1),
       lwd=0,lwd.ticks=0.5,cex.axis=0.6,las=1)
  text(x=log(c(0.02,0.05,0.1,0.2,0.5,1)),
       c(0.02,0.05,0.1,0.2,0.5,1),y=1.1*ymin,adj=1,srt=45,xpd=T,cex=0.6)
  #mtext(side=2,expression(paste(1-hat(VE),', %',sep='')),cex=0.6,line=1.25)
  mtext(side=2,expression(hat(OR)),cex=0.6,line=1.75,las=1)
}

for (j in 1:2){
  for (i in 1:3){
    plot.fn(objTCC=log(1-estsTCC1),objTND=log(1-estsTND1),i=i,j=j,lab='S',ymin=-5.5,ymax=0)
    plot.fn(objTCC=log(1-estsTCC2),objTND=log(1-estsTND2),i=i,j=j,lab='P',ymin=-5.5,ymax=0)
  }
}

par(mar=c(0,0,0,0))
plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
ys = c(0.6,0.35,0.1)
for (i in 1:3){
  lines(x=c(0,0.05),y=rep(ys[i],2),col=cols[i])
  lines(x=c(0.075,0.125),y=rep(ys[i],2),col=cols[i],lty=5)
}
text(x=c(0.025,0.1),y=0.85,c('TCC','TND'),cex=0.7)
text(x=0.15,y=ys,c('Any symptoms','Hospitalization','ICU admission'),adj=0,cex=0.7)
text(x=0.275,y=ys+0.01,c(expression(italic(d)==0.5),expression(italic(d)==0.029),expression(italic(d)==0.012)),adj=0,cex=0.7)
text(x=0.375,y=ys+0.01,c(expression(italic(k)==0.184),expression(italic(k)==0.0018),expression(italic(k)==0.00021)),adj=0,cex=0.7)

for (i in 1:2){
  plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
  text(x=-0.03,y=0.5,expression(paste('Low prevalence of infection, ',italic(a)==0.01,sep='')),adj=0,cex=1,xpd=T)
  plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
  text(x=-0.03,y=0.5,expression(paste('High prevalence of infection, ',italic(a)==0.1,sep='')),adj=0,cex=1,xpd=T)
}

dev.off()







pdf(file='fig s1.pdf',height=4,width=6.5)
layout(matrix(c(rep(14,3),18,rep(16,3),
                1,3,5,18,2,4,6,
                rep(15,3),18,rep(17,3),
                7,9,11,18,8,10,12,
                rep(13,7)),nrow=5,ncol=7,byrow=T),heights=c(0.1,1,0.1,1,0.3),widths=c(rep(1,3),0.25,rep(1,3)))
xs = 1-theta2
par(lwd=0.5)
par(tck=-0.04)
par(mgp=c(3,0.3,0))
par(mar=c(3.5,3,1,0.5))
cols = c('darkorchid4','darkslategray4','goldenrod3')

plot.fn = function(objTCC,objTND,lab,i,j,ymin,ymax){
  plot(y=trueVE[i,],x=(1-xs),type='n',axes=F,ann=F,ylim=c(ymin,ymax))
  lines(y=(1-trueVE[i,]),x=(1-xs),col='red')
  for (k in 1:3){
    lines(y=objTCC[i,j,k,],x=(1-xs),col=cols[k])
    lines(y=objTND[i,j,k,],x=(1-xs),col=cols[k],lty=5)
  }
  if (lab=='S'){
    if (i==1){
      mtext(expression(theta[italic(S)]==1),adj=0,cex=0.6)    
    } else{
      if (i==2){
        mtext(expression(theta[italic(S)]==0.5),adj=0,cex=0.6)
      } else{
        mtext(expression(theta[italic(S)]==0.2),adj=0,cex=0.6)
      }
    }
    mtext(side=1,line=1.5,expression(theta[P]),cex=0.6)
  } else{
    if (i==1){
      mtext(expression(theta[italic(P)]==1),adj=0,cex=0.6)    
    } else{
      if (i==2){
        mtext(expression(theta[italic(P)]==0.5),adj=0,cex=0.6)
      } else{
        mtext(expression(theta[italic(P)]==0.2),adj=0,cex=0.6)
      }
    }
    mtext(side=1,line=1.5,expression(theta[S]),cex=0.6)
  }
  box(bty='l')
  axis(side=1,at=seq(0,1,0.2),labels=NA,lwd=0,lwd.ticks=0.5)
  axis(side=2,at=seq(0,1,0.1),
       labels=seq(0,1,0.1),
       lwd=0,lwd.ticks=0.5,cex.axis=0.6,las=1)
  text(x=seq(0,1,0.2),
       seq(0,1,0.2),y=ymin-0.1*ymax,adj=1,srt=45,xpd=T,cex=0.6)
  #mtext(side=2,expression(paste(1-hat(VE),', %',sep='')),cex=0.6,line=1.25)
  mtext(side=2,expression(hat(OR)),cex=0.6,line=1.75,las=1)
}

for (j in 1:2){
  for (i in 1:3){
    plot.fn(objTCC=(1-estsTCC1),objTND=(1-estsTND1),i=i,j=j,lab='S',ymin=0,ymax=1)
    plot.fn(objTCC=(1-estsTCC2),objTND=(1-estsTND2),i=i,j=j,lab='P',ymin=0,ymax=1)
  }
}

par(mar=c(0,0,0,0))
plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
ys = c(0.6,0.35,0.1)
for (i in 1:3){
  lines(x=c(0,0.05),y=rep(ys[i],2),col=cols[i])
  lines(x=c(0.075,0.125),y=rep(ys[i],2),col=cols[i],lty=5)
}
text(x=c(0.025,0.1),y=0.85,c('TCC','TND'),cex=0.7)
text(x=0.15,y=ys,c('Any symptoms','Hospitalization','ICU admission'),adj=0,cex=0.7)
text(x=0.275,y=ys+0.01,c(expression(italic(d)==0.5),expression(italic(d)==0.029),expression(italic(d)==0.012)),adj=0,cex=0.7)
text(x=0.375,y=ys+0.01,c(expression(italic(k)==0.184),expression(italic(k)==0.0018),expression(italic(k)==0.00021)),adj=0,cex=0.7)

for (i in 1:2){
  plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
  text(x=-0.03,y=0.5,expression(paste('Low prevalence of infection, ',italic(a)==0.01,sep='')),adj=0,cex=1,xpd=T)
  plot(1,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(0,1))
  text(x=-0.03,y=0.5,expression(paste('High prevalence of infection, ',italic(a)==0.1,sep='')),adj=0,cex=1,xpd=T)
}

dev.off()

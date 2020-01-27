########Models evaluating whether older males increase their probability of siring the next infant by being a primary associate (mating effort)

#ID: male ID
#ID_score: Elo score
#age: male age
#ageest: male age known or estimated
#PHG: group membership variables (ENK is reference)

#####################
library(rethinking)

d$s_age=(d$prtage-mean(d$prtage))/sd(d$prtage)
d$agesq=d$s_age^2
d$dyad=paste(d$male,d$mom)
d$dyad=droplevels(as.factor(d$dyad))
d$dyad=as.integer(d$dyad)

d$male_index=droplevels(as.factor(d$male))
d$male_index=as.integer(as.factor(d$male))
d$mom_index=droplevels(as.factor(d$mom))
d$mom_index=as.integer(as.factor(d$mom))

#Model code for model without interaction
nextdad<- map2stan(  
  alist(
    nextdad ~ dbinom( 1 , p ) ,
    
    logit(p) ~  a  + b_age*s_age + b_agesq*agesq  + b_phg*PHG + b_pa*primassoc +
				m_id[male_index] + sigma_fid*f_id[mom_index] + sigma_dyad*dyad_id[dyad] +
      
				ba[male_index]*s_age + bp[male_index]*primassoc + basq[male_index]*agesq ,
    
    c(a,b_age,b_phg,b_pa,b_agesq) ~ dnorm(0,2),
    c(m_id,ba,bp,basq)[male_index] ~ dmvnormNC( sigma_id, Rho_id ),
    f_id[mom_index] ~ dnorm(0,2),
    dyad_id[dyad] ~ dnorm(0,2),
    c(Rho_id) ~ dlkjcorr(3),
    c(sigma_id,sigma_dyad,sigma_fid) ~ dcauchy(0,2)
  ),
  data=d, cores=3, chains=4, warmup=3500, iter=7000, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)
  
#Model code for model with interaction 
 nextdadX <-  map2stan(
	alist(
		nextdad ~ dbinom(1, p),
		logit(p) ~ a +  b_age * s_age + b_agesq * agesq + b_phg * PHG + b_pa * primassoc + 
					b_paa * s_age * primassoc + b_paagesq * agesq * primassoc + 
					m_id[male_index] + sigma_fid*f_id[mom_index] + sigma_dyad * dyad_id[dyad] +
					ba[male_index]*s_age + bp[male_index]*primassoc + basq[male_index]*agesq + bpa[male_index]*s_age*primassoc + bpaasq[male_index]*agesq*primassoc,
	
		c(a,b_age, b_phg, b_pa, b_agesq, b_paa, b_paagesq) ~ dnorm(0,2),
		c(m_id, ba, bp, bpa, basq, bpaasq)[male_index] ~ dmvnormNC(sigma_id,Rho_id),
		f_id[mom_index] ~ dnorm(0, 2),
		dyad_id[dyad] ~ dnorm(0, 2),
		c(Rho_id) ~ dlkjcorr(3),
		c(sigma_id, sigma_dyad,sigma_fid) ~ dcauchy(0, 2)
		),	
		data = d, iter = 7000, warmup = 3500, chains = 4, WAIC = TRUE, cores = 3, control = list(adapt_delta = 0.99,max_treedepth = 15))

compare(nextdad,nextdadX)  

####Construct posterior predictions  
a_mom_z <- matrix(0,1000,length(unique(d$mom_index)))
a_male_z<- matrix(0,1000,length(unique(d$male_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad)))

age.seq=seq(min(d$s_age),max(d$s_age),length=200)

pred0 <- list(
  mom_index=rep(1,length(age.seq)),
  male_index=rep(1,length(age.seq)),
  dyad=rep(1,length(age.seq)),
  primassoc=rep(0,length(age.seq)),
  s_age=age.seq,
  agesq=age.seq^2,
  PHG=rep(mean(d$PHG),length(age.seq))
)	

pred1 <- list(
  mom_index=rep(1,length(age.seq)),
  male_index=rep(1,length(age.seq)),
  dyad=rep(1,length(age.seq)),
  primassoc=rep(1,length(age.seq)),
  s_age=age.seq,
  agesq=age.seq^2,
  PHG=rep(mean(d$PHG),length(age.seq))
)	

#Model averaging
Pa0=ensemble(nextdadX,nextdad,data=pred0,replace= list(f_id=a_mom_z,dyad_id=a_dyad_z, m_id=a_male_z,ba=a_male_z,bp=a_male_z,bpa=a_male_z,basq=a_male_z, bpaasq=a_male_z))
Pa1=ensemble(nextdadX,nextdad,data=pred1,replace= list(f_id=a_mom_z,dyad_id=a_dyad_z, m_id=a_male_z,ba=a_male_z,bp=a_male_z,bpa=a_male_z,basq=a_male_z, bpaasq=a_male_z))


median0=apply(Pa0$link , 2 , median )
HPDI0=apply(Pa0$link , 2 , HPDI )
median1=apply(Pa1$link , 2 , median )
HPDI1=apply(Pa1$link , 2 , HPDI )

#############Plotting
yy0=as.vector(t(Pa0$link[1:1000,]))
yy1=as.vector(t(Pa1$link[1:1000,]))

yy0=as.vector(t(Pa0[1:1000,]))
yy1=as.vector(t(Pa1[1:1000,]))


par(mar=c(2.3,2.3,0.5,0.5))

smoothScatter(rep(age.seq,1000),yy1,xlim=c(min(d$s_age),max(d$s_age)),ylim=c(0,1),bandwidth=0.01,colramp =colorRampPalette(c(adjustcolor(rgb(1, 1, 1, 0),alpha.f=0.5), adjustcolor(rgb(0, 0, 255,maxColorValue=255),alpha.f=0.5))),nbin=500,transformation = function(x) x^.7,ylab='',xlab='',cex=1.2,xaxt='n',yaxt='n',yaxs='r',nrpoints=0)

smoothScatter(rep(age.seq,1000),yy0,xlim=c(min(d$s_age),max(d$s_age)),ylim=c(0,1),bandwidth=0.01,colramp =colorRampPalette(c(adjustcolor(rgb(1, 1, 1, 0),alpha.f=0.5), adjustcolor(rgb(255, 230, 0, maxColorValue=255),alpha.f=0.5)),alpha=TRUE),nbin=500,transformation = function(x) x^.4,ylab='',xlab='',cex=1.2,xaxt='n',yaxt='n',yaxs='r',nrpoints=0,add=TRUE)

lines( age.seq ,median1,lwd=2,lty=1,col="blue")
lines(age.seq,HPDI1[1,],lty=2,lwd=1,col="blue")
lines(age.seq,HPDI1[2,],lty=2,lwd=1,col="blue")
lines( age.seq ,median0,lwd=2,lty=1,col="orange")
lines(age.seq,HPDI0[1,],lty=2,lwd=1,col="orange")
lines(age.seq,HPDI0[2,],lty=2,lwd=1,col="orange")

points(d$s_age[d$primassoc==1&d$nextdad==1],rep(1.027,length(d$s_age[d$primassoc==1&d$nextdad==1])), col=alpha("blue",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$primassoc==1&d$nextdad==0],rep(-.0277,length(d$s_age[d$primassoc==1&d$nextdad==0])), col=alpha("blue",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$primassoc==0&d$nextdad==1],rep(1.027,length(d$s_age[d$primassoc==0&d$nextdad==1])) , col=alpha("orange",0.8),pch=-as.hexmode("007C"),cex=0.25)
points(d$s_age[d$primassoc==0&d$nextdad==0],rep(-0.0277,length(d$s_age[d$primassoc==0&d$nextdad==0])) , col=alpha("orange",0.8),pch=-as.hexmode("007C"),cex=0.25)

lab.a=seq(6,16,by=2)
lab.sa=(lab.a-mean(d$prtage))/sd(d$prtage)
axis(1,labels=NA,at=lab.sa,tck=-0.01)
mtext(lab.a,at=lab.sa,line=0.1,side=1,cex=0.7)

lab.s=pretty(c(0,1))
axis(2,labels=NA,at=lab.s,tck=-0.01)
mtext(lab.s,at=lab.s,line=0.2,side=2,cex=0.7)
mtext("Male age", side=1, line=1, cex=0.8,tck=-0.01)
mtext("Probability of siring the next offspring", side=2, line=1.1, cex=0.8)

legend(x=0.3,y=0.9, legend = c("Primary associates","Other males"),
       col=c(col.alpha("blue", 1) ,col.alpha("gold", 1) ) ,
       pt.cex=2 , bty="n",y.intersp=1.3,x.intersp=1, lty=c(1,1) ,cex=0.7,lwd=c(2,2))


#######Model assessing the relationship of age and rank##########

#ID: male ID
#ID_score: Elo_score
#PHG: group membership variable (ENK reference)
#age: male age
#ageest: male age estimated or known

#########################
d$male=as.integer(as.factor(d$ID))
d$s_age=(d$age-mean(d$age))/sd(d$age)
d$agesq=d$s_age^2

library(rethinking)

########'Null' model code
model0<- map2stan(  

alist(
ID_score ~ dnorm(mu,sigma),

	mu <- a + a_mid[male] + b_phg*PHG ,
	
	a ~ dnorm(0,100),
	b_phg ~ dnorm(0,100),
	a_mid[male] ~ dnorm(0,sigma_ID),
	sigma_ID ~ dcauchy(0,2),
	sigma ~ dunif(0,100)
	),

data=d, cores=2 , chains=4, warmup=3500, iter=7000,control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 

######Model code
rank_model<- map2stan(  

alist(
ID_score ~ dnorm(mu,sigma),

	mu <- a + a_mid[male] + (b_age + b_age_m[male])*s_age + (b_agesq + b_agesq_m[male])*agesq + b_phg*PHG ,
	
	a ~ dnorm(0,100),
	c(b_age,b_agesq,b_phg) ~ dnorm(0,100),
	c(a_mid,b_age_m,b_agesq_m)[male] ~ dmvnormNC(sigma_ID,Rho),
	sigma_ID ~ dcauchy(0,2),
	Rho ~ dlkjcorr(3),
	sigma ~ dunif(0,100)
	),

data=d, cores=2 , chains=4, warmup=3500, iter=7000,control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 

res=precis(rank_model, depth=2 , digits=2)@output
plot(rank_model)
compare(model0,rank_model)

##########Construct posterior predictions
a_male_z <- matrix(0,1000,length(unique(d$male)))

age.seq=seq(min(d$s_age),max(d$s_age),length=100)

pred_a1 <- list(
	male=rep(1,length(age.seq)),
    agesq=age.seq^2,
	s_age=age.seq,
	PHG=rep(mean(d$PHG),length(age.seq))
	)

RMa1 <- link(rank_model, n=1000 , data=pred_a1, replace=list(a_mid=a_male_z,b_age_m=a_male_z,b_agesq_m=a_male_z), WAIC=TRUE)

pred.mean.a=apply(RMa1 , 2 , mean )
pred.HPDI.a=apply(RMa1 , 2 , HPDI )


#########Plotting
xx=age.seq
yy=as.vector(t(RMa1[1:1000,]))
par(mar=c(2.3,2.3,0.5,0.5))
smoothScatter(rep(xx,1000),yy,xlim=c(min(d$s_age),max(d$s_age)),ylim=c(min(d$ID_score),max(d$ID_score)),colramp = colorRampPalette(c("white", "steelblue1")),nbin=200,transformation = function(x) x^.7,ylab='',xlab='',cex=1.2,xaxt='n',yaxt='n',yaxs='r',nrpoints=0)
points(ID_score[d$ageest=="known"]~ s_age[d$ageest=="known"] , data=d , col=alpha("blue",0.6),pch=16 ,cex=0.9)
points(ID_score[d$ageest=="estimate"] ~ s_age[d$ageest=="estimate"]  , data=d , col=alpha("darkturquoise",0.6),pch=16,cex=0.9)
lines( age.seq , pred.mean.a ,lwd=1)
lines(age.seq[-(99:100)],pred.HPDI.a[1,-(99:100)],lty=2,lwd=0.7)
lines(age.seq,pred.HPDI.a[2,],lty=2,lwd=0.7)

lab.a=seq(6,18,by=2)
lab.sa=(lab.a-mean(d$age))/sd(d$age)
axis(1,labels=NA,at=lab.sa,tck=-0.01)
mtext(lab.a,at=lab.sa,line=0.1,side=1,cex=0.7)
lab.r=seq(-200,300,by=100)

axis(2,labels=NA,at=lab.r,tck=-0.01)
mtext(side=2,text=lab.r,at=lab.r,line=0.2,cex=0.7)
mtext("Male age", side=1,  cex=0.8,line=1)
mtext("Male rank (Elo-score)", side=2, cex=0.8,line=1.1)






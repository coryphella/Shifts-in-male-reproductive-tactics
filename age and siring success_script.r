######Model assessing the relationship of age and siring success

#inf_id: infant ID
#male_id: male ID
#mom_id: mom ID
#prtdobknown: male age known or estimated
#age_doc: male age at estimated time of conception
#sire: is male the sire (yes/no)
#PHG: group membership variables (ENK is reference)

####
d$male=as.integer(d$male_id)

d$s_age=(d$age_doc-mean(d$age_doc))/sd(d$age_doc)
d$agesq=d$s_age^2
library(rethinking)

#########'Null' model code


S0<- map2stan(  
alist(
sire ~ dbinom( 1 , p ) ,

	logit(p) <- a + a_mid[male] + b_group*PHG, 
			
    c(a,b_group) ~ dnorm(0,2),
	a_mid[male] ~ dnorm( 0,sigma_id ),
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2)
		),
data=d, cores=2 , chains=4, warmup=3500, iter=7000 , control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 

##########Model code
Sa<- map2stan(  
alist(
sire ~ dbinom( 1 , p ) ,

	logit(p) <- a  + b_age*s_age + b_agesq*agesq + b_group*PHG + a_mid[male] + b_age_m[male]*s_age + b_agesq_m[male]*agesq  ,  
			
    c(a,b_age,b_group,b_agesq) ~ dnorm(0,2),
	c(a_mid,b_age_m,b_agesq_m)[male] ~ dmvnormNC( sigma_mid, Rho_mid ),
	Rho_mid ~ dlkjcorr(3),
	sigma_mid ~ dcauchy(0,2)
		),
data=d, cores=2 , chains=2, warmup=3500, iter=7000, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 


output_Sa=precis(Sa, depth=2 , digits=2)@output
plot(precis(Sa, pars=c("a","b_rank","b_group","b_age","b_agesq"),depth=2))
plot(Sa)
compare(S0,Sa)

#########Constructing posterior predictions
a_male_z <- matrix(0,1000,length(unique(d$male)))
age.seq=seq(min(d$s_age),max(d$s_age),length=50)

pred_a1 <- list(
    male=rep(1,length(age.seq)),
	agesq=age.seq^2,
	s_age=age.seq,
	PHG=rep(mean(d$PHG),length(age.seq))	
)

Pa <- link(Sa, n=1000 , data=pred_a1, replace= 	list(a_mid=a_male_z,b_age_m=a_male_z,b_agesq_m=a_male_z), WAIC=TRUE)

pred.mean.a=apply(Pa , 2 , mean )
pred.HPDI.a=apply(Pa , 2 , HPDI )


##############Plotting
xxa=age.seq
yya=as.vector(t(Pa[1:1000,]))
par(mar=c(2.3,2.3,0.5,0.5))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_age),max(d$s_age)),colramp = colorRampPalette(c("white", "steelblue1")),nbin=200,transformation = function(x) x^.5,ylab='',xlab='',cex=1.2,xaxt='n',yaxt='n',ylim=c(0,1),nrpoints=0)
points(d$s_age[d$prtdobknown=="KNOWN"&d$sire==1],rep(1.027,length(d$s_age[d$prtdobknown=="KNOWN"&d$sire==1])), col=alpha("blue",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$prtdobknown=="EST"&d$sire==1],rep(1.027,length(d$s_age[d$prtdobknown=="EST"&d$sire==1])), col=alpha("darkturquoise",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$prtdobknown=="EST"&d$sire==0],rep(-0.027,length(d$s_age[d$prtdobknown=="EST"&d$sire==0])) , col=alpha("darkturquoise",0.8),pch=-as.hexmode("007C"),cex=0.25)
points(d$s_age[d$prtdobknown=="KNOWN"&d$sire==0],rep(-0.027,length(d$s_age[d$prtdobknown=="KNOWN"&d$sire==0])) , col=alpha("blue",0.8),pch=-as.hexmode("007C"),cex=0.25)

lines(age.seq,pred.mean.a ,lwd=1)
lines(age.seq,pred.HPDI.a[1,],lty=2,lwd=0.7)
lines(age.seq,pred.HPDI.a[2,],lty=2,lwd=0.7)

lab.a=seq(6,18,2)
lab.sa=(lab.a-mean(d$age_doc))/sd(d$age_doc)
axis(1,labels=NA,at=lab.sa,tck=-0.01)

mtext(lab.a,at=lab.sa,line=0.1,side=1,cex=0.7)

lab.s=pretty(c(0,1))

axis(2,labels=NA,at=lab.s,tck=-0.01)
mtext(lab.s,at=lab.s,line=0.2,side=2,cex=0.7)
mtext("Male age", side=1, line=1, cex=0.8,tck=-0.01)
mtext("Probability of siring", side=2, line=1.1, cex=0.8)

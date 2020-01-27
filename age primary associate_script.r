#####Models assesssing the relationship between age and primary association
#inf: infant ID
#male: male ID
#mom: mom ID
#prt1: indexing column for varying effects specification for non-directional dyadic variable
#prt2: indexing column for varying effects specification for non-directional dyadic variable
#prtdobknown: male age known or estimated
#prtage: male age
#primassoc: primary associate (yes/no)
#realdad: sire (yes/no)
#PHG, YNT: group membership variables (ENK is reference)

######################
library(rethinking)

#for indexing to work each ID has to be present at least once in each column prt1 and prt2
unique(d$mom[!(d$prt1 %in% d$prt2)]) #check whether each ID is present in each column
unique(d$male[!(d$prt2 %in% d$prt1)]) 

d$prt1=as.integer(droplevels(as.factor(d$prt1))) 
d$prt2=as.integer(droplevels(as.factor(d$prt2)))

d$dyad=paste(d$male,d$mom)
d$dyad=droplevels(as.factor(d$dyad))
d$dyad=as.integer(d$dyad)

d$s_age=(d$prtage-mean(d$prtage))/sd(d$prtage)

############'Null' model code
PA0<- map2stan(  
alist(
primassoc ~ dbinom( 1 , p ) ,

	logit(p) ~  a  + b_phg*PHG + b_ynt*YNT + a_id[prt2] + a_id[prt1] + sigma_dyad*dyad_id[dyad],
	
	
    c(a,b_phg,b_ynt) ~ dnorm(0,2),
	c(a_id)[prt1] ~ dnorm( 0, sigma_id ),
	dyad_id[dyad] ~ dnorm(0,sigma_dyad),
	sigma_id ~ dcauchy(0,2),
	sigma_dyad ~ dcauchy(0,2)
				),
data=d, cores=2, chains=4, warmup=3500, iter=7000, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 

##########Code for model including age
PA<- map2stan(  
alist(
primassoc ~ dbinom( 1 , p ) ,

	logit(p) ~  a + b_age*s_age  + b_phg*PHG + b_ynt*YNT + b_sire*realdad +
				a_id[prt2] + a_id[prt1] + sigma_dyad*dyad_id[dyad] + BA*s_age,
		
				BA ~ ba[prt1] + ba[prt2],
			
    c(a,b_age,b_phg,b_ynt,b_sire) ~ dnorm(0,2),
	c(a_id,ba)[prt1] ~ dmvnormNC( sigma_id, Rho_id ),
	dyad_id[dyad] ~ dnorm(0,10),
	c(Rho_id) ~ dlkjcorr(3),
	sigma_id ~ dcauchy(0,2),
	sigma_dyad ~ dcauchy(0,2)
			),
data=d, cores=2, chains=4, warmup=3500, iter=7000, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 


##########Code for model without age
PA_noage<- map2stan(  
alist(
primassoc ~ dbinom( 1 , p ) ,

	logit(p) ~  a    + b_phg*PHG + b_ynt*YNT + b_sire*realdad + a_id[prt2] + a_id[prt1] + sigma_dyad*dyad_id[dyad],

    c(a,b_phg,b_ynt,b_sire) ~ dnorm(0,2),
	c(a_id)[prt1] ~ dnorm( 0,sigma_id ),
	sigma_id ~ dcauchy(0,2),
	dyad_id[dyad] ~ dnorm(0,10),
	sigma_dyad ~ dcauchy(0,2)
				),
data=d, cores=2, chains=4, warmup=3500, iter=7000, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE) 

compare(PA,PA0,PA_noage)
res=precis(PA, depth=2 , digits=2)@output

###Constructing posterior predictions
a_prt_z <- matrix(0,1000,length(unique(d$prt1)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad)))

age.seq=seq(min(d$s_age),max(d$s_age),length=200)
age.seq2=seq(min(d$s_age),max(d$s_age),length=100)


pred_pa <- list(
    prt1=rep(1,length(age.seq)),
	prt2=rep(1,length(age.seq)),
	dyad=rep(1,length(age.seq)),
	realdad=rep(mean(d$realdad),length(age.seq)),
	s_age=age.seq,
	PHG=rep(mean(d$PHG),length(age.seq)),
	YNT=rep(mean(d$YNT),length(age.seq))
)	


L <- link(PA, n=1000 , data=pred_pa,replace= list(a_id=a_prt_z,ba=a_prt_z,dyad_id=a_dyad_z,bas=a_prt_z),  WAIC=TRUE)

pred.mean=apply(L$p , 2 , mean )
pred.HPDI=apply(L$p , 2 , HPDI )

#######Plotting
xxa=age.seq
yya=as.vector(t(Pa$p[1:1000,]))

par(mar=c(2.3,2.3,0.5,0.5))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_age),max(d$s_age)),colramp = colorRampPalette(c("white", "steelblue1")),nbin=200,transformation = function(x) x^.5,ylab='',xlab='',cex=1.2,xaxt='n',yaxt='n',ylim=c(0,1),nrpoints=0,bandwidth=0.01)
points(d$s_age[d$prtdobknown=="KNOWN"&d$primassoc==1],rep(1.027,length(d$s_age[d$prtdobknown=="KNOWN"&d$primassoc==1])), col=alpha("blue",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$prtdobknown=="EST"&d$primassoc==1],rep(1.027,length(d$s_age[d$prtdobknown=="EST"&d$primassoc==1])), col=alpha("darkturquoise",0.8),pch=-as.hexmode("007C") ,cex=0.25)
points(d$s_age[d$prtdobknown=="EST"&d$primassoc==0],rep(-0.027,length(d$s_age[d$prtdobknown=="EST"&d$primassoc==0])) , col=alpha("darkturquoise",0.8),pch=-as.hexmode("007C"),cex=0.25)
points(d$s_age[d$prtdobknown=="KNOWN"&d$primassoc==0],rep(-0.027,length(d$s_age[d$prtdobknown=="KNOWN"&d$primassoc==0])) , col=alpha("blue",0.8),pch=-as.hexmode("007C"),cex=0.25)

lines( age.seq , pred.mean,lwd=1)
lines(age.seq,pred.HPDI[1,],lty=2,lwd=0.7)
lines(age.seq,pred.HPDI[2,],lty=2,lwd=0.7)

lab.a=seq(6,16,2)
lab.sa=(lab.a-mean(d$prtage))/sd(d$prtage)
axis(1,labels=NA,at=lab.sa,tck=-0.01)
mtext(lab.a,at=lab.sa,line=0.1,side=1,cex=0.7)
lab.s=pretty(c(0,1))

axis(2,labels=NA,at=lab.s,tck=-0.01)
mtext(lab.s,at=lab.s,line=0.2,side=2,cex=0.7)
mtext("Male age", side=1, line=1, cex=0.8)
mtext("Probability of being the primary associate", side=2, line=1.1, cex=0.8)

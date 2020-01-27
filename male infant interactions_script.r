#Code for models assessing the relationship between primary association and male infant interactions 

#inf: infant ID
#male: male ID
#mom: mom ID
#prtdobknown: male age known or estimated
#prtage: male age
#primassoc: primary associate (yes/no)
#babyrate: male-infant interaction index
#PHG, YNT: group membership variables (ENK is reference)

############
library(rethinking)

d$male_index=as.integer(as.factor(d$male))
d$s_age=(d$prtage-mean(d$prtage))/sd(d$prtage)

#null model
null <- map2stan(
alist(

babyrate ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_phg*PHG + bp_ynt*YNT + ap_id[male_index],
				
  log(mu)  ~	am + bm_phg*PHG + bm_ynt*YNT + am_id[male_index],
  
    c(ap,am,bp_phg,bm_phg,bp_ynt,bm_ynt) ~ dnorm(0,2),
	c(ap_id,am_id)[male_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=3 , chains=4, warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

#########Code for model without age
inf_model <- map2stan(
alist(

babyrate ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_PA*primassoc +  bp_phg*PHG + bp_ynt*YNT + ap_id[male_index],
				
  log(mu)  ~	am + bm_PA*primassoc +  bm_phg*PHG + bm_ynt*YNT + am_id[male_index],
  
    c(ap,am,bm_PA,bp_PA,bp_phg,bm_phg,bp_ynt,bm_ynt) ~ dnorm(0,2),
	c(ap_id,am_id)[male_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=3 , chains=4, warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)


###########Code for model with age
inf_model_age <- map2stan(
alist(

babyrate ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_PA*primassoc + bp_age*s_age +  bp_phg*PHG + bp_ynt*YNT + bp_a[male_index]*s_age + ap_id[male_index],
				
  log(mu)  ~	am + bm_PA*primassoc + bm_age*s_age +  bm_phg*PHG + bm_ynt*YNT + bm_a[male_index]*s_age + am_id[male_index],
  
    c(ap,am,bm_PA,bp_PA,bp_age,bm_age,bp_phg,bm_phg,bp_ynt,bm_ynt) ~ dnorm(0,2),
	c(ap_id,am_id,bp_a,bm_a)[male_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2, warmup=3500, iter=7000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)


##############Constructing posterior predictions
a_male_z <- matrix(0,1000,length(unique(d$male)))

pred_PA <- list(
    male_index=1,
	primassoc=1,
	s_age=mean(d$s_age),
	YNT=mean(d$YNT),
	PHG=mean(d$PHG)	
)

pred_noPA <- list(
    male_index=1,
	primassoc=0,
	s_age=mean(d$s_age),
	YNT=mean(d$YNT),
	PHG=mean(d$PHG)	
)


L_PA <- link(inf_model_age, n=1000 , data=pred_PA, replace=	list(ap_id=a_male_z, am_id=a_male_z, bp_a=a_male_z,bm_a=a_male_z), WAIC=TRUE)
L_noPA <- link(inf_model_age, n=1000 , data=pred_noPA, replace=	list(ap_id=a_male_z, am_id=a_male_z, bp_a=a_male_z,bm_a=a_male_z), WAIC=TRUE)
L_PA_noage <- link(inf_model, n=1000 , data=pred_PA, replace=	list(ap_id=a_male_z, am_id=a_male_z, bp_a=a_male_z,bm_a=a_male_z), WAIC=TRUE)
L_noPA_noage <- link(inf_model, n=1000 , data=pred_noPA, replace=	list(ap_id=a_male_z, am_id=a_male_z, bp_a=a_male_z,bm_a=a_male_z), WAIC=TRUE)

P_PA <- (1-L_PA$p)*L_PA$mu
P_noPA <- (1-L_noPA$p)*L_noPA$mu
P_PA_noage <- (1-L_PA_noage$p)*L_PA_noage$mu
P_noPA_noage <- (1-L_noPA_noage$p)*L_noPA_noage$mu

######Averaging predictions
w <- compare(inf_model,inf_model_age.s,sort=FALSE)@output$weight  
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 ) ###round weights to nearest integer so samples we extract sum up to 1000


pred_PA <- P_PA 
pred_PA[1:idw[1],] <- P_PA_noage[1:idw[1],] 
pred_noPA <- P_noPA 
pred_noPA[1:idw[1],] <- P_noPA_noage[1:idw[1],] 

pred_PA <- sample(pred_PA)
median_PA=median(pred_PA)
HPDI_noage=HPDI(pred_PA)
pred_noPA <- sample(pred_noPA)
median_noPA=median(pred_noPA)
HPDI_noPA=HPDI(pred_noPA)


#####Plotting
#function 'dens' from 'rethinking' automatically sets and overrides margins, modified 'densm' code below to prevent that
densm(pred_PA, xlim=c(0,1.8) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-1.5,36),mar=c(0,0,0,0))
axis(1, tck=-0.02,cex=1,labels=NA,at=seq(0,1.8,by=0.1))
axis(1, cex.axis=0.8,at= seq(0,1.8,by=0.2),labels=seq(0,1.8,by=0.2),line=-0.8,col=NA)
ll <- d$babyrate[d$primassoc==1]
points(ll, rep(-2.2,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_PA) , lim= as.vector(HPDI(pred_PA, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_noPA) , lim= as.vector(HPDI(pred_noPA, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$babyrate[d$primassoc==0]
points(ll, rep(-0.8,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(P_PA) , lty=1,lwd=1)
abline(v=median(P_noPA) , lty=2,lwd=1)

mtext(side=1,line=1.3,cex=0.9,text="Composite index of male infant relationship")

legend(x=1.2,y=35, legend = c("Primary associates","Other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("orange", 0.5) ) , pch=c(15,16),
       pt.cex=2 , bty="n",y.intersp=1.3,x.intersp=1, lty=c(0,0) ,cex=0.9,lwd=c(1,1))
legend(x=1.2,y=35, legend = c("Primary associates","Other males"),
       col=1 , pch=c(NA,NA),
       pt.cex=2 , bty="n",y.intersp=1.3,x.intersp=1, lty=c(1,2) ,cex=0.9,lwd=c(1,1))




#modified 'dens' function
densm=function (x, adj = 0.5, norm.comp = FALSE, main = "", show.HPDI = FALSE, 
    show.zero = FALSE, rm.na = TRUE, add = FALSE, ...) 
{
    the.class <- class(x)[1]
    if (the.class == "data.frame") {
        n <- ncol(x)
        cnames <- colnames(x)
        set_nice_margins()
        par(mfrow = make.grid(n))
        for (i in 1:n) {
            dens(x[, i], adj = adj, norm.comp = norm.comp, show.HPDI = show.HPDI, 
                show.zero = TRUE, xlab = cnames[i], ...)
        }
    }
    else {
        if (rm.na == TRUE) 
            x <- x[!is.na(x)]
        thed <- density(x, adjust = adj)
        if (add == FALSE) {
            plot(thed, main = main, ...)
        }
        else lines(thed$x, thed$y, ...)
        if (show.HPDI != FALSE) {
            hpd <- HPDI(x, prob = show.HPDI)
            shade(thed, hpd)
        }
        if (norm.comp == TRUE) {
            mu <- mean(x)
            sigma <- sd(x)
            curve(dnorm(x, mu, sigma), col = "white", lwd = 2, 
                add = TRUE)
            curve(dnorm(x, mu, sigma), add = TRUE)
        }
        if (show.zero == TRUE) {
            lines(c(0, 0), c(0, max(thed$y) * 2), lty = 2)
        }
    }
}
<bytecode: 0x000000002862bf40>
<environment: namespace:rethinking>
> set_nice_margins
function () 
{
    par_mf <- par("mfrow", "mfcol")
    if (all(unlist(par_mf) == 1)) {
        par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
            tck = -0.02)
    }
}

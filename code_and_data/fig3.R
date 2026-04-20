rm(list=ls(all=TRUE))
###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

###################################################
## Install package
###################################################
{
library(dplyr)
library(tidyr)
library(FME)
library(ggplot2)
library(GillespieSSA)
library(patchwork)
}
###################################################
## Basic configuration
###################################################

pltsetb <- theme( plot.background = element_blank() ) +
theme( axis.text = element_text(colour = "black", size = 20)) +
theme( axis.ticks = element_line(colour = "black", size = 1) ) +
theme( axis.ticks.length = unit(0.5,"cm") ) +
theme( axis.line.x = element_line(colour = "black", size = 1) ) +
theme( axis.line.y = element_line(colour = "black", size = 1) ) +
theme( axis.title = element_text(colour = "black", size = 20)) +
theme( panel.background=element_blank() ) +
theme( panel.grid.major=element_blank() ) +
theme( panel.grid.minor=element_blank() )

caka <- "#CC3366"; cao <- "#3366CC";
cakao <- "#CC336650"; caoo <- "#3366CC50";
cmidori <- "#009966"; cmidorio <- "#00996650";
orange <- "#ff8c00"
orangekage <- "#ff8c0050"
midori <- "#008000"
midorikage <- "#00800050"

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###################################################
## Data reading
###################################################
load("final.Rdata")
print(fit_mixture[c(1,4,6:10)])

fitpar <- fit_mixture[c(1,4,6:10)]
load("kakunin/MCMC_mixture.Rdata")
parset <- data.frame(MCMC_mixture$pars)
av_par <- parset %>%
  summarise(x0=mean(x0),
            dx=mean(dx),
            lambda=mean(lambda),
            f=mean(f),
            alpha=mean(alpha),
            rho=mean(rho),
            rr=mean(rr)) %>% as.numeric()
mode_par <- parset %>%
  summarise(x0=Mode(x0),
            dx=Mode(dx),
            lambda=Mode(lambda),
            f=Mode(f),
            alpha=Mode(alpha),
            rho=Mode(rho),
            rr=Mode(rr)) %>% as.numeric()
###################################################
## SSR Setting
###################################################

Tmin <- 0.0
Tmax <- 100.0
step_size <- 1

stime <- seq(Tmin,Tmax,step_size)

## fatal ##
system("R CMD SHLIB simulate.c")
dyn.load( paste("simulate",.Platform$dynlib.ext,sep="") )

## ODE model ##
ODEs <- function(pars) {

    #initial value is squared parameters
    rhs <- c(x=as.numeric(pars[1]),y=0)
    times <- seq(Tmin,Tmax,step_size)
    pars <- pars
    out <- ode(y=rhs,parms=pars,times=times,func="derivs",initfunc="initparms",nout=1,outnames=c(""),dllname="simulate",method="rk4")
    as.data.frame(out)

}

sensrfun <- function(pars) {
  return( ODEs(pars) )
}

###
out <- ODEs(mode_par)%>%
  filter(time<=50)

### auc_integrated
integrated_hbv <- function(lambda,r,t){
  val <- lambda*(1-exp(-r*t))
  return(val)
}

mode_par

lambda <- mode_par[3]
r <- mode_par[7]
int.df <- data.frame(time=seq(0,50),
                     integrate=integrated_hbv(lambda,r,seq(0,50)))
rep.df <- left_join(out,int.df,by="time") %>%
  mutate(intra=mode_par[5]*x) %>%
  mutate(rate.integ=integrate/(integrate+intra)) %>%
  filter(time>0)

max_integrate <- max(rep.df$integrate)
scale_factor <- 100 / log10(max_integrate + 1)
plt <- ggplot(rep.df, aes(x=time)) +
  geom_bar(aes(y=rate.integ * 100), stat="identity", fill="grey", alpha=.25) +
  geom_line(aes(y=log10(integrate + 1) * scale_factor), colour=caka,lwd=2) +
  geom_line(aes(y=log10(intra + 1) * scale_factor), colour=cao,lwd=2) +
  scale_y_continuous(
    name="Contribution of integrated HBV DNA (%)",
    sec.axis = sec_axis(~(. / scale_factor),
                        name="Intracellular HBV DNA originated from\n cccDNA / integrated HBV DNA (copies/well)",
                        breaks = seq(0,8,2),labels = c(expression(10^0,10^2,10^4,10^6,10^8))))+
  xlab("Days")+
  pltsetb
summary(rep.df)
sapply(rep.df, function(x) sum(!is.finite(x)))

# ggsave("fig/contribution_hbv2.png",plt,w=11,h=6)

fitpar
#積分
auc_result <- integrate(integrated_hbv, lower = 0, upper = 30, lambda = lambda, r = r)
ccchbv <- (sum(out$x)*as.numeric(fitpar[7]))/30
integ <- ((auc_result$value)/30)
sprintf("%.2e",ccchbv)
sprintf("%.2e",integ)
sprintf("%.2e",integ/(integ+ccchbv))

###MCMC ver
#cccDNA
sR <- sensRange(func=sensrfun,parms=NULL,parInput=MCMC_mixture$pars[,c(1,4,6:10)])
name_lst <- colnames(sR)
gids <- grep("x+[[:digit:]]+",name_lst)
mat <- sR[,gids]
mat <- mat[,-1]
alpha <- MCMC_mixture$pars[,8]
#alphaの値をかける
for (i in 1:length(mat)) {
  mat[, i] <- mat[, i] * alpha[i]
}
times <- seq(Tmin,Tmax,step_size)
cmean <- log10( apply(mat,2,mean) )
yCIlow <- log10( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- log10( apply(mat,2,function(x){quantile(x,0.975)}) )
clow <- as.numeric(yCIlow)
cup  <- as.numeric(yCIhigh)

clow[clow<1.0] <- 1.0

xrange <- c(Tmin,Tmax)
yrange <- c(min(clow),max(cup))
xn <- length(times)
yn <- length(cmean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
x <- rep(times,3)
c.df <- data.frame(x=times,mean=cmean,ymin=clow,ymax=cup)

#integrate DNA
name_lst <- colnames(sR)
gids <- grep("V4+[[:digit:]]+",name_lst)
i.mat <- sR[,gids]

imean <- log10( apply(i.mat,2,mean) )
yCIlow <- log10( apply(i.mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- log10( apply(i.mat,2,function(x){quantile(x,0.975)}) )
ilow <- as.numeric(yCIlow)
iup  <- as.numeric(yCIhigh)

ilow[ilow<1.0] <- 1.0

xrange <- c(Tmin,Tmax)
yrange <- c(min(ilow),max(iup))
xn <- length(times)
yn <- length(imean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
x <- rep(times,3)
i.df <- data.frame(x=times,mean=imean,ymin=ilow,ymax=iup) %>%
  filter(x>0)

test.mat <- data.frame(time=times)
for (i in 1:nrow(mat)) {
  test.mat <- as.numeric(mat[i,]/(mat[i,]+i.mat[i,])) %>%
    as.data.frame() %>%
    cbind(test.mat,.)
  kariname <- paste0("test",i)
  colnames(test.mat)[i+1] <- kariname
}

summary.df <- data.frame()
for (i in 2:nrow(test.mat)) {
  mean.val <- mean(as.numeric(test.mat[i,-1]))
  sd.val <- sd(test.mat[i,-1])
  karidf <- data.frame(time=times[i],
                       mean=mean.val,
                       sd=sd.val)
  summary.df <- rbind(summary.df,karidf)
}

max_integrate <- max(i.df$mean)
scale_factor <- 100 / (max_integrate)
summary.df$mean <- summary.df$mean*100
summary.df$sd <- summary.df$sd*100
summary.df <- summary.df %>% filter(time<=50)

summary.df

plt <- ggplot() +
  geom_bar(data=summary.df,aes(x=time,y=mean), stat="identity", fill="gray", alpha=.75) +
  geom_errorbar(data=summary.df,aes(x=time, ymax=mean+(0.434*sd),ymin=mean-(0.434*sd) ),
                color="black", size=0.7, width=0.8)+
  geom_ribbon(data=i.df,aes(x=x,ymin=ymin* scale_factor,ymax=ymax* scale_factor),fill=orangekage) +
  geom_line(data=i.df,aes(x=x,y=mean * scale_factor), colour=orange,lwd=2) +
  geom_ribbon(data=c.df,aes(x=x,ymin=ymin* scale_factor,ymax=ymax* scale_factor),fill=midorikage) +
  geom_line(data=c.df,aes(x=x,y=mean * scale_factor), colour=midori,lwd=2) +
  scale_y_continuous(
    name="Contribution of cccDNA (%)",
    sec.axis = sec_axis(~(. / scale_factor),
                        name="Intracellular HBV DNA originated from\n cccDNA / integrated HBV DNA (copies/well)",
                        breaks = seq(0,8,2),labels = c(expression(10^0,10^2,10^4,10^6,10^8))))+
  scale_x_continuous(limits = c(0,51))+
  xlab("Days post infection")+
  pltsetb
plt

ggsave("fig/fig3.png",plt,w=9.75,h=6.5)


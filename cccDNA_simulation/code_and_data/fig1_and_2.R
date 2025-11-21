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
  library(tidyverse)
  library(patchwork)
  library(reshape2)
  library(FME)
  library(ggplot2)
  library(scales)
  library(ggpmisc)
  library(DescTools)
  library(HDInterval)
  library(tidyr)
}

###################################################
## Basic configuration
###################################################
{
mcmc_iter <- 200000
Tmin <- 0.0
Tmax <- 60.0
Tmnew <- 60.0
step_size <- 0.1
stime <- seq(Tmin,Tmax,step_size)
lbt <- 0.001 # lower bound threshold
ubt <- 1000.0 # upper bound threshold
pars <- c(x0=1.829e5,x40=8.464e4,y0 = 1.0,dx=0.193,epsilon=0.999,lambda=1.115e14,f=0.0877,alpha=1,rho=0.1,rr=1.4e-7,rr1=2.8e-7)

datatypes <- c("condition1_ccc","condition1_dna",
               "condition2_ccc","condition2_dna",
               "condition3_ccc","condition3_dna",
               "condition4_ccc","condition4_dna")
datakakunin <- c("ccc4","dna4","ccc5","dna5","ccc6","dna6","ccc7","dna7")
erbc <- c("datc","datc1","datc2","datc3","datc_ctr2","datc4","datc5","datc6","datc7")
erbd <- c("datd","datd1","datd2","datd3","datd_ctr2","datd4","datd5","datd6","datd7")
erbe <- c("date","date1","date2","date3","date_ctr2","date4","date5","date6","date7")

tme30 <- 30
tme45 <- 45
tme9 <- 9
tme18 <- 18

pltsetb <- theme(legend.position="none") +
  theme( plot.background = element_blank() ) +
  theme( text=element_text("Helvetica") ) +
  theme( axis.text = element_text(colour = "black", size = 15)) +
  theme( axis.ticks = element_line(colour = "black", size = 1) ) +
  theme( axis.ticks.length = unit(0.2,"cm") ) +
  theme( axis.line.x = element_line(colour = "black", size = 1) ) +
  theme( axis.line.y = element_line(colour = "black", size = 1) ) +
  theme( axis.title = element_text(colour = "black", size = 15) ) +
  theme( panel.background=element_blank() ) +
  theme( panel.grid.major=element_blank() ) +
  theme( panel.grid.minor=element_blank() )

caka <- "#CC3366"; cakao <- "#CC336650";
cao <- "#3366CC"; caoo <- "#3366CC50";
cmidori <- "#009966"; cmidorio <- "#00996650";
corange <- "#FF6600"; corangeo <- "#FFCC33";
cmurasaki <- "#6600FF"; cmurasakio <- "#CC66FF"
cpink <- "deeppink3"; cpinko <- "deeppink1"
ccha <- "#BB4513"; cchao <- "#CD661D"
paka <- "#ff7f7f"; pakao <- "#ff9e9e"
pao <- "#8484ff"; paoo <- "#a3a3ff"
pmidori <- "#7fff7f"; pmidorio <- "#9eff9e"
}


obs_data <- vector(mode="list",length=length(datatypes))
for( i in 1:length(datatypes) ) {
    rfn <- paste("data/",datatypes[i],".txt",sep="")
    dsf <- read.csv(rfn,sep="\t",header=F)
    colnames(dsf) <- c("time","value")
    dtime <- dsf[,1]
    tn <- length(dtime)
    ids <- rep(0,length=tn)
    for( j in 1:tn ) {
        ids[j] <- which( abs(stime-dtime[j])==min(abs(stime-dtime[j])) )
    }
    obs_data[[i]] <- cbind(ids,dsf[,2])
}


###################################################
## functions
###################################################
{
## fatal ##
  if ("final" %in% names(getLoadedDLLs())) {
    dyn.unload(paste0("final", .Platform$dynlib.ext))
  }
  suppressWarnings(file.remove("final.o", paste0("final", .Platform$dynlib.ext)))  # 古いビルドを消す
  system("R CMD SHLIB final.c")                        # 再コンパイル
  dyn.load(paste0("final", .Platform$dynlib.ext))
## ODE model ##
ODEs <- function(pars) {

    # initial conditions
    x00 <- as.numeric(pars[1])
    y00 <- as.numeric(pars[3])
    x10 <- x00
    y10 <- y00
    x20 <- x00
    y20 <- y00
    x30 <- x00
    y30 <- y00
    x40 <- as.numeric(pars[2])
    y40 <- 0.0
    rhs <- c(x=x00,y=y00,x1s=x00,y1s=y00,x2s=x00,y2s=y00,x3s=x00,y3s=y00,x4s=x40,y4s=y40)
    times <- seq(Tmin,Tmax,step_size)
    ## C compiled version ##
    pars <- c(pars,tme30,tme45,tme9,tme18)
    out <- ode(y=rhs,parms=pars,times=times,func="derivs",initfunc="initparms",nout=1,outnames=c(""),dllname="final")
    as.data.frame(out)
}

ODEk <- function(pars) {

    # initial conditions
    x00 <- as.numeric(pars[1])
    y00 <- as.numeric(pars[3])
    x10 <- x00
    y10 <- y00
    x20 <- x00
    y20 <- y00
    x30 <- x00
    y30 <- y00
    x40 <- as.numeric(pars[2])
    y40 <- y00
    rhs <- c(x=x00,y=y00,x1s=x00,y1s=y00,x2s=x00,y2s=y00,x3s=x30,y3s=y30,x4s=x40,y4s=y40)
    times <- seq(Tmin,Tmax,step_size)
    ## C compiled version ##
    pars <- c(pars,tme30,tme45,tme9,tme18)
    out <- ode(y=rhs,parms=pars,times=times,func="derivs1",initfunc="initparms",nout=1,outnames=c(""),dllname="final")
    as.data.frame(out)
}

SSR <- function(pars) {
  val <- 0.0
  out <- try(ODEs(pars), silent = TRUE)
  if (inherits(out, "try-error")) return(Inf)
  if (any(!is.finite(as.matrix(out)))) return(Inf)

  for (i in 1:length(datatypes)) {
    ids <- obs_data[[i]][,1]
    datavalues <- obs_data[[i]][,2]
    for (j in 1:length(ids)) {
      simval <- pmax(out[ids[j], i+1], 1e-8)
      datval <- pmax(datavalues[j], 1e-8)
      val <- val + (log(simval) - log(datval))^2
    }
  }
  return(val)
}

## log-SSR ##
SSRlog <- function(lpars) {
    SSR(exp(lpars))
}

set_constr <- function(low_vec,up_vec) {
    n_vec <- length(low_vec)
    stopifnot( length(up_vec)==n_vec )
    ## for constrOptim(): we can use the same rule to interexchange components so that we do not have to do it.##
    AL <- diag(n_vec)
    AU <- -diag(n_vec)
    Amat <- rbind(AL,AU) # 'ui'
    bvec <- c(low_vec,-up_vec) # 'ci'
    return(list(low.val=low_vec,up.val=up_vec,Amat=Amat,bvec=bvec))
}

## sensitivity function ##
sensrfun <- function(pars) {
    return( ODEs(pars) )
}

sensrfun1 <- function(pars) {
    return ( ODEk(pars) )
}
}
###################################################
## Parameter fitting
###################################################

{
  lbtvec <- lbt * pars
  ubtvec <- ubt * pars
  lbtvec["y0"] <- 1e-3
  ubtvec["y0"] <- 1e2
  loglbtvec <- log(lbtvec)
  logubtvec <- log(ubtvec)
  cbox <- set_constr(loglbtvec, logubtvec)

  safe_SSR <- function(pars) {
    out <- try(ODEs(pars), silent = TRUE)
    if (inherits(out, "try-error")) return(Inf)
    if (any(!is.finite(as.matrix(out)))) return(Inf)
    val <- 0.0
    for (i in 1:length(datatypes)) {
      ids <- obs_data[[i]][,1]
      datavalues <- obs_data[[i]][,2]
      for (j in 1:length(ids)) {
        simval <- pmax(out[ids[j], i+1], 1e-6)
        datval <- pmax(datavalues[j], 1e-6)
        val <- val + (log(simval) - log(datval))^2
      }
    }
    return(val)
  }

  safe_SSRlog <- function(lpars) safe_SSR(exp(lpars))

  theta0 <- log(pars)
  theta0 <- pmin(pmax(theta0, loglbtvec + 1e-6), logubtvec - 1e-6)

  est <- constrOptim(
    theta = theta0,
    f = safe_SSRlog,
    grad = NULL,
    method = "Nelder-Mead",
    ui = cbox$Amat,
    ci = cbox$bvec,
    control = list(trace = 1, maxit = 500, reltol = 1e-16)
  )

  fit_mixture <- exp(est$par)
  print(pars)
  print(safe_SSR(pars))
  print(fit_mixture)
  print(safe_SSR(fit_mixture))

  save(fit_mixture, file = "final.Rdata")
}

###################################################
## Plot fitted curve
###################################################

fitted <- ODEs(pars=fit_mixture)

data_list <- vector(mode="list",length=length(datatypes))
for( i in 1:length(datatypes) ) {
    ## experimental data ##
    rfn <- paste("data/",datatypes[i],".txt",sep="")
    dsf <- read.csv(rfn,sep="\t",header=F)
    data_list[[i]] <- data.frame( time=dsf[,1],value=log10(dsf[,2]) )
    ## simulation adjustment ##
    zeroid <- which (fitted[,i]==0.0 )
    fitted[zeroid,i] <- 0.1
}

## control ##
dd_CCC <- data_list[[1]]
dd_DNA <- data_list[[2]]
ds_CCC <- data.frame(time=fitted$time,value=log10(fitted[,2]),type=rep(datatypes[1],length=nrow(fitted)))
ds_DNA <- data.frame(time=fitted$time,value=log10(fitted[,3]),type=rep(datatypes[2],length=nrow(fitted)))
## ETV0 ##
dd_CCC1 <- data_list[[3]]
dd_DNA1 <- data_list[[4]]
ds_CCC1 <- data.frame(time=fitted$time,value=log10(fitted[,4]),type=rep(datatypes[3],length=nrow(fitted)))
ds_DNA1 <- data.frame(time=fitted$time,value=log10(fitted[,5]),type=rep(datatypes[4],length=nrow(fitted)))
#ds_eDNA <- data.frame(time=fitted$time,value=log(fitted[,4]),type=rep(datatypes[4],length=nrow(fitted)))
## ETV30 ##
dd_CCC2 <- data_list[[5]]
dd_DNA2 <- data_list[[6]]
ds_CCC2 <- data.frame(time=fitted$time,value=log10(fitted[,6]),type=rep(datatypes[5],length=nrow(fitted)))
ds_DNA2 <- data.frame(time=fitted$time,value=log10(fitted[,7]),type=rep(datatypes[6],length=nrow(fitted)))
## TET30 ##
dd_CCC3 <- data_list[[7]]
dd_DNA3 <- data_list[[8]]
ds_CCC3 <- data.frame(time=fitted$time,value=log10(fitted[,8]),type=rep(datatypes[7],length=nrow(fitted)))
ds_DNA3 <- data.frame(time=fitted$time,value=log10(fitted[,9]),type=rep(datatypes[8],length=nrow(fitted)))

## plotfinc ##
fig_plot_func <- function(sccc,sdna,dccc,ddna){
  plt <- ggplot() +
    geom_path(data=sccc,aes(x=time,y=value),color=cao,lwd=1) +
    geom_path(data=sdna,aes(x=time,y=value),color=caka,lwd=1) +
    geom_point(data=dccc,aes(x=time,y=value),color=cao,size=6) +
    geom_point(data=ddna,aes(x=time,y=value),color=caka,size=6) +
    labs(x="days",y="ccc/HBV DNA (copies/well)")+
    scale_x_continuous(breaks=seq(0,60,by=10),labels=c("0","10","20","30","40","50","60")) +
    scale_y_continuous(limits=c(0,12)) + pltsetb
  return(plt)
}


## plot control ##
plt <- fig_plot_func(ds_CCC,ds_DNA,dd_CCC,dd_DNA)
plt

# plot ETV t=0 ##
plt <- fig_plot_func(ds_CCC1,ds_DNA1,dd_CCC1,dd_DNA1)
plt

# plot ETV t=30 ##
plt <- fig_plot_func(ds_CCC2,ds_DNA2,dd_CCC2,dd_DNA2)
plt

# plot ETV t=0 ##
plt <- fig_plot_func(ds_CCC3,ds_DNA3,dd_CCC3,dd_DNA3)
plt


##################################################
# MCMC
##################################################

# {
#   mcmc_pars <- fit_mixture
#   var0 <- NULL
#   cov0 <- 0.001 * diag(length(pars))
#   MCMC_mixture <- modMCMC(f=SSRlog,p=log(mcmc_pars),lower=log(lbt*(mcmc_pars)),upper=log(ubt*(mcmc_pars)),niter=mcmc_iter,jump=cov0,var0=var0,wvar0=0.1,updatecov=20)
#   MCMC_mixture$pars <- exp(MCMC_mixture$pars)
#   save(MCMC_mixture,file="kakunin/MCMC_mixture.Rdata")
# }

{
  load("kakunin/MCMC_mixture.Rdata")
  png("kakunin/MCMC_mixture.png",width=600,height=600)
  plot(MCMC_mixture,Full=TRUE)
  pairs(MCMC_mixture,nsample=1000)
  dev.off()
}

{
  library(coda)
  load("kakunin/MCMC_mixture.Rdata")
  cobj <- as.mcmc(MCMC_mixture$pars)
  pdf("kakunin/CMC_mixture_conv.pdf")
  plot(cobj)
  dev.off()
}

###################################################
## Sensitivity ranges
###################################################

{
    load("kakunin/MCMC_mixture.Rdata")

    # fatal #
    mx0 <- MCMC_mixture$pars[,1]
    mx40 <- MCMC_mixture$pars[,2]
    mdx <- MCMC_mixture$pars[,3]
    mep <- MCMC_mixture$pars[,4]
    mlam <- MCMC_mixture$pars[,5]
    mf <- MCMC_mixture$pars[,6]
    mal <- MCMC_mixture$pars[,7]
    mrho <- MCMC_mixture$pars[,8]
    mrr <- MCMC_mixture$pars[,9]
    mrr1 <- MCMC_mixture$pars[,10]

    sR <- sensRange(func=sensrfun,parms=NULL,parInput=MCMC_mixture$pars)
    save(sR,file="kakunin/MCMC_globsens.Rdata")
    png("MCMC_globsens.png",width=600,height=600)
    pdf("kakunin/MCMC_globsens.pdf")
    plot(summary(sR),xlab="time")
    dev.off()
}


  ## data reading ##
  data_list <- vector(mode="list",length=length(datatypes))
  for( i in 1:length(datatypes) ) {
      ## experimental data ##
      rfn <- paste("data/",datatypes[i],".txt",sep="")
      dsf <- read.csv(rfn,sep="\t",header=F)
      data_list[[i]] <- data.frame( time=dsf[,1],value=log10(dsf[,2]) )
  }

  # data_list
  for (i in 1:length(datatypes)) {
      if (1==i%%2) {
          assign(sprintf("ccc%d",(i-1)/2),data_list[[i]])
      } else if (0==i%%2) {
          assign(sprintf("dna%d",(i-2)/2),data_list[[i]])
      }
  }

  ## load global sensitivity analysis result ##
  load("kakunin/MCMC_globsens.Rdata")

  dataframeC.make <- function(mat){
    c.list <- list()

    cmean <- log10( apply(mat,2,mean) )
    yCIlow <- log10( apply(mat,2,function(x){quantile(x,0.025)}) )
    yCIhigh <- log10( apply(mat,2,function(x){quantile(x,0.975)}) )
    clow <- as.numeric(yCIlow)
    cup  <- as.numeric(yCIhigh)

    # set 0 for un-detectable virusload #
    clow[clow<1.0] <- 1.0

    xrange <- c(Tmin,Tmax)
    yrange <- c(min(clow),max(cup))
    xn <- length(times)
    yn <- length(cmean)

    labels <- gl(3,xn,label=c("mean","5%","95%"))
    x <- rep(times,3)
    c.list[[1]] <- data.frame(x=times,mean=cmean,ymin=clow,ymax=cup)

    c.list[[2]] <- data.frame(x=times,y=clow)
    c.list[[3]] <- data.frame(x=times,y=cmean)
    c.list[[4]] <- data.frame(x=times,y=cup)

    return(c.list)
  }

  dataframeD.make <- function(mat){
    d.list <- list()

    dmean <- log10( apply(mat,2,mean) )
    yCIlow <- log10( apply(mat,2,function(x){quantile(x,0.025)}) )
    yCIhigh <- log10( apply(mat,2,function(x){quantile(x,0.975)}) )
    dlow <- as.numeric(yCIlow)
    dup  <- as.numeric(yCIhigh)

    # set 0 for un-detectable virusload #
    dlow[dlow<1.0] <- 1.0

    xrange <- c(Tmin,Tmax)
    yrange <- c(min(dlow),max(dup))
    xn <- length(times)
    yn <- length(dmean)

    labels <- gl(3,xn,label=c("mean","2.5%","97.5%"))
    x <- rep(times,3)
    d.list[[1]] <- data.frame(x=times,mean=dmean,ymin=dlow,ymax=dup)

    d.list[[2]] <- data.frame(x=times,y=dlow)
    d.list[[3]] <- data.frame(x=times,y=dmean)
    d.list[[4]] <- data.frame(x=times,y=dup)

    return(d.list)
  }

  mcmc_plot_func <- function(c.list,d.list,ccc,dna){
    plt <- ggplot() +
      geom_ribbon(data=c.list[[1]],aes(x=x,ymin=ymin,ymax=ymax),fill=caoo) +
      geom_line(data=c.list[[2]],aes(x=x,y=y),lwd=1,color=caoo) +
      geom_line(data=c.list[[3]],aes(x=x,y=y),lwd=1,color=cao) +
      geom_line(data=c.list[[4]],aes(x=x,y=y),lwd=1,color=caoo) +
      geom_ribbon(data=d.list[[1]],aes(x=x,ymin=ymin,ymax=ymax),fill=cakao) +
      geom_line(data=d.list[[2]],aes(x=x,y=y),lwd=1,color=cakao) +
      geom_line(data=d.list[[3]],aes(x=x,y=y),lwd=1,color=caka) +
      geom_line(data=d.list[[4]],aes(x=x,y=y),lwd=1,color=cakao) +
      geom_point(data=ccc,aes(x=time,y=value),colour=cao,size=4) +
      geom_point(data=dna,aes(x=time,y=value),colour=caka,size=4) +
      labs(x = "Days post infection",y = "ccc / HBV DNA (copies/well)")+
      scale_x_continuous(breaks=seq(0,60,by=10),labels=c("0","10","20","30","40","50","60")) +
      scale_y_continuous(limits = c(0,12),breaks = seq(0,12,2),labels = c(expression(10^0,10^2,10^4,10^6,10^8,10^10,10^12))) +
      pltsetb

    return(plt)
  }

  Tmin_plot <- 0
  Tmax_plot <- Tmax

  plot_list <- list()

  load("final.Rdata")

  parset <- as.data.frame(MCMC_mixture$pars)


  ### condition1
  name_lst <- colnames(sR)
  gids <- grep("x[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  c.list <- dataframeC.make(mat)

  name_lst <- colnames(sR)
  gids <- grep("y[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))

  d.list <- dataframeD.make(mat)

  d.list


  check.Cdf0 <- c.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,ccc0,by="time") %>% na.omit()
  check.Ddf0 <- d.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,dna0,by="time") %>% na.omit()


  ## plot ##
  plt <- mcmc_plot_func(c.list,d.list,ccc0,dna0)+
    scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
  plot_list[[1]] <- plt
  ggsave("fig/condition1.png",plt,w=5,h=5)

  ### condition2
  name_lst <- colnames(sR)
  gids <- grep("x1s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  c.list <- dataframeC.make(mat)

  name_lst <- colnames(sR)
  gids <- grep("y1s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  d.list <- dataframeD.make(mat)

  check.Cdf1 <- c.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,ccc1,by="time") %>% na.omit()
  check.Ddf1 <- d.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,dna1,by="time") %>% na.omit()

  ## plot ##
  plt <- mcmc_plot_func(c.list,d.list,ccc1,dna1)+
    geom_vline(xintercept = 0,color=cmurasaki,lwd=2,alpha=.5)+
    scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))

  plot_list[[2]] <- plt
  ggsave("fig/condition2.png",plt,w=5,h=5)

  ### condition3
  name_lst <- colnames(sR)
  gids <- grep("x2s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  c.list <- dataframeC.make(mat)

  name_lst <- colnames(sR)
  gids <- grep("y2s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  d.list <- dataframeD.make(mat)

  check.Cdf2 <- c.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,ccc2,by="time") %>% na.omit()
  check.Ddf2 <- d.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,dna2,by="time") %>% na.omit()

  ## plot ##
  plt <- mcmc_plot_func(c.list,d.list,ccc2,dna2)+
    geom_vline(xintercept = 30,color=cmurasaki,lwd=2,alpha=.5)+
    scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
  plot_list[[3]] <- plt
  ggsave("fig/condition3.png",plt,w=5,h=5)

  ### condition4
  name_lst <- colnames(sR)
  gids <- grep("x3s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  c.list <- dataframeC.make(mat)

  name_lst <- colnames(sR)
  gids <- grep("y3s[[:digit:]]+\\.[[:digit:]]+",name_lst)
  mat <- sR[,gids]
  times <- seq(Tmin,Tmax,length=length(gids))
  d.list <- dataframeD.make(mat)

  check.Cdf3 <- c.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,ccc3,by="time") %>% na.omit()
  check.Ddf3 <- d.list[[3]] %>%
    mutate(x=round(x,1)) %>%
    rename("time"="x") %>%
    left_join(.,dna3,by="time") %>% na.omit()

  ## plot ##
  plt <- mcmc_plot_func(c.list,d.list,ccc3,dna3)+
    geom_vline(xintercept = 30,color="black",lwd=2,alpha=.5)+
    scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
  plot_list[[4]] <- plt
  ggsave("fig/condition4.png",plt,w=5,h=5)

#################################################
# fig1/2
#################################################

{
  load("kakunin/MCMC_mixture.Rdata")

  sR <- sensRange(func=sensrfun1,parms=NULL,parInput=MCMC_mixture$pars)
  save(sR,file="kakunin/MCMC_globsens_kakunin2.Rdata")
  pdf("kakunin/MCMC_globsens_kakunin.pdf")
  plot(summary(sR),xlab="time")
  dev.off()
}

erb_c <- vector(mode="list",length=length(erbc))
erb_d <- vector(mode="list",length=length(erbd))
erb_e <- vector(mode="list",length=length(erbe))
for (i in 1:length(erb_c)) {
  rfnc <- paste("plot/",erbc[i],".csv",sep="")
  dsfc <- read.csv(rfnc, sep=",", header=F)
  # grpc <- rep(erbc[i],nrow(dsfc))
  erb_c[[i]] <- data.frame( time=dsfc[,1], value=dsfc[,2], sd=dsfc[,3])

  rfnd <- paste("plot/",erbd[i],".csv",sep="")
  dsfd <- read.csv(rfnd, sep=",", header=F)
  grpd <- rep(erbd[i],nrow(dsfd))
  erb_d[[i]] <- data.frame( time=dsfd[,1], value=dsfd[,2], sd=dsfd[,3])

  rfne <- paste("plot/",erbe[i],".csv",sep="")
  dsfe <- read.csv(rfne, sep=",", header=F)
  grpe <- rep(erbe[i],nrow(dsfe))
  erb_e[[i]] <- data.frame( time=dsfe[,1], value=dsfe[,2], sd=dsfe[,3])
}

load("kakunin/MCMC_globsens_kakunin.Rdata")

kakunin_list <- vector(mode="list",length=length(datakakunin))
for( i in 1:length(datakakunin) ) {
    ## experimental data ##
    rfn <- paste("data_kakunin/",datakakunin[i],".csv",sep="")
    dsf <- read.csv(rfn,sep=",",header=F)
    kakunin_list[[i]] <- data.frame( time=dsf[,1],value=log10(dsf[,2]) )
}
# data_list
for (i in 1:length(datakakunin)) {
  if (1==i%%2) {
    assign(sprintf("ccc%d",(i-1)/2),kakunin_list[[i]])
  } else if (0==i%%2) {
    assign(sprintf("dna%d",(i-2)/2),kakunin_list[[i]])
  }
}

##condition5
name_lst <- colnames(sR)
gids <- grep("x[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
c.list <- dataframeC.make(mat)

name_lst <- colnames(sR)
gids <- grep("y[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]

times <- seq(Tmin,Tmax,length=length(gids))
d.list <- dataframeD.make(mat)

check.Cdf5 <- c.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,ccc0,by="time") %>% na.omit()
check.Ddf5 <- d.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,dna0,by="time") %>% na.omit()

plt <- mcmc_plot_func(c.list,d.list,ccc0,dna0)+
  geom_vline(xintercept = 0,color=cmurasaki,lwd=2,alpha=.5)+
  geom_vline(xintercept = 30,color="black",lwd=2,alpha=.5)+
  scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
plot_list[[5]] <- plt
ggsave("fig/condition5.png",plt,w=5,h=5)

## condition6
name_lst <- colnames(sR)
gids <- grep("x1s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
c.list <- dataframeC.make(mat)

name_lst <- colnames(sR)
gids <- grep("y1s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
d.list <- dataframeD.make(mat)

check.Cdf6 <- c.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,ccc1,by="time") %>% na.omit()
check.Ddf6 <- d.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,dna1,by="time") %>% na.omit()

plt <- mcmc_plot_func(c.list,d.list,ccc1,dna1)+
  geom_vline(xintercept = 30,color=cmurasaki,lwd=2,alpha=.5)+
  geom_vline(xintercept = 30,color="black",lwd=2,alpha=.5,lty=2)+
  scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
plot_list[[6]] <- plt
ggsave("fig/condition6.png",plt,w=5,h=5)

## condition7
name_lst <- colnames(sR)
gids <- grep("x2s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
c.list <- dataframeC.make(mat)

name_lst <- colnames(sR)
gids <- grep("y2s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
d.list <- dataframeD.make(mat)

check.Cdf7 <- c.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,ccc2,by="time") %>% na.omit()
check.Ddf7 <- d.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,dna2,by="time") %>% na.omit()

plt <- mcmc_plot_func(c.list,d.list,ccc2,dna2)+
  geom_vline(xintercept = 30,color=cmurasaki,lwd=2,alpha=.5)+
  geom_vline(xintercept = 45,color="black",lwd=2,alpha=.5)+
  scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
plot_list[[7]] <- plt
ggsave("fig/condition7.png",plt,w=5,h=5)

## condition8
name_lst <- colnames(sR)
gids <- grep("x3s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
c.list <- dataframeC.make(mat)

name_lst <- colnames(sR)
gids <- grep("y3s[[:digit:]]+\\.[[:digit:]]+",name_lst)
mat <- sR[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
d.list <- dataframeD.make(mat)

check.Cdf8 <- c.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,ccc3,by="time") %>% na.omit()
check.Ddf8 <- d.list[[3]] %>%
  mutate(x=round(x,1)) %>%
  rename("time"="x") %>%
  left_join(.,dna3,by="time") %>% na.omit()

plt <- mcmc_plot_func(c.list,d.list,ccc3,dna3)+
  geom_vline(xintercept = 45,color=cmurasaki,lwd=2,alpha=.5)+
  geom_vline(xintercept = 30,color="black",lwd=2,alpha=.5)+
  scale_x_continuous(limits = c(Tmin_plot, Tmax_plot))
plot_list[[8]] <- plt
ggsave("fig/condition8.png",plt,w=5,h=5)


fitfig <- plot_list[[1]]+plot_list[[2]]+plot_list[[3]]+plot_list[[4]]+plot_layout(ncol=4)
simufig <- plot_list[[5]]+plot_list[[6]]+plot_list[[7]]+plot_list[[8]]+plot_layout(ncol=4)
ggsave("fig/fig1.png",fitfig,w=16,h=4)
ggsave("fig/fig2A.png",simufig,w=16,h=4)

check.Cdf0$type <- "condition1"
check.Cdf1$type <- "condition2"
check.Cdf2$type <- "condition3"
check.Cdf3$type <- "condition4"
check.Cdf5$type <- "condition5"
check.Cdf6$type <- "condition6"
check.Cdf7$type <- "condition7"
check.Cdf8$type <- "condition8"

check.Ddf0$type <- "condition1"
check.Ddf1$type <- "condition2"
check.Ddf2$type <- "condition3"
check.Ddf3$type <- "condition4"
check.Ddf5$type <- "condition5"
check.Ddf6$type <- "condition6"
check.Ddf7$type <- "condition7"
check.Ddf8$type <- "condition8"


concat.Cdf <- rbind(check.Cdf5,check.Cdf6,check.Cdf7,check.Cdf8) %>%
  mutate(rsq=(y-value)^2)
cdf.rsm <- sqrt(sum(concat.Cdf$rsq)/nrow(concat.Cdf))
concat.Ddf <- rbind(check.Ddf5,check.Ddf6,check.Ddf7,check.Ddf8) %>%
  mutate(rsq=(y-value)^2)
ddf.rsm <- sqrt(sum(concat.Ddf$rsq)/nrow(concat.Ddf))

# --- cccDNA ----
tmp.ccc <- CCC(concat.Cdf$y, concat.Cdf$value, ci = "z-transform", conf.level = 0.95)
ccc_val <- tmp.ccc$rho.c[, "est"]
ccc_ci <- tmp.ccc$rho.c[, c("lwr.ci", "upr.ci")]

plt <- ggplot(concat.Cdf, aes(x = value, y = y)) +

  geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
  geom_point(size = 2, alpha = .9, color = cao) +
  # geom_smooth(method = "lm", formula = y ~ x, lwd = 2, color = cao, fill = caoo) +
  stat_correlation(use_label(c("r", "P")), method = "pearson", size = 4.5) +
  annotate("text", x = 4.63, y = 9,
           label = sprintf("CCC = %.3f (95%% CI: %.3f–%.3f)",
                           ccc_val, ccc_ci[1], ccc_ci[2]),
           hjust = 0, size = 4.5) +
  labs(y = "Estimated value", x = "Observed value") +
  scale_y_continuous(limits = c(4.6, 9.4),
                     breaks = seq(5, 9, 1),
                     labels = c(expression(10^5, 10^6, 10^7, 10^8, 10^9))) +
  scale_x_continuous(limits = c(4.6, 9.4),
                     breaks = seq(5, 9, 1),
                     labels = c(expression(10^5, 10^6, 10^7, 10^8, 10^9))) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  pltsetb

plt
ggsave("fig/fig2B1.png", plt, w = 5, h = 5)


## --- HBV DNA ----
tmp.ccc_dna <- CCC(concat.Ddf$y, concat.Ddf$value, ci = "z-transform", conf.level = 0.95)
ccc_val_dna <- tmp.ccc_dna$rho.c[, "est"]
ccc_ci_dna <- tmp.ccc_dna$rho.c[, c("lwr.ci", "upr.ci")]

plt <- ggplot(concat.Ddf, aes(x = value, y = y)) +

  geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
  geom_point(size = 2, alpha = .9, color = caka) +
  # geom_smooth(method = "lm", formula = y ~ x, lwd = 2, color = caka, fill = cakao) +
  stat_correlation(use_label(c("r", "P")), method = "pearson", size = 4.5) +
  annotate("text", x = 4.63, y = 9,
           label = sprintf("CCC = %.3f (95%% CI: %.3f–%.3f)",
                           ccc_val_dna, ccc_ci_dna[1], ccc_ci_dna[2]),
           hjust = 0, size = 4.5) +
  labs(y = "Estimated value", x = "Observed value") +
  scale_y_continuous(limits = c(4.6, 9.4),
                     breaks = seq(5, 9, 1),
                     labels = c(expression(10^5, 10^6, 10^7, 10^8, 10^9))) +
  scale_x_continuous(limits = c(4.6, 9.4),
                     breaks = seq(5, 9, 1),
                     labels = c(expression(10^5,10^6, 10^7, 10^8, 10^9))) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  pltsetb

ggsave("fig/fig2B2.png", plt, w = 5, h = 5)

###################################################
## Table1
###################################################

load("final.Rdata")
print(fit_mixture)

parset <- as.data.frame(MCMC_mixture$pars)

ci_summary <- parset %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean  = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    lower = hdi(value, credMass = 0.95)[1],
    upper = hdi(value, credMass = 0.95)[2]
  )

print(ci_summary) #table1


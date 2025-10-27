graphics.off()
rm(list=ls())

library(fda)
library(rworldmap)
library(KernSmooth)
library(roahd)

sentiments <- c('positive','neutral','negative')
sbp <- c('polarity','positivity')

year_all <- as.character(2020:2024)
date_month <- as.Date(paste('01', rep(as.character(1:12),length(year_all)), rep(year_all,rep(12,length(year_all)))), format='%d %m %Y')


#### Data upload ####

load("AI_data.RData")


#### Baseline: monthly proportions ####

monthly_baseline_prop <- array(data=0, dim=c(27,60,3), dimnames=list(paesi,1:60,c("pos","neu","neg")))

for(i in 1:27)
  for(j in 1:60)
    monthly_baseline_prop[i,j,] <- monthly_baseline[i,j,]/sum(monthly_baseline[i,j,])

x11(width=21,height=7)
par(mfrow=c(1,3),oma=c(0,0,2,0))

for(j in 1:3){
  plot(date_month[1:36],monthly_baseline_prop[1,1:36,j],type='l',xlab='date',ylab='proportions',ylim=c(0,1),col='white',main=sentiments[j],cex.main=2,cex.lab=2,cex.axis=2,xaxt='n')
  for(i in 1:27)
    lines(date_month[1:36],monthly_baseline_prop[i,1:36,j],type='l',col=i,lwd=2)
  axis(1,at=as.numeric(date_month[c(1,13,25,37)]),labels = c(2020,2021,2022,2023),cex.axis=2)
}
mtext("Baseline", side=3, line=0, outer=TRUE, cex=2)


#### Baseline: polarity and positivity ####

monthly_baseline_sbp <- array(data=0, dim=c(27,60,2), dimnames=list(paesi,1:60,c("pol","pos")))

for(i in 1:27)
  for(j in 1:60){
    monthly_baseline_sbp[i,j,1] <- log(sqrt(monthly_baseline[i,j,1]*monthly_baseline[i,j,3])/monthly_baseline[i,j,2])
    monthly_baseline_sbp[i,j,2] <- log(monthly_baseline[i,j,1]/monthly_baseline[i,j,3])
  }


#### Baseline: smoothing ####

smoothed_monthly_baseline <- array(data = NA, dim = c(401,27,2))

for(i in 1:2)
  for(j in 1:27){
    fit <- locpoly(x=as.numeric(date_month), y=monthly_baseline_sbp[j,,i], bandwidth=50, degree=0, range.x=range(as.numeric(date_month)))
    smoothed_monthly_baseline[,j,i] <- fit$y
  }
abscissa <- fit$x

x11(width=14,height=7)
par(mfrow=c(1,2),oma=c(0,0,2,0))

for(j in 1:2){
  plot(abscissa,smoothed_monthly_baseline[,1,j],type='l',xlab='date',ylab=bquote('p'[.(j)]),ylim=c(-3,1),col='white',main=sbp[j],cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n')
  for(i in 1:27)
    lines(abscissa,smoothed_monthly_baseline[,i,j],type='l',col=i,lwd=2)
  axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels = c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)
}
mtext("Baseline", side=3, line=0, outer=TRUE, cex=2)


#### Baseline: functional PCA ####

data_fd_baseline <- Data2fd(argvals=abscissa, y=smoothed_monthly_baseline)
pca_baseline <- pca.fd(data_fd_baseline,nharm=5,centerfns=TRUE)

# scree plot
x11()
plot(cumsum(pca_baseline$values)[1:20]/sum(pca_baseline$values),xlab='number of PC',ylab='CPV',ylim=c(0.3,1),type='b',main="Scree plot",cex.main=2,cex.lab=1.5)

cumsum(pca_baseline$values)[1]/sum(pca_baseline$values)

# plot of the FPCs as perturbation of the mean
harmfd <- pca_baseline[[1]]
basisfd <- harmfd$basis
rangex <- basisfd$rangeval
x <- seq(rangex[1], rangex[2], length = 50)
fdmat <- eval.fd(x, harmfd)
meanmat <- drop(eval.fd(x, pca_baseline$meanfd))
matharm <- fdmat[, 1, ]
fac <- 2 * sqrt(pca_baseline$values[1])
x11(width=14,height=7)
par(mfrow=c(1,2),oma=c(0,0,2,0))
for (jvar in 1:2) {
  pcmat <- cbind(meanmat[, jvar] + fac * matharm[,jvar], meanmat[, jvar] - fac * matharm[,jvar])
  plot(x, meanmat[, jvar], type = "p", ylab = bquote('p'[.(jvar)]), main = sbp[jvar], xlab = 'date', cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=2, xaxt='n', ylim=c(-3,1))
  axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels=c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)
  points(x, pcmat[, 1], pch = "+",cex=1.2)
  points(x, pcmat[, 2], pch = "-",cex=1.2)
}
mtext('Loadings of first PC (percentage of variability 66.5%)', side=3, line=0, outer=TRUE, cex=2)

rm(harmfd,basisfd,rangex,x,fdmat,meanmat,matharm,fac,jvar,pcmat)

# scatterplot of the scores
x11()
plot(pca_baseline$scores[,1,1],pca_baseline$scores[,1,2],col='white',xlab='Score PC1 polarity',ylab='Score PC1 positivity',cex.lab=1.5)
text(pca_baseline$scores[,1,1],pca_baseline$scores[,1,2],paesi,cex = 1.5)

# geographical map visualization of the scores
labelDF <- data.frame(country = paesi, label = pca_baseline$scores[,1,1])
malMap <- joinCountryData2Map(labelDF, joinCode = "ISO2", nameJoinColumn = "country")
x11()
mapCountryData(malMap, nameColumnToPlot="label", catMethod = "pretty",
               colourPalette="heat", missingCountryCol = gray(.8),
               addLegend=T,mapTitle='Score PC1 polarity',lwd=0.01, mapRegion='europe')

labelDF <- data.frame(country = paesi, label = pca_baseline$scores[,1,2])
malMap <- joinCountryData2Map(labelDF, joinCode = "ISO2", nameJoinColumn = "country")
x11()
mapCountryData(malMap, nameColumnToPlot="label", catMethod = "pretty",
               colourPalette="heat", missingCountryCol = gray(.8),
               addLegend=T,mapTitle='Score PC1 positivity',lwd=0.01, mapRegion='europe')


#### AI: total counts ####

x11(width=10,height=7)
par(oma=c(0,0,2,0))

plot(date_month,monthly_AI[1,,1],type='l',xlab='date',ylab='counts',ylim=c(0,500),col='white',cex.main=2,cex.lab=1.5,cex.axis=1.5)
for(i in 1:27)
  lines(date_month,monthly_AI[i,,1]+monthly_AI[i,,2]+monthly_AI[i,,3],type='l',col=i,lwd=2)
mtext("AI", side=3, line=0, outer=TRUE, cex=2)
axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels = c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)


#### AI: polarity and positivity ####

monthly_AI_sbp <- array(data=0, dim=c(27,60,2), dimnames=list(paesi,1:60,c("pol","pos")))

for(i in 1:27)
  for(j in 1:60){
    monthly_AI_sbp[i,j,1] <- log(sqrt(monthly_AI[i,j,1]*monthly_AI[i,j,3])/monthly_AI[i,j,2])
    monthly_AI_sbp[i,j,2] <- log(monthly_AI[i,j,1]/monthly_AI[i,j,3])
  }


#### AI: smoothing ####

smoothed_monthly_AI <- array(data = NA, dim = c(401,27,2))

for(i in 1:2)
  for(j in 1:27){
    ind <- !is.infinite(monthly_AI_sbp[j,,i]) & !is.nan(monthly_AI_sbp[j,,i])
    fit <- locpoly(x=as.numeric(date_month)[ind], y=monthly_AI_sbp[j,ind,i], bandwidth=50, degree=0, range.x=range(as.numeric(date_month)))
    smoothed_monthly_AI[,j,i] <- fit$y
  }

x11(width=14,height=7)
par(mfrow=c(1,2),oma=c(0,0,2,0))

for(j in 1:2){
  plot(abscissa,smoothed_monthly_AI[,1,j],type='l',xlab='date',ylab=bquote('p'[.(j)]),ylim=c(-3,3),col='white',main=sbp[j],cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n')
  axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels = c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)
  for(i in 1:27)
    lines(abscissa,smoothed_monthly_baseline[,i,j],type='l',col='grey',lwd=2)
  for(i in 1:27)
    lines(abscissa,smoothed_monthly_AI[,i,j],type='l',col=i,lwd=2)
}
mtext("AI", side=3, line=0, outer=TRUE, cex=2)


#### AI: functional PCA ####

data_fd_AI <- Data2fd(argvals=abscissa, y=smoothed_monthly_AI)
pca_AI <- pca.fd(data_fd_AI,nharm=5,centerfns=TRUE)

# scree plot
x11()
plot(cumsum(pca_AI$values)[1:20]/sum(pca_AI$values),xlab='number of PC',ylab='CPV',ylim=c(0.3,1),type='b',main="Scree plot",cex.main=2,cex.lab=1.5)

cumsum(pca_AI$values)[1]/sum(pca_AI$values)

# plot of the FPCs as perturbation of the mean
harmfd <- pca_AI[[1]]
basisfd <- harmfd$basis
rangex <- basisfd$rangeval
x <- seq(rangex[1], rangex[2], length = 50)
fdmat <- eval.fd(x, harmfd)
meanmat <- drop(eval.fd(x, pca_AI$meanfd))
matharm <- fdmat[, 1, ]
fac <- 2 * sqrt(pca_AI$values[1])
x11(width=14,height=7)
par(mfrow=c(1,2),oma=c(0,0,2,0))
for (jvar in 1:2) {
  pcmat <- cbind(meanmat[, jvar] + fac * matharm[,jvar], meanmat[, jvar] - fac * matharm[,jvar])
  plot(x, meanmat[, jvar], type = "p", ylab = bquote('p'[.(jvar)]), main = sbp[jvar], xlab = 'date', cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=2, xaxt='n', ylim=c(-2,2))
  axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels=c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)
  points(x, pcmat[, 1], pch = "+",cex=1.2)
  points(x, pcmat[, 2], pch = "-",cex=1.2)
}
mtext('Loadings of first PC (percentage of variability 47.8%)', side=3, line=0, outer=TRUE, cex=2)

rm(harmfd,basisfd,rangex,x,fdmat,meanmat,matharm,fac,jvar,pcmat)


#### AI: first derivative  ####

smoothed_monthly_AI_d1 <- array(data = NA, dim = c(401,27,2))

for(i in 1:2)
  for(j in 1:27){
    ind <- !is.infinite(monthly_AI_sbp[j,,i]) & !is.nan(monthly_AI_sbp[j,,i])
    fit <- locpoly(x=as.numeric(date_month)[ind], y=monthly_AI_sbp[j,ind,i], bandwidth=100, degree=1, range.x=range(as.numeric(date_month)), drv = 1)
    smoothed_monthly_AI_d1[,j,i] <- fit$y
  }


#### Baseline: first derivative ####

smoothed_monthly_baseline_d1 <- array(data = NA, dim = c(401,27,2))

for(i in 1:2)
  for(j in 1:27){
    ind <- !is.infinite(monthly_baseline_sbp[j,,i])
    fit <- locpoly(x=as.numeric(date_month)[ind], y=monthly_baseline_sbp[j,ind,i], bandwidth=100, degree=1, range.x=range(as.numeric(date_month)), drv = 1)
    smoothed_monthly_baseline_d1[,j,i] <- fit$y
  }


#### Ordering the energy thematic curves with respect to the baseline: centrality and outlyingness through the notion of depth ####
graphics.off()

nrep <- 100   # number of re-samplings

BD_d0 <- array(NA,c(27,nrep+1,2))

for(i in 1:27){
  
  smoothed_monthly_realizations <- array(data = NA, dim = c(401,nrep,2))
  
  for(irep in 1:nrep){
    for(k in 1:2){
      ind <- !is.infinite(monthly_resampling_baseline[irep,i,,k]) & !is.nan(monthly_resampling_baseline[irep,i,,k])
      fit <- locpoly(x=as.numeric(date_month)[ind], y=monthly_resampling_baseline[irep,i,ind,k], bandwidth=50, degree=0, range.x=range(as.numeric(date_month)))
      smoothed_monthly_realizations[,irep,k] <- fit$y
    }

  }
  rm(irep,k,ind,fit)
  
  x11(width=14,height=7)
  par(mfrow=c(1,2),oma=c(0,0,2,0))
  
  for(j in 1:2){
    plot(abscissa,smoothed_monthly_realizations[,1,j],type='l',xlab='date',ylab=bquote('p'[.(j)]),ylim=c(-4,2),col='white',main=sbp[j],cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n')
    axis(1, c(as.numeric(date_month[seq(from = 1, by = 12, length = 5)]), max(as.numeric(date_month)) + 30), 
         c(year_all, "2025"), cex.lab = 1.5, cex.axis = 1.5)
    
    for(irep in 1:nrep)
      lines(abscissa,smoothed_monthly_realizations[,irep,j],type='l',col='gray')
    
    lines(abscissa,smoothed_monthly_AI[,i,j],type='l',col='red',lwd=2)
    lines(abscissa,smoothed_monthly_baseline[,i,j],type='l',col='black',lwd=2)
    
  }
  mtext(paste("Resampling of baseline VS AI in",paesi[i]), side=3, line=-1, outer=TRUE, cex=2)
  
  for(j in 1:2){
    BD_d0[i,,j] <- MBD(rbind(t(smoothed_monthly_realizations[,,j]),smoothed_monthly_AI[,i,j]))
  }
  
}

graphics.off()

x11(width=14,height=8)
par(mfrow=c(2,1))

plot(1,1,xlim=c(1,27),ylim=c(0,0.5),xlab='countries',ylab='depths',xaxt='n')
axis(1,1:27,paesi,cex.lab=1,cex.axis=1)
for(i in 1:27)
  points(rep(i,nrep+1),BD_d0[i,,1],col=c(rep('gray',nrep),'red'),pch=20)
title("Depth of original curves - polarity")

plot(1,1,xlim=c(1,27),ylim=c(0,0.5),xlab='countries',ylab='depths',xaxt='n')
axis(1,1:27,paesi,cex.lab=1,cex.axis=1)
for(i in 1:27)
  points(rep(i,nrep+1),BD_d0[i,,2],col=c(rep('gray',nrep),'red'),pch=20)
title("Depth of original curves - positivity")


# depth of AI: how many standard deviations from the mean?
score_d0 <- array(data = NA,dim = c(27,2))
for(i in 1:27)
  for(j in 1:2){
    m <- mean(BD_d0[i,-101,j])
    s <- sd(BD_d0[i,-101,j])
    score_d0[i,j] <- (m-BD_d0[i,101,j])/s
  }

x11()
plot(score_d0[,1],score_d0[,2],col='white',xlab='polarity',ylab='positivity',asp=1,xlim=c(-1,11),ylim=c(-1,11),xaxt='n',yaxt='n',cex.lab=1.5,cex.axis=1.5)
axis(1,-1:11,-1:11)
axis(2,-1:11,-1:11)
text(score_d0[,1],score_d0[,2],paesi,cex=1.5)
abline(h=0,v=0,col='grey')


#### Statistically testing the significance of temporal trends: boostrap test on the first derivatives of positivity curves ####
graphics.off()

p_ecdf <- function(x,y){
  emp_cdf <- ecdf(x)
  val0 <- emp_cdf(y)
  if(val0<0.5)
    return(val0*2)
  else
    return((1-val0)*2)
}

level <- 0.01

nrep <- 100   # number of re-samplings

pval <- array(NA,c(27,401))

for(i in 1:27){
  
  smoothed_monthly_realizations_d1 <- array(data = NA, dim = c(401,nrep))
  
  for(irep in 1:nrep){
    ind <- !is.infinite(monthly_resampling_baseline[irep,i,,2]) & !is.nan(monthly_resampling_baseline[irep,i,,2])
    fit_d1 <- locpoly(x=as.numeric(date_month)[ind], y=monthly_resampling_baseline[irep,i,ind,2], bandwidth=100, degree=1, range.x=range(as.numeric(date_month)), drv = 1)
    smoothed_monthly_realizations_d1[,irep] <- fit_d1$y
  }
  rm(irep,ind,fit_d1)
  
  x11(width=14,height=7)
  par(mfrow=c(1,2),oma=c(0,0,2,0))
  
  plot(abscissa,smoothed_monthly_realizations_d1[,1],type='l',xlab='date',ylab='',ylim=range(smoothed_monthly_AI_d1),col='white',cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n',main='First derivative')
  axis(1, c(as.numeric(date_month[seq(from = 1, by = 12, length = 5)]), max(as.numeric(date_month)) + 30), 
       c(year_all, "2025"), cex.lab = 1.5, cex.axis = 1.5)
  for(irep in 1:nrep)
    lines(abscissa,smoothed_monthly_realizations_d1[,irep],type='l',col='gray')
  
  lines(abscissa,smoothed_monthly_baseline_d1[,i,2],type='l',lwd=2)
  lines(abscissa,smoothed_monthly_AI_d1[,i,2],type='l',col='red',lwd=2)
  abline(h=0)
  
  mtext(paste('Result of test on first derivative of positivity for',paesi[i]), side=3, line=-1, outer=TRUE, cex=2)
  
  for(j in 1:401){
    pval[i,j] <- p_ecdf(x = smoothed_monthly_realizations_d1[j,], y = smoothed_monthly_AI_d1[j,i,2])
  }
  
  semistep <- diff(abscissa)[1]/2
  regions_signif_1 <- NULL
  for(j in 1:401)
    if(pval[i,j]<level)
      regions_signif_1 <- rbind(regions_signif_1,c(abscissa[j]-semistep,abscissa[j]+semistep))
  
  if(!is.null(regions_signif_1))
    for(ii in 1:dim(regions_signif_1)[1])
      rect(regions_signif_1[ii,1],par("usr")[3],regions_signif_1[ii,2],par("usr")[4],col=rgb(0,0,0,alpha=0.3),border=NA) 
  
  plot(abscissa,smoothed_monthly_baseline[,i,2],type='l',xlab='date',ylab='',ylim=range(smoothed_monthly_AI[,i,2],smoothed_monthly_baseline[,i,2]),col='black',lwd=2,cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n',main='Original curves')
  axis(1,as.numeric(date_month[seq(from=1,by=12,length=5)]),year_all,cex.lab=1.5,cex.axis=1.5)
  lines(abscissa,smoothed_monthly_AI[,i,2],type='l',col='red',lwd=2)
  
  if(!is.null(regions_signif_1))
    for(ii in 1:dim(regions_signif_1)[1])
      rect(regions_signif_1[ii,1],par("usr")[3],regions_signif_1[ii,2],par("usr")[4],col=rgb(0,0,0,alpha=0.3),border=NA) 
  
}

graphics.off()

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1,cex.axis=1.5)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

sum_significant <- colSums(pval<level)
colorscale <- gray.colors (max(sum_significant),start = 1, end = 0.3)

j <- 2

x11(width=8.5,height=7)
layout(mat = matrix(1:2, nrow = 1, ncol = 2), heights = 7, widths = c(7, 1.5))

plot(abscissa,smoothed_monthly_AI[,1,j],type='l',xlab='date',ylab='',ylim=c(-2,3),col='white',main='',cex.main=2,cex.lab=1.5,cex.axis=1.5,xaxt='n')
axis(1,at=as.numeric(date_month[c(1,13,25,37,49,60)]),labels = c(2020,2021,2022,2023,2024,2025),cex.axis=1.5)
for(ii in 1:400)
  rect(abscissa[ii]-semistep,par("usr")[3],abscissa[ii]+semistep,par("usr")[4],col=colorscale[sum_significant[ii]],border=NA)
for(i in 1:27)
  lines(abscissa,smoothed_monthly_AI[,i,j],type='l',col=i,lwd=2)

color.bar(colorscale, 1,max(sum_significant),max(sum_significant))

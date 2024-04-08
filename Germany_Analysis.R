list.of.packages <- c("devtools","lubridate","reshape2","zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(devtools)
library(lubridate)
library(reshape2)
library(zoo)  

source_url("https://github.com/lcyjames/srtaCFR/blob/main/CoreFunctions.R?raw=TRUE")
dat<-read.csv("https://raw.githubusercontent.com/lcyjames/srtaCFR/main/Germanybyage.csv")

dat<-subset(dat,(is.na(dat$Cases)==FALSE)&(is.na(dat$AgeGroup)==FALSE))
dat$Date<-as.Date(dat$Date, format="%Y-%m-%d")
dat$Date1<-as.numeric(dat$Date)
month.labels<-rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),4)

(max(dat$Date1)-min(dat$Date1)+1)==nrow(dat)/3

start.date<-min(dat$Date)
end.date<-as.Date("2022-06-01")

date.seq<-c()#labelled all months' start date
for(i in c(2020:2023)){
  for(j in c(1:12)){
    date.seq<-c(date.seq, as.Date(paste0(i,"-",j,"-01")))
  }}

mar.seq<-c()
for(i in c(2020:2023)){
  j<-3
    mar.seq<-c(mar.seq, as.Date(paste0(i,"-",j,"-01")))
  }

sub    <-which((dat$Date1>=as.numeric(start.date))&(dat$Date1<=as.numeric(end.date)))
dat    <-dat[sub,] 

length(unique(dat$Date))
length(unique(dat$Date1))

case_cum_ts <-dcast(dat,Date1~AgeGroup,value.var = "Cases", fun.aggregate = sum)
death_cum_ts<-dcast(dat,Date1~AgeGroup,value.var = "Deaths", fun.aggregate = sum)

case_cum_ts_ma<-death_cum_ts_ma<-c()
for(i in c(1:3)){
  case_cum_ts_ma <-cbind(case_cum_ts_ma,round(rollmean(case_cum_ts[,i+1], k = 7)))
  death_cum_ts_ma<-cbind(death_cum_ts_ma,round(rollmean(death_cum_ts[,i+1], k = 7)))
}

case_ts <-apply(case_cum_ts_ma,2,diff)
death_ts<-apply(death_cum_ts_ma,2,diff)
death_ts[which(death_ts<0,arr.ind = T)]<-0

ct_sum<-apply(case_ts,1,sum)
dt_sum<-apply(death_ts,1,sum)

CT<-apply(case_cum_ts_ma,1,sum)
DT<-apply(death_cum_ts_ma,1,sum)

png(paste0("CaseDeaths.png"),res=300,width = 4500,height=3000)
par(mar = c(7,5,4,5))
plot(ct_sum, ylab="Confirmed Cases",xaxt="n", xlab="", type="l",lwd=2.5, ylim=c(0,max(ct_sum)),
     cex.lab=2,cex.axis = 1.5)
par(new = TRUE)
plot(rollmean(dt_sum, k = 7), axes = FALSE, bty = "n", ylim = c(0,2*max(dt_sum)),
     xlab = "", ylab = "", type="l", lty=2, lwd=2.5, col="red", cex.lab=2, cex.axis = 1.5)
axis(side=4, at = seq(0,4000,100),cex.axis=1.5)
axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
legend("topleft",col=c("black","red"),lty=c(1,2),
       lwd = c(2.5,2.5),bty = "n",x.intersp=0.5,text.width=150,legend = c("Confirmed case","Death"),ncol=2,cex=1.5)
mtext(text = "Deaths", side=4, line=3,cex = 2)
dev.off()

#from https://www.statista.com/statistics/454349/population-by-age-group-germany/
pop2021<-c(0.2,0.51,0.29)
fit_G<-rtaCFR.EST(ct = ct_sum, dt = dt_sum)
fit_G1<-srtaCFR.EST(ct_mat = case_ts, dt_mat = death_ts, q_mat = t(replicate(nrow(case_ts),pop2021)))

############End of estimation procedure#################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set pathway for storing figures
#Figure 1 sensitivity of rtaCFR
png(paste0("Sensitivity1.png"),res=300,width = 2500,height=1800)
par(mar = c(7,5,2,1))
#par(mfrow=c(2,2),mar = c(7,5,2,1))

colfunc <-colorRampPalette(c("black", "purple"))
colgrad <-colfunc(100)
lamb_grid<-seq(20000,60000,length.out=100)
legend_image <- as.raster(rev(matrix(colgrad, ncol=1)))

for(r in c(1:100)){
  cd_l         <-which.min(abs(fit_G$sol_path$lambda-lamb_grid[r]))
  rtaCFR_lambda<-as.numeric(fit_G$p_matrix[,cd_l])
  
  rtaCFR_lambda<-ifelse(dt_sum<=1,0,rtaCFR_lambda)
  
  ks0<-ksmooth(c(1:length(ct_sum)), rtaCFR_lambda, kernel = "normal",bandwidth = 15)
  if(r==1){
    plot(ks0, ylab="rtaCFR", xlab="",xaxt="n",xlim=c(0,801), ylim=c(0,0.15), type="l", col=colgrad[r],
         lwd=2.5,cex.lab=2,cex.axis = 1.5,cex.main=2)
    axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
    axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
  }
  else{lines(ks0, type="l", lwd=1.5, col=colgrad[r])
  }
  abline(h=0)
}

par(fig = c(0.7, 0.9, 0.45, 0.95), new = T)  

plot(c(0, 1.5), c(0, 1), type='n', axes=F, xlab='', ylab='')      ## ini plot2
rasterImage(legend_image, 0, 0, 1, 1)                           ## the gradient
lbsq <- seq.int(0, 1, l=5)                                      ## seq. for labels
axis(4, at=lbsq, pos=1, labels=F, col=0, col.ticks=1, tck=-.1)  ## axis ticks
mtext(seq(20000,60000,10000), 4, -.5, at=lbsq, las=2, cex=1.2)  ## tick labels
mtext(expression(lambda), 3, -.125, cex=2, adj=.1, font=2)    ## title
#par(op)                                                         ## reset par
on.exit(par(op))
dev.off()

#Sensitivity of srtaCFR
png(paste0("Sensitivity2.png"),res=300,width = 2500,height=1800)
par(mar = c(7,5,2,1))

colfunc <-colorRampPalette(c("black", "purple"))
colgrad <-colfunc(100)
lamb_grid1<-seq(10000,50000,length.out=100)
legend_image <- as.raster(rev(matrix(colgrad, ncol=1)))

for(r in c(1:100)){
  p_hat_by_age<-c()
  for(j in c(1:3)){
    cd_l<-which.min(abs(fit_G1$sol_path_tensor[[j]]$lambda-lamb_grid1[r]))
    p_hat_by_age<-cbind(p_hat_by_age,as.numeric(fit_G1$p_tensor[[j]][,cd_l]))
  }
  p_hat_by_age[,1]<-ifelse(dt_sum<=1,0,p_hat_by_age[,1])
  p_hat_by_age[,2]<-ifelse(dt_sum<=1,0,p_hat_by_age[,2])
  p_hat_by_age[,3]<-ifelse(dt_sum<=1,0,p_hat_by_age[,3])
  
  ks1  <-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,1], kernel = "normal",bandwidth = 15)$y
  ks2  <-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,2], kernel = "normal",bandwidth = 15)$y
  ks3  <-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,3], kernel = "normal",bandwidth = 15)$y
  ks4  <-apply(cbind(ks1,ks2,ks3)*t(replicate(nrow(case_ts),pop2021)),1,sum)
  
  if(r==1){
    plot(x = c(1:length(ct_sum)),y = ks4, ylab="srtaCFR", xlab="",xaxt="n",xlim=c(0,801), ylim=c(0,0.15), 
         type="l", col=colgrad[r], lwd=2.5,cex.lab=2,cex.axis = 1.5,cex.main=2)
    axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
    axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
  }
  else{lines(x = c(1:length(ct_sum)),y = ks4, type="l", lwd=1.5, col=colgrad[r])
  }
  abline(h=0)
}

par(fig = c(0.7, 0.9, 0.45, 0.95), new = T)  

plot(c(0, 1.5), c(0, 1), type='n', axes=F, xlab='', ylab='')      ## ini plot2
rasterImage(legend_image, 0, 0, 1, 1)                           ## the gradient
lbsq <- seq.int(0, 1, l=5)                                      ## seq. for labels
axis(4, at=lbsq, pos=1, labels=F, col=0, col.ticks=1, tck=-.1)  ## axis ticks
mtext(seq(10000,50000,10000), 4, -.5, at=lbsq, las=2, cex=1.2)  ## tick labels
mtext(expression(lambda), 3, -.125, cex=2, adj=.1, font=2)    ## title
#par(op)                                                         ## reset par
on.exit(par(op))
dev.off()

#Age-specific rtaCFRs
png(paste0("AgeSpecificFatality.png"),res=300,width = 2500,height=1800)
par(mar = c(7,5,2,1))

p_hat_by_age<-fit_G1$p_gp_spec
p_hat_by_age[,1]<-ifelse(dt_sum<=1,0,p_hat_by_age[,1])
p_hat_by_age[,2]<-ifelse(dt_sum<=1,0,p_hat_by_age[,2])
p_hat_by_age[,3]<-ifelse(dt_sum<=1,0,p_hat_by_age[,3])

ks1_30k<-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,1], kernel = "normal",bandwidth = 15)$y
ks2_30k<-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,2], kernel = "normal",bandwidth = 15)$y
ks3_30k<-ksmooth(c(1:length(ct_sum)), p_hat_by_age[,3], kernel = "normal",bandwidth = 15)$y
ks4_30k<-apply(cbind(ks1,ks2,ks3)*t(replicate(nrow(case_ts),pop2021)),1,sum)

plot(ks1_30k, ylab="Age-group-specific fatality rates", xlab="", type="l", xaxt="n", 
     lwd=2.5, cex.lab=2,cex.axis = 1.5,cex.main=2, col="darkgreen", xlim=c(0,801), ylim=c(0,0.55)) #ylim=c(0,max(ks3)))
lines(ks2_30k, lty=2, type="l", lwd=2.5, col="red")
lines(ks3_30k, lty=3, type="l", lwd=2.5, col="blue")
axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
legend("topleft",col=c("darkgreen","red","blue"),lty=c(1,2,3),
       lwd = c(2.5,2.5,2.5),bty = "n",x.intersp=0.5,text.width=100,
       legend = c(expression(paste("rtaCFR"^"(1)")),
                  expression(paste("rtaCFR"^"(2)")),
                  expression(paste("rtaCFR"^"(3)"))),ncol = 3,cex=1.5)
dev.off()


#Both rtaCFR and srtaCFR
png(paste0("BothFatality.png"),res=300,width = 2500,height=1800)
par(mar = c(7,5,2,1))

rtaCFR_lambda<-fit_G$p_hat 
rtaCFR_lambda<-ifelse(dt_sum<=1,0,rtaCFR_lambda)

ks0_40k<-ksmooth(c(1:length(ct_sum)), rtaCFR_lambda, kernel = "normal",bandwidth = 15)

CFR    <-ksmooth(c(1:length(ct_sum)), ifelse(dt_sum<=1,0,DT/CT), kernel = "normal",bandwidth = 15)

plot(ks0_40k, ylab="Fatality rates", xlab="", type="l",xaxt="n",xlim=c(0,801), ylim=c(0,0.15),
     lwd=2.5,cex.lab=2,cex.axis = 1.5,cex.main=2,col="black")
lines(ks4_30k, type="l", lwd=2.5, col="red",lty=2)
lines(CFR, type="l", lwd=2.5, col="blue",lty=3)
axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
legend("topleft",col=c("black","red","blue"),lty=c(1,2,3),
       lwd = c(2.5,2.5,2.5),bty = "n",x.intersp=0.5,text.width=100,legend = c("rtaCFR","srtaCFR","CFR"),ncol=3,cex=1.5)

dev.off()

#Temporal age distribution
png(paste0("AgeDist.png"),res=300,width = 4000,height=2500)
par(mar = c(7,5,4,1), oma=c(0,0,0,12.5))
colfunc2 <-colorRampPalette(c("coral", "cornflowerblue"))
prop1<-t(case_ts/replicate(3,apply(case_ts,1,sum)))[c(3:1),][1,]
prop2<-t(case_ts/replicate(3,apply(case_ts,1,sum)))[c(3:1),][2,]+prop1
prop3<-t(case_ts/replicate(3,apply(case_ts,1,sum)))[c(3:1),][3,]+prop2
plot(prop3, type = "l", ylim = c(0, 1), ylab="Proportion", xlab="",
     cex.lab=2,cex.axis = 1.5,cex.main=2,xaxt='n')
polygon(c(1, seq(prop3), length(prop3)), c(0, prop3, 0), col = colfunc2(3)[1])
lines(prop2, type = "l")
polygon(c(1, seq(prop2), length(prop2)), c(0, prop2, 0), col = colfunc2(3)[2])
lines(prop1, type = "l")
polygon(c(1, seq(prop1), length(prop1)), c(0, prop1, 0), col = colfunc2(3)[3])
axis(side=1, at = (date.seq-as.numeric(start.date)), labels = month.labels,las=2,hadj = 0.8,cex.axis=1.5)
axis(side=1, at = (mar.seq-as.numeric(start.date)), labels = c(2020:2023),padj = 2.5,cex.axis=1.5,lwd.tick=0)
legend(par('usr')[2], par('usr')[4],xpd=NA, inset=c(-0.2,0),fill=colfunc2(3),
       bty = "n",x.intersp=0.5,text.width=150,
       legend = c("0 - 20 years old","21 - 60 years old","60+ years old"),ncol=1,cex=1.5)
dev.off()

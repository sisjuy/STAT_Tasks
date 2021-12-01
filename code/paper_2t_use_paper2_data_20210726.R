#setwd("C:\\Users\\HCresearch\\Dropbox\\computing\\R_program")
setwd("C:\\Users\\news\\Dropbox\\computing\\R_program")

####I. read in cleaned data
#restricting to California, interstates accident, whole blood tests
#if later people ask 
bac_dat=read.csv("data.csv",header=TRUE)
#bac_dat=read.csv("FARS_2009_BACdata_subset.csv",header=TRUE) #Texas
head(bac_dat)

age=bac_dat$AGE
alc_res=bac_dat$ALC_RES/100 #noted 20150130!
bac_dat$ALC_RES
hist(age)
table(age)
table(alc_res) #already restricted to (0,96) #take out 0 because we look at 'drunker drivers'

length(age)
#[1] 125
length(alc_res)
#[1] 125
?hist
#California
par(mfrow=c(2,1))
hist(alc_res[age<30],prob=TRUE,xlim=c(0,40),ylim=c(0,0.06))
hist(alc_res[age>=30],prob=TRUE,xlim=c(0,40),ylim=c(0,0.06))

####II. California
set.seed(111)
young=alc_res[age<30]#+runif(length(alc_res[age<30]),-0.1,0.1) #20140613: no need to break ties
old=alc_res[age>=30]#+runif(length(alc_res[age>=30]),-0.1,0.1) #20140613: no need to break ties
data=cbind(c(young,old),c(young^0.5,old),c(rep(1,times=length(young)),rep(2,times=length(old))))
#data=cbind(c(old,young),c(old,young^0.5),c(rep(1,times=length(old)),rep(2,times=length(young))))

length(young)
#67
length(old)
#58
young
par(mfrow=c(2,1))
hist(young,prob=TRUE,xlim=c(0,0.9),ylim=c(0,0.05))
hist(old,prob=TRUE,xlim=c(0,40),ylim=c(0,0.05))

#observed ecdfs
par(mfrow=c(1,1))
plot(ecdf(young))
lines(ecdf(old),lty=3,col="red")

#weighted ecdfs
library(Hmisc)
wtcdf1=wtd.Ecdf(sort(young), weights=1/sort(young^0.5)/(sum(1/sort(young^0.5))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)

wtcdf1

plot(stepfun(wtcdf1$x[-1],wtcdf1$ecdf), verticals = TRUE, do.points = FALSE, xlim=c(0,40),ylab="")
wtcdf2=wtd.Ecdf(sort(old), weights=1/sort(old)/(sum(1/sort(old))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)
lines(stepfun(wtcdf2$x[-1],wtcdf2$ecdf), verticals = TRUE, do.points = FALSE, lty=2, col="red")
title(ylab=expression(paste(tilde(F),"(x)")))

#check wtcdf1, then
> abline(v=1)
> abline(h=0.1015373)
> abline(v=2)
> abline(h=0.1374361)
> abline(v=30)

#check wtcdf2, then 
> abline(h=0.1611889)
> abline(h=0.3223778)
> abline(v=38)

####III. analysis
#first read functions from file "\paper2_dataanalysis_functions_20140525.R"
#which is based on file "\paper2_addsump_F_app1_dM_p2nd_2_beta_20140523.R" functions
#but with unnecessary arguments and truth-related calculations removed

source("paper2_dataanalysis_functions_20140525.R")
Sys.time()
set.seed(111) #20140619
result=teststat(data=data,c(8,30),r=c(0.5,1),nboot=1000) #young=1 old=2
#result=teststat(data=data,c(8,30),r=c(1,0.5),nboot=1000) #old=1 young=2
Sys.time()
result$suptest
#4.458534
result$suptest_pval
#[1] 0.11 #[1]0.109 #new in 20140629 after from largest lambda and setseed
result$suptest_Wald
#3.538191
result$suptest_pval_Wald
#[1] 0.157 #[1] 0.168 #new in 20140629 after from largest lambda and setseed
data_noweight=cbind(c(young,old),rep(1,times=length(alc_res)),c(rep(1,times=length(young)),rep(2,times=length(old))))
result_now=teststat(data=data_noweight,c(8,30),r=c(0,0),nboot=1000)
result_now$suptest
#0.3432606 #[1] 0.3432584
result_now$suptest_pval
#[1] 0.844 #[1] 0.841

#trunc=c(8,30)
#r=c(0.5,1)
#nboot=1000


####III.1 try other weight functions

#(1) young weight 0.25
data=cbind(c(young,old),c(young^0.25,old),c(rep(1,times=length(young)),rep(2,times=length(old))))

Sys.time()
result=teststat(data=data,c(8,30),r=c(0.25,1),nboot=1000) #young=1 old=2
Sys.time()
result$suptest
#8.324941
result$suptest_pval
#[1] 0.016
result$suptest_Wald
#5.764885
result$suptest_pval_Wald
#[1] 0.047

#(2) young weight 0.75
data=cbind(c(young,old),c(young^0.75,old),c(rep(1,times=length(young)),rep(2,times=length(old))))

Sys.time()
result=teststat(data=data,c(8,30),r=c(0.75,1),nboot=1000) #young=1 old=2
Sys.time()
result$suptest
#1.484372
result$suptest_pval
#[1] 0.401
result$suptest_Wald
#1.360724
result$suptest_pval_Wald
#[1] 0.428

#(3) a series of weights
wtset=seq(0.1,1,by=0.05)
p_ours=1:length(wtset)*0
p_Wald=1:length(wtset)*0
Sys.time()
for (wi in 1:length(wtset)){
	data=cbind(c(young,old),c(young^wtset[wi],old),c(rep(1,times=length(young)),rep(2,times=length(old))))
	result=teststat(data=data,c(8,30),r=c(wtset[wi],1),nboot=1000) #young=1 old=2
	p_ours[wi]=result$suptest_pval
	p_Wald[wi]=result$suptest_pval_Wald
}
Sys.time()
plot(wtset,p_ours,type="o",ylim=c(0,1))
points(wtset,p_Wald,type="o",lty=2)
abline(h=0.05,lty=3)
abline(h=0.1,lty=3)

####################IV. Output plots####################


############20150124 subfigure IV.3:  ecdfs---observed data and weighted ecdfs (using the one-sample formula)
#20150130: /100 for certain xaxis quantities

#####inputs: fonts of labels
label_size=1.8*1.8
tcl_size=0.3*1.8
axis_size=1.5*1.6
position_lab=3.7*1.8
bot_left_mar=6*1.8
lwd=1.5*1.8
right_mar=0.5 #added 20150127
#####

#setwd("C:\\Users\\HCresearch\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")
setwd("C:\\Users\\hchang\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")

pdf("Figure3a.pdf",width=10.4,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,right_mar),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

#ecdf:
young_ecdf=ecdf(young)
old_ecdf=ecdf(old)

plot(old_ecdf, 
verticals = TRUE, do.points = FALSE,main="",xlim=c(0,40/100),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
lines(young_ecdf, verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,40/100,by=10/100),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC (g/dL)",cex.lab=label_size) #added 20150124
title(ylab="ecdf",cex.lab=label_size,line=position_lab)

dev.off()

pdf("Figure3b.pdf",width=10.4,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,right_mar),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)


#weighted ecdfs
library(Hmisc)
young_wtcdf=wtd.Ecdf(sort(young), weights=1/sort(young^0.5)/(sum(1/sort(young^0.5))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)
old_wtcdf=wtd.Ecdf(sort(old), weights=1/sort(old)/(sum(1/sort(old))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)

plot(stepfun(old_wtcdf$x[-1],old_wtcdf$ecdf), 
verticals = TRUE, do.points = FALSE,main="",xlim=c(0,40/100),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
abline(h=1,lty=2,col="gray",lwd=lwd)
lines(stepfun(young_wtcdf$x[-1],young_wtcdf$ecdf), verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,40/100,by=10/100),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC (g/dL)",cex.lab=label_size)
title(ylab=expression(paste(tilde(F))),cex.lab=label_size,line=position_lab)


dev.off()



############IV.1: histograms

#####inputs: fonts of labels
label_size=1.8
tcl_size=0.3
axis_size=1.5
position_lab=3.7
bot_left_mar=5
lwd=1.5
#####

setwd("C:\\Users\\HCresearch\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")

pdf("paper2_fig_dat.pdf",width=5.2,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(2,1),mar=c(bot_left_mar,bot_left_mar,1,0),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

hist(young,prob=TRUE,xlim=c(0,40),ylim=c(0,0.05),ann=FALSE,axes=FALSE,main="",xaxs="i",yaxs="i")
axis(1, at=seq(0,40,by=10),las=1,tcl=-tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,0.05,by=0.01),las=1,tcl=-tcl_size,cex.axis=axis_size) #side=1: bottom
title(ylab="Density",cex.lab=label_size,line=position_lab)

hist(old,prob=TRUE,xlim=c(0,40),ylim=c(0,0.05),ann=FALSE,axes=FALSE,main="",xaxs="i",yaxs="i")
axis(1, at=seq(0,40,by=10),las=1,tcl=-tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,0.05,by=0.01),las=1,tcl=-tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC values",cex.lab=label_size)
title(ylab="Density",cex.lab=label_size,line=position_lab)


dev.off()

############IV.2: ecdfs---observed data 

#####inputs: fonts of labels
label_size=1.8
tcl_size=0.3
axis_size=1.5
position_lab=3.7
bot_left_mar=5
lwd=1.5
#####

setwd("C:\\Users\\HCresearch\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")

pdf("paper2_fig_ecdf.pdf",width=5.2,height=5.2)
#setEPS()
#postscript("paper1_fig3.eps", height=5.2, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,0),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

#ecdf:
young_ecdf=ecdf(young)
old_ecdf=ecdf(old)

plot(old_ecdf, verticals = TRUE, do.points = FALSE,main="",xlim=c(0,40),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
lines(young_ecdf, verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,40,by=10),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC values",cex.lab=label_size)
title(ylab="ecdf",cex.lab=label_size,line=position_lab)

dev.off()


############IV.3:  ecdfs---observed data and weighted ecdfs (using the one-sample formula)

#####inputs: fonts of labels
label_size=1.8
tcl_size=0.3
axis_size=1.5
position_lab=3.7
bot_left_mar=6  #5
lwd=1.5
#####

setwd("C:\\Users\\HCresearch\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")

pdf("paper2_fig_wecdf.pdf",width=5.2,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(2,1),mar=c(bot_left_mar,bot_left_mar,1,0),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

#ecdf:
young_ecdf=ecdf(young)
old_ecdf=ecdf(old)

plot(old_ecdf, 
verticals = TRUE, do.points = FALSE,main="",xlim=c(0,40),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
lines(young_ecdf, verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,40,by=10),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
#title(xlab="BAC values",cex.lab=label_size)
title(ylab="ecdf",cex.lab=label_size,line=position_lab)

#weighted ecdfs
library(Hmisc)
young_wtcdf=wtd.Ecdf(sort(young), weights=1/sort(young^0.5)/(sum(1/sort(young^0.5))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)
old_wtcdf=wtd.Ecdf(sort(old), weights=1/sort(old)/(sum(1/sort(old))), 
         #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
         normwt=FALSE, na.rm=TRUE)

plot(stepfun(old_wtcdf$x[-1],old_wtcdf$ecdf), 
verticals = TRUE, do.points = FALSE,main="",xlim=c(0,40),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
abline(h=1,lty=2,col="gray",lwd=lwd)
lines(stepfun(young_wtcdf$x[-1],young_wtcdf$ecdf), verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,40,by=10),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC values",cex.lab=label_size)
title(ylab=expression(paste(tilde(F))),cex.lab=label_size,line=position_lab)


dev.off()

############IV.4: p values vs different weight functions (20150124 update size and name)

#####inputs: fonts of labels
label_size=1.8
tcl_size=0.3
axis_size=1.5
position_lab=3.7
bot_left_mar=5
lwd=1.5
#####

setwd("C:\\Users\\HCresearch\\Dropbox\\papers_mine\\paper_ELSOb_20130702\\figures")

#pdf("paper2_fig_pvswt.pdf",width=5.2,height=5.2) #before 20150124
pdf("Figure4.pdf",width=10.4,height=10.4)

#setEPS()
#postscript("paper1_fig3.eps", height=5.2, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,0),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

wtset=seq(0.1,1,by=0.05)
p_ours=1:length(wtset)*0
p_Wald=1:length(wtset)*0
Sys.time()
for (wi in 1:length(wtset)){
	data=cbind(c(young,old),c(young^wtset[wi],old),c(rep(1,times=length(young)),rep(2,times=length(old))))
	result=teststat(data=data,c(8,30),r=c(wtset[wi],1),nboot=1000) #young=1 old=2
	p_ours[wi]=result$suptest_pval
	p_Wald[wi]=result$suptest_pval_Wald
}
Sys.time()

#tried ,type="o",pch=19, but look strange in latex even when figure very large... so will just do lines...
plot(wtset,p_ours, 
,xlim=c(0,1),ylim=c(0,1),type="l",lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,tcl=tcl_size,frame.plot=FALSE)#,cex=0.7
abline(h=0.05,lty=3,lwd=lwd)
abline(h=0.1,lty=3,lwd=lwd)
points(wtset,p_Wald,type="l",lwd=lwd,lty=2)
axis(1, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=c(0,0.1, seq(0.2,1,by=0.2)),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="r" #expression(paste("r in ",w[y](x)==x^r))
,cex.lab=label_size)
title(ylab="p-value",cex.lab=label_size,line=position_lab)

dev.off()



########################OLD###########################
####I. read in cleaned data
bac_dat=read.csv("C:\\Users\\HCresearch\\Dropbox\\computing\\R_program\\FARS_2009_BACdata.csv",header=TRUE)
age=bac_dat[,1]
alc_res=bac_dat[,2]

####II. find subset such that
that separates ALC_RES (ALCOHOL TEST RESULTS) best
####(because maybe after all these years, drinking age decrease or increase...)


#(1) a function that calculates p-value based on bootstrap sample
#see simulate_H0_boot in file "paper2_test_20130913.R" bootstrap variable names
mean(test_star>test[i])


####II. find age cut that separates ALC_RES (ALCOHOL TEST RESULTS) best
####(because maybe after all these years, drinking age decrease or increase...)

#take bulks of ages (>90 people)
ages=20:50
med_diffs=cbind(ages,1:length(ages)*0)
tstats=cbind(ages,1:length(ages)*0)

for (ai in 1:length(ages)){
	tstats[ai,2]=abs(t.test(alc_res[age<ages[ai]],alc_res[age>=ages[ai]])$statistic) #,var.equal=TRUE
	med1=median(alc_res[age<ages[ai]])
	med2=median(alc_res[age>=ages[ai]])
	med_diffs[ai,2]=abs(med1-med2)
}

med_diffs[which.max(med_diffs[,2]),] #largest difference: age 17

tstats[sort(tstats[,2],index.return=TRUE)$ix,]
tstats[which.max(tstats[,2]),] #largest difference: equalvar: 24, unequal var: 21 (but 24 is second)
#-->will do >=25 then... (25 ranked 5th largest, ok)
#indeed drinking age decreases!

par(mfrow=c(2,1))
range(alc_res)
hist(alc_res[age<24],prob=TRUE,xlim=c(1,79),ylim=c(0,0.05))
hist(alc_res[age>=24],prob=TRUE,xlim=c(1,79),ylim=c(0,0.05))


#############COULDN'T CONVERGE BECAUSE DIDN'T READ IN PACKAGES###############
#############The following investigations are in vain#####################

#Texas
par(mfrow=c(2,1))
range(alc_res) #[1]  2 41
hist(alc_res[age<=30],prob=TRUE,xlim=c(1,41),ylim=c(0,0.06))
hist(alc_res[age>30],prob=TRUE,xlim=c(1,41),ylim=c(0,0.06))
####II. Texas: two groups
mean(alc_res)
sd(alc_res)

set.seed(111)
young=alc_res[age<=30]+rnorm(length(alc_res[age<=30]),0,0.2) #+unif to break the ties...
old=alc_res[age>30]+rnorm(length(alc_res[age>30]),0,0.2)

par(mfrow=c(2,1))
range(alc_res) #[1]  2 41
hist(young,prob=TRUE,xlim=c(1,41),ylim=c(0,0.06))
hist(old,prob=TRUE,xlim=c(1,41),ylim=c(0,0.06))

data=cbind(c(young,old),c(young^0.5,old),c(rep(1,times=length(young)),rep(2,times=length(old))))

####II.1 two groups switched

data=cbind(c(old,young),c(old,young^0.5),c(rep(1,times=length(old)),rep(2,times=length(young))))

####II.2 same sample size?
set.seed(111)
young=alc_res[age<=30]+rnorm(length(alc_res[age<=30]),0,0.2) #+unif to break the ties...
old=alc_res[age>30]+rnorm(length(alc_res[age>30]),0,0.2)
old=old[1:48]

data=cbind(c(old,young),c(old,young^0.5),c(rep(1,times=length(old)),rep(2,times=length(young))))

####II.3 same sample size+standardize?
set.seed(111)
young=alc_res[age<=30]+rnorm(length(alc_res[age<=30]),0,0.2) #+unif to break the ties...
old=alc_res[age>30]+rnorm(length(alc_res[age>30]),0,0.2)
old=old[1:49]

data_std=(c(old[-14],young)-min(c(old[-14],young))+0.1)/(max(c(old[-14],young))+1.2)

data=cbind(data_std,c(data_std[1:48],data_std[49:96]^0.5),c(rep(1,times=length(old[-14])),rep(2,times=length(young))))

####II.4 restrict to >=8 & <=30 (to-do: find out why taking out extreme values
still does not work)
set.seed(111)
young0=alc_res[age<=30]
old0=alc_res[age>30] #[-c(8,14)]

young=young0[young0>=8 & young0<=30]
old=old0[old0>=8 & old0 <=30]

young=young+runif(length(young),-0.01,0.01) #+unif to break the ties...
#don't use rnorm... might end up with more extreme values
old=old+runif(length(old),-0.01,0.01)

#data_std=(c(old,young)-min(c(old,young))+0.1)/(max(c(old,young))+1.2)

#data=cbind(data_std,c(data_std[1:47],data_std[-(1:47)]^0.5),c(rep(1,times=length(old)),rep(2,times=length(young))))

data=cbind(c(old,young),c(old,young^0.5),c(rep(1,times=length(old)),rep(2,times=length(young))))

sort((1/old)/sum(1/old),index.return=TRUE)
sort((1/young)/sum(1/young),index.return=TRUE)

####II.5 restrict to >8+standardize within sample?





library(Hmisc)

filepath1 = "data/match_condition_data/west_malefemale_35/2018.csv"
filepath2 = "data/match_condition_data/northeast+midwest_malefemale_35/2018.csv"
bac_dat1=read.csv(filepath1,header=TRUE)
age1=bac_dat1$AGE
alc_res1=bac_dat1$ALC_RES/1000
alc_res1=round(alc_res1,2)

young1=alc_res1[age1<35]
old1=alc_res1[age1>=35]

bac_dat2=read.csv(filepath2,header=TRUE)
age2=bac_dat2$AGE
alc_res2=bac_dat2$ALC_RES/1000
alc_res2=round(alc_res2,2)

young2=alc_res2[age2<35]
old2=alc_res2[age2>=35]


label_size=1.8*1.8
tcl_size=0.3*1.8
axis_size=1.5*1.6
position_lab=3.7*1.8
bot_left_mar=6*1.8
lwd=1.5*1.8
right_mar=0.5 #added 20150127
#####

pdf("Figure_west_ecdf.pdf",width=10.4,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,right_mar),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)

#ecdf:
young_ecdf1=ecdf(young1)
old_ecdf1=ecdf(old1)

young_ecdf2=ecdf(young2)
old_ecdf2=ecdf(old2)

plot(old_ecdf1, 
     verticals = TRUE, do.points = FALSE,main="",xlim=c(0,50/100),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
lines(young_ecdf1, verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,50/100,by=10/100),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC (g/dL)",cex.lab=label_size) #added 20150124
title(ylab="ecdf",cex.lab=label_size,line=position_lab)

dev.off()

pdf("Figure_northeast+midwest_NPMLE.pdf",width=10.4,height=10)
#setEPS()
#postscript("paper1_fig3.eps", height=10, width=5.2) #for SJS
par(mfrow=c(1,1),mar=c(bot_left_mar,bot_left_mar,1,right_mar),oma=c(0,0,0,1),las=1) #mar: bottom, left.. affect title(...line=?)!
require(graphics)


#weighted ecdfs
young = young2
old = old2
young_wtcdf=wtd.Ecdf(sort(young), weights=1/sort(young^0.5)/(sum(1/sort(young^0.5))), 
                     #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
                     normwt=FALSE, na.rm=TRUE)
old_wtcdf=wtd.Ecdf(sort(old), weights=1/sort(old)/(sum(1/sort(old))), 
                   #type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
                   normwt=FALSE, na.rm=TRUE)

plot(stepfun(old_wtcdf$x[-1],old_wtcdf$ecdf), 
     verticals = TRUE, do.points = FALSE,main="",xlim=c(0,50/100),ylim=c(0,1),lty=2,lwd=lwd,cex.axis=axis_size,xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",las=1,cex=0.7,tcl=tcl_size,frame.plot=FALSE)
abline(h=1,lty=2,col="gray",lwd=lwd)
lines(stepfun(young_wtcdf$x[-1],young_wtcdf$ecdf), verticals = TRUE, do.points = FALSE,lwd=lwd,lty=1)
axis(1, at=seq(0,50/100,by=10/100),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
axis(2, at=seq(0,1,by=0.2),las=1,tcl=tcl_size,cex.axis=axis_size) #side=1: bottom
title(xlab="BAC (g/dL)",cex.lab=label_size)
title(ylab=expression(paste(hat(F))),cex.lab=label_size,line=position_lab)


dev.off()





for(i in c(filepath1,filepath2)){
  bac_dat=read.csv(i,header=TRUE)
  age=bac_dat$AGE
  alc_res=bac_dat$ALC_RES/1000
  alc_res=round(alc_res,2)
  
  young=alc_res[age<35]
  old=alc_res[age>=35]
  data = list(young,old)
  r=c(0.5,1)
  nboot=1000
  group1 = young
  group2 = old
  
  young_ecdf=ecdf(young)
  old_ecdf=ecdf(old)
  plot(old_ecdf,verticals = TRUE, do.points = FALSE,ylim=c(0,1),main=i)
  lines(young_ecdf, verticals = TRUE, do.points = FALSE,col="red")
}

bac_dat=read.csv(filepath1,header=TRUE)
age=bac_dat$AGE
alc_res=bac_dat$ALC_RES/1000
alc_res=round(alc_res,2)

young=alc_res[age<35]
old=alc_res[age>=35]
data = list(young,old)
r=c(0.5,1)
nboot=1000
group1 = young
group2 = old
main = "west"
young_ecdf=ecdf(young)
old_ecdf=ecdf(old)

plot(old_ecdf,verticals = TRUE, do.points = FALSE,xlim=c(0,40/100),ylim=c(0,1),main=main)
lines(young_ecdf, verticals = TRUE, do.points = FALSE,col="red")



plotdata <- function(year,group1,group2){
  group1 = young
  group2 = old
  wtcdf1=wtd.Ecdf(sort(group1), weights=1/sort(group1^0.5)/(sum(1/sort(group1^0.5))), 
                  normwt=FALSE, na.rm=TRUE)
  wtcdf2=wtd.Ecdf(sort(group2), weights=1/sort(group2)/(sum(1/sort(group2))), 
                  normwt=FALSE, na.rm=TRUE)
  r = crosscheck(wtcdf1,wtcdf2)
  main = ifelse(r==T,paste0(as.character(year)," *cross"),as.character(year))
  
  plot(stepfun(wtcdf1$x[-1],wtcdf1$ecdf),verticals = T,main=main,
       do.points = F, xlim=c(0,max(wtcdf1$x,wtcdf2$x)),ylab="")
  lines(stepfun(wtcdf2$x[-1],wtcdf2$ecdf), verticals = T, do.points = FALSE, col="red")
}


crosscheck <- function(c1,c2){
  #length cdf1 < cdf2
  result = c()
  if(length(c1$x)<length(c2$x)){
    cdf1 = c1
    cdf2 = c2
  }else{
    cdf1 = c2
    cdf2 = c1
  }
  for(i in c(2:length(cdf1$x))){
    if(i == length(cdf1$x)){
      break
    }else{
      a = cdf1$ecdf[i] - cdf2$ecdf[i]
      b = cdf1$ecdf[i+1] - cdf2$ecdf[i+1]
      if(a*b<0 || a*b==0){
        result = c(result,cdf1$x[i])
      }
    }
    
  }
  if(length(result)==0){
    return(F)
  }else{
    return(T)
  }
}
BS_2sample_revised <- function(data, r, nboot,year,young,old){
  data = data
  T1 = data[[1]]
  T2 = data[[2]]
  T12 <- sort(unique(c(T1,T2)))
  
  nboot=nboot
  n1 = length(T1)
  n2 = length(T2)
  n = n1 + n2 
  nn = length(T12)
  
  weights1=T1 ^ r[1]
  weights2=T2 ^ r[2] 
  weights = c(weights1,weights2)
  WWeights1 = n1/(sum(1/weights1)) 
  WWeights2 = n2/(sum(1/weights2)) 
  WWeights  = c(WWeights1,WWeights2)
  
  kk = c(n1/n,n2/n)
  
  p_hat1 = WWeights1/(n1*weights1)
  p_hat2 = WWeights2/(n2*weights2)
  
  ECDF1 <- apply((matrix(rep(p_hat1,nn-1), nn-1, n1, byrow=TRUE) * (matrix(rep(T1,nn-1), nn-1, n1, byrow=TRUE) <= matrix(rep(T12[-nn],n1), nn-1, n1, byrow=FALSE))), 1, sum) #F_hat1
  ECDF2 <- apply((matrix(rep(p_hat2,nn-1), nn-1, n2, byrow=TRUE) * (matrix(rep(T2,nn-1), nn-1, n2, byrow=TRUE) <= matrix(rep(T12[-nn],n2), nn-1, n2, byrow=FALSE))), 1, sum) #F_hat2
  
  H_hat =  ( n1*ECDF1 + n2*ECDF2 )/n
  
  dH_hatvec = c(0, H_hat)
  
  theta_hat_t_t = theta_func("t","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)
  
  AD6_o<- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat*(1-H_hat)) * diff(dH_hatvec))
  AD6_t<- n * sum( ((ECDF1 - ECDF2)^2) / (theta_hat_t_t) * diff(dH_hatvec))
  
  #AD7
  theta_hat_s_t = theta_func("s","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)
  #theta_hat_s_t1 = theta_func1("s","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)
  
  b1  <- ECDF1 - ECDF2 #A2 first matrix[1,1]
  b2  <- b1 #A2 first matrix[1,2]
  a11_t <- theta_hat_t_t #A2 second matrix[1,1]
  a11_o <- H_hat*(1-H_hat) 
  a21_t <- theta_hat_s_t #A2 second matrix[2,1]&[1,2]
  a21_o <- H_hat%*%t(1-H_hat) 
  a22_t <- a11_t #A2 second matrix[2,2]
  a22_o <- a11_o #A2 second matrix[2,2]
  
  AD7_unit_t <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
  AD7_unit_o <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
  for(k in 2:(nn-1)){  #t k=2 nn=5, 2~4
    for(j in 1:(k-1)){ #s j=1
      #k = 2 
      #j = 1
      
      Vinv_t <- solve(matrix(c(a11_t[j], a21_t[j,k], a21_t[j,k], a22_t[k]),2,2)) #A2 second matrix inv
      Vinv_o <- solve(matrix(c(a11_o[j], a21_o[j,k], a21_o[j,k], a22_o[k]),2,2)) 
      
      dH_hat_t_vec <- c(0, H_hat[k])
      dH_hat_s_vec <- c(0, H_hat[j])
      
      AD7_unit_t[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv_t%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
      AD7_unit_o[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv_o%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
      
    }
  }
  
  AD7_o <- (n1*n2/n)*sum(AD7_unit_o, na.rm=TRUE)  # NA bc for j=1, AD7_unit[1, 1] not defined
  AD7_t <- n*sum(AD7_unit_t, na.rm=TRUE)  
  
  AD7p_t <- AD6_t + AD7_t
  AD7p_o <- AD6_o + AD7_o
  
  
  AD6_star_t <- rep(NA,nboot)
  AD7_star_t <- rep(NA,nboot)
  AD6_star_o <- rep(NA,nboot)
  AD7_star_o <- rep(NA,nboot)
  
  #b = 1 #for test
  for(b in 1:nboot){
    
    xi1<- rnorm(n1, 0, 1)
    xi2<- rnorm(n2, 0, 1)
    
    D1_star <- apply(matrix(rep((xi1 * p_hat1),nn-1), nn-1, n1, byrow=TRUE) * ((matrix(rep(T1,nn-1), nn-1, n1, byrow=TRUE) <= matrix(rep(T12[-nn],n1), nn-1, n1, byrow=FALSE)) - H_hat[-nn]), 1, sum)
    D2_star <- apply(matrix(rep((xi2 * p_hat2),nn-1), nn-1, n2, byrow=TRUE) * ((matrix(rep(T2,nn-1), nn-1, n2, byrow=TRUE) <= matrix(rep(T12[-nn],n2), nn-1, n2, byrow=FALSE)) - H_hat[-nn]), 1, sum)
    
    #AD6
    # AD_s_unit_6_I1 <- apply(matrix(rep(xi1,n-1), n-1, n1, byrow=TRUE) * ((matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE)) - H_hat[-n]), 1, mean)
    # AD_s_unit_6_I2 <- apply(matrix(rep(xi2,n-1), n-1, n2, byrow=TRUE) * ((matrix(rep(T2,n-1), n-1, n2, byrow=TRUE) <= matrix(rep(T12[-n],n2), n-1, n2, byrow=FALSE)) - H_hat[-n]), 1, mean)
    
    AD_s_unit_6 <- (D1_star - D2_star) ^ 2
    AD6_star_o[b] <- (n1*n2/n)*sum(AD_s_unit_6/(H_hat[-nn]*(1-H_hat[-nn]))*diff(dH_hatvec), na.rm = TRUE)
    AD6_star_t[b] <- n*sum(AD_s_unit_6/(theta_hat_t_t)*diff(dH_hatvec), na.rm = TRUE)
    
    
    #AD7
    AD7_unit_star_t <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
    AD7_unit_star_o <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
    
    # b1s_I1 <- apply(matrix(rep(xi1,n-1), n-1, n1, byrow=TRUE) * ((matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE)) - H_hat[-n]), 1, mean)
    # b1s_I2 <- apply(matrix(rep(xi2,n-1), n-1, n2, byrow=TRUE) * ((matrix(rep(T2,n-1), n-1, n2, byrow=TRUE) <= matrix(rep(T12[-n],n2), n-1, n2, byrow=FALSE)) - H_hat[-n]), 1, mean)
    b1s <- (D1_star - D2_star)
    
    b2s<- b1s # v2
    
    for(k in 2:(nn-1)){  #t
      for(j in 1:(k-1)){ #s
        dH_hat_t_vec = c(0, H_hat[k])
        dH_hat_s_vec = c(0, H_hat[j])
        
        Vinv_t <- solve(matrix(c(a11_t[j], a21_t[j,k], a21_t[j,k], a22_t[k]),2,2))
        Vinv_o <- solve(matrix(c(a11_o[j], a21_o[j,k], a21_o[j,k], a22_o[k]),2,2)) 
        
        AD7_unit_star_t[j,k]<- (t(c(b1s[j], b2s[k]))%*%Vinv_t%*%c(b1s[j], b2s[k]))*diff(dH_hat_s_vec)*diff(dH_hat_t_vec)
        AD7_unit_star_o[j,k]<- (t(c(b1s[j], b2s[k]))%*%Vinv_o%*%c(b1s[j], b2s[k]))*diff(dH_hat_s_vec)*diff(dH_hat_t_vec)
        
      }
    }
    AD7_star_t[b]<- n*sum(AD7_unit_star_t,na.rm=TRUE)
    AD7_star_o[b]<- (n1*n2/n)*sum(AD7_unit_star_o,na.rm=TRUE)
    
    
  }
  AD7p_star_t  <- AD6_star_t + AD7_star_t
  rej.rate6_t <- AD6_t  > quantile(AD6_star_t,0.95)
  quantile(AD6_star_t,0.95,type=1)
  rej.rate7_t <- AD7_t  > quantile(AD7_star_t,0.95)
  rej.rate7p_t <- AD7p_t > quantile(AD7p_star_t,0.95)
  
  AD7p_star_o  <- AD6_star_o + AD7_star_o
  rej.rate6_o <- AD6_o  > quantile(AD6_star_o,0.95)
  rej.rate7_o <- AD7_o  > quantile(AD7_star_o,0.95)
  rej.rate7p_o <- AD7p_o > quantile(AD7p_star_o,0.95) 
  
  p6_t <- mean(AD6_star_t >= AD6_t,na.rm=TRUE)
  p7_t <- mean(AD7_star_t >= AD7_t,na.rm=TRUE)
  p7p_t <- mean(AD7p_star_t >= AD7p_t, na.rm=TRUE)
  
  p6_o <- mean(AD6_star_o >= AD6_o,na.rm=TRUE)
  p7_o <- mean(AD7_star_o >= AD7_o,na.rm=TRUE)
  p7p_o <- mean(AD7p_star_o >= AD7p_o, na.rm=TRUE)
  
  out_array=list(
    rej.rate6_o=rej.rate6_o,
    rej.rate7_o=rej.rate7_o,
    rej.rate7p_o=rej.rate7p_o,
    rej.rate6_t=rej.rate6_t,
    rej.rate7_t=rej.rate7_t,
    rej.rate7p_t=rej.rate7p_t
  )
  p = list(
    p6_o = p6_o,
    p7_o = p7_o,
    p7p_o = p7p_o,
    p6_t = p6_t,
    p7_t = p7_t,
    p7p_t = p7p_t
  )
  if(p$p7_t <= p$p6_t){
    plotdata(year,young,old)
    print(c(lapply(out_array, function(i)ifelse(i==T,"reject","not to reject")),
            p))
  }
}

BS_2sample_revised(data,r,nboot,2009,young,old)
'''
filepath = "data/match_condition_data/midwest+west_male_30/2009.csv"
bac_dat=read.csv(filepath,header=TRUE)
age=bac_dat$AGE
i=2009
a=30
index = 1
if(i<2015){
  alc_res=bac_dat$ALC_RES/100  
}else{
  alc_res=bac_dat$ALC_RES/1000
  alc_res=round(alc_res,2)
}
young=alc_res[age<a]
old=alc_res[age>=a]
data = list(young,old)
seed = index*10000 + i
set.seed(seed)
BS_2sample_revised(data,r,nboot)
plotdata(2009,young,old)
'''
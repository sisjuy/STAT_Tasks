library(Hmisc)

filepath = "data/match_condition_data/midwest+west_male_30/2009.csv"



bac_dat=read.csv(filepath,header=TRUE)
age=bac_dat$AGE
alc_res=bac_dat$ALC_RES/100  

young=alc_res[age<30]
old=alc_res[age>=30]
data = list(young,old)
r=c(0.5,1)
nboot=1000

plotdata <- function(year,group1,group2){
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
  #data = data
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
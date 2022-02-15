
library(tidyverse)

# datafilter
#for 2009~2014
datafilter <- function(accident_path, person_path, output, state,sex){
  beforedata1  = read_csv(accident_path, col_names=T, skip = 0)
  beforedata1 = filter(beforedata1, STATE == state) #California(CA)
  
  beforedata2  = read_csv(person_path, col_names=T, skip = 0)
  beforedata2 = filter(beforedata2, STATE == state) #California(CA)
  
  step2a = filter(beforedata1, ROAD_FNC == 1) #Rural Principal Arterial-Interstate
  step2b = filter(beforedata1, ROAD_FNC == 11) #Urban Principal Arterial-Interstate
  step2 = rbind(step2a,step2b)
  step2 = step2[order(step2[,44]),] #order by state-case number(ST_CASE)
  
  step3 = data.frame() #combine two tables by state-case number(ST_CASE)
  for (i in c(1:nrow(step2))){
    for (j in c(1:nrow(beforedata2))){
      if (beforedata2$ST_CASE[j] == step2$ST_CASE[i]){
        tempdata = beforedata2[j,]
        step3 = rbind(step3, tempdata)
      }
      if (beforedata2$ST_CASE[j] > step2$ST_CASE[i]){
        break
      }
    }
  }
  
  step4 = filter(step3, SEX == sex) #sex = male
  step5 = filter(step4, ATST_TYP == 1) #alcohol test = Blood
  step6 = filter(step5, ALC_RES >= 1, ALC_RES < 96) #BAC>=0.01%
  
  
  
  write_csv(step6,file=output)
  
  
}
#for 2015~
datafilter1 <- function(accident_path, person_path, output, state,sex){
  beforedata1  = read_csv(accident_path, col_names=T, skip = 0)
  beforedata1 = filter(beforedata1, STATE == state) #California(CA)
  
  beforedata2  = read_csv(person_path, col_names=T, skip = 0)
  beforedata2 = filter(beforedata2, STATE == state) #California(CA)
  
  step2a = filter(beforedata1, RUR_URB == 1) #Rural Principal Arterial-Interstate
  step2b = filter(beforedata1, RUR_URB == 2) #Urban Principal Arterial-Interstate
  step2 = rbind(step2a,step2b)
  step2 = filter(step2, FUNC_SYS == 1)
  step2 = step2[order(step2$ST_CASE),]
  
  step3 = data.frame() #combine two tables by state-case number(ST_CASE)
  for (i in c(1:nrow(step2))){
    for (j in c(1:nrow(beforedata2))){
      if (beforedata2$ST_CASE[j] == step2$ST_CASE[i]){
        tempdata = beforedata2[j,]
        step3 = rbind(step3, tempdata)
      }
      if (beforedata2$ST_CASE[j] > step2$ST_CASE[i]){
        break
      }
    }
  }
  
  step4 = filter(step3, SEX == sex) #sex = male
  step5 = filter(step4, ATST_TYP == 1) #alcohol test = Blood
  step6 = filter(step5, ALC_RES >= 10, ALC_RES < 960) #BAC>=0.01%
  
  
  
  write_csv(step6,file=output)
  
  
}

#dir create state and sex
dir.create(paste0("data/","STATE/",as.character(state),"/",
                  as.character(ifelse(sex==1,"male","female"))))

for(i in c(2015:2019)){
  accident_path = paste0("data/",as.character(i),"/ACCIDENT.csv")
  person_path = paste0("data/",as.character(i),"/PERSON.csv")
  output = paste0("data/","STATE/",as.character(state),
                  "/",as.character(ifelse(sex==1,"male","female")),
                  "/",as.character(i),"data.csv")
  
  datafilter(accident_path,person_path,output,state,sex)
}


filterdata <- function(state,sex){
  for(i in c(2015:2019)){
    accident_path = paste0("data/",as.character(i),"/ACCIDENT.csv")
    person_path = paste0("data/",as.character(i),"/PERSON.csv")
    output = paste0("data/","STATE/",as.character(state),
                    "/",as.character(ifelse(sex==1,"male","female")),
                    "/",as.character(i),"data.csv")
    
    datafilter1(accident_path,person_path,output,state,sex)
  }
}

for(ii in x[-c(1:9)]){
  #state 1:56
  #dir.create(paste0("data/","STATE/",as.character(ii)))
  for(jj in c(1:2)){
    dir.create(paste0("data/","STATE/",as.character(ii),"/",
                      as.character(ifelse(jj==1,"male","female"))))
    filterdata(ii,jj)
  }
}
x[-(1:5)]
#unique STATE
x = {}
for(i in c(2015:2019)){
  accident_path = paste0("data/",as.character(i),"/ACCIDENT.csv")
  person_path = paste0("data/",as.character(i),"/PERSON.csv")
  d = read_csv(accident_path, col_names=T, skip = 0)
  x = unique(d$STATE)
}
x

# plotdata
library(Hmisc)
plotdata <- function(year,group1,group2){
  
  
  wtcdf1=wtd.Ecdf(sort(group1), weights=1/sort(group1^0.5)/(sum(1/sort(group1^0.5))), 
                  normwt=FALSE, na.rm=TRUE)
  wtcdf2=wtd.Ecdf(sort(group2), weights=1/sort(group2)/(sum(1/sort(group2))), 
                  normwt=FALSE, na.rm=TRUE)
  
  plot(stepfun(wtcdf1$x[-1],wtcdf1$ecdf),verticals = T,main=as.character(year),
         do.points = F, xlim=c(0,max(wtcdf1$x,wtcdf2$x)),ylab="")
  lines(stepfun(wtcdf2$x[-1],wtcdf2$ecdf), verticals = T, do.points = FALSE, col="red")
  
  
}

for(i in c(2009:2014)){
  
  filepath = paste0("data/",as.character(i),"data.csv")
  bac_dat=read.csv(filepath,header=TRUE)
  age=bac_dat$AGE
  alc_res=bac_dat$ALC_RES/100 
  young=alc_res[age<30]
  old=alc_res[age>=30]
  data = list(young,old)
  plotdata(i,young,old)
  cat("year:",i,"\n")
  BS_2sample(data,c(0.5,1),1000)
}

# hypothesis
BS_2sample <- function(data, r, nboot){
  T1 = data[[1]]
  T2 = data[[2]]
  #T12 = sort(c(T1,T2))
  T12 = sort(unique(c(T1,T2)))
  
  nboot=nboot
  n1 = length(T1)
  n2 = length(T2)
  n = n1 + n2 
  nn = length(T12)
  
  weights1=T1 ^ r[1]
  weights2=T2 ^ r[2]
  
  WWeights1 = n1/(sum(1/weights1))
  WWeights2 = n2/(sum(1/weights2))
  
  p_hat1 = WWeights1/(n1*weights1)
  p_hat2 = WWeights2/(n2*weights2)
  
  #F hat
  ECDF1 <- apply((matrix(rep(p_hat1,nn-1), nn-1, n1, byrow=TRUE) * (matrix(rep(T1,nn-1), nn-1, n1, byrow=TRUE) <= matrix(rep(T12[-nn],n1), nn-1, n1, byrow=FALSE))), 1, sum) #F_hat1
  ECDF2 <- apply((matrix(rep(p_hat2,nn-1), nn-1, n2, byrow=TRUE) * (matrix(rep(T2,nn-1), nn-1, n2, byrow=TRUE) <= matrix(rep(T12[-nn],n2), nn-1, n2, byrow=FALSE))), 1, sum) #F_hat2
  
  H_hat =  ( n1*ECDF1 + n2*ECDF2 )/n
  
  dH_hatvec = c(0, H_hat)
  AD6<- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat*(1-H_hat)) * diff(dH_hatvec))
  
  b1 <- ECDF1 - ECDF2      
  b2 <- b1                 
  a11 <- H_hat*(1-H_hat)   
  a21 <- H_hat%*%t(1-H_hat)
  a22 <- a11
  AD7_unit <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
  for(k in 2:(nn-1)){ #t k=2
    for(j in 1:(k-1)){ #s j=1
      Vinv <- solve(matrix(c(a11[j], a21[j,k], a21[j,k], a22[k]),2,2)) 
      dH_hat_t_vec <- c(0, H_hat[k])
      dH_hat_s_vec <- c(0, H_hat[j])
      
      AD7_unit[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    }
  }
  AD7 <- (n1*n2/n)*sum(AD7_unit, na.rm=TRUE)  
  AD7p <- AD6 + AD7
  
  AD6_star <- rep(NA,nboot)
  AD7_star <- rep(NA,nboot)
  for(b in 1:nboot){
    xi1<- rnorm(n1, 0, 1)
    xi2<- rnorm(n2, 0, 1)
    
    D1_star <- apply(matrix(rep((xi1 * p_hat1),nn-1), nn-1, n1, byrow=TRUE) * ((matrix(rep(T1,nn-1), nn-1, n1, byrow=TRUE) <= matrix(rep(T12[-nn],n1), nn-1, n1, byrow=FALSE)) - H_hat[-nn]), 1, sum)
    D2_star <- apply(matrix(rep((xi2 * p_hat2),nn-1), nn-1, n2, byrow=TRUE) * ((matrix(rep(T2,nn-1), nn-1, n2, byrow=TRUE) <= matrix(rep(T12[-nn],n2), nn-1, n2, byrow=FALSE)) - H_hat[-nn]), 1, sum)
    
    AD_s_unit_6 <- (D1_star - D2_star) ^ 2
    AD6_star[b] <- (n1*n2/n)*sum(AD_s_unit_6/(H_hat[-nn]*(1-H_hat[-nn]))*diff(dH_hatvec), na.rm = TRUE)
    
    AD7_unit_star <- matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1)
    
    b1s <- (D1_star - D2_star)
    
    b2s<- b1s # v2
    
    for(k in 2:(nn-1)){   #t
      for(j in 1:(k-1)){ #s
        dH_hat_t_vec = c(0, H_hat[k])
        dH_hat_s_vec = c(0, H_hat[j])
        
        Vinv<- solve(matrix(c(a11[j], a21[j,k], a21[j,k], a22[k]),2,2))
        AD7_unit_star[j,k]<- (t(c(b1s[j], b2s[k]))%*%Vinv%*%c(b1s[j], b2s[k]))*diff(dH_hat_s_vec)*diff(dH_hat_t_vec)
      }
    }
    AD7_star[b]<- (n1*n2/n)*sum(AD7_unit_star,na.rm=TRUE)
    
    
  }
  AD7p_star  <- AD6_star + AD7_star
  
  rej.rate6<- AD6  > quantile(AD6_star,0.95,na.rm=TRUE)
  rej.rate7<- AD7  > quantile(AD7_star,0.95,na.rm=TRUE)
  rej.rate7p<- AD7p > quantile(AD7p_star,0.95,na.rm=TRUE)
  l1 = list(rej.rate6,rej.rate7,rej.rate7p)
  l2 = {}
  for(i in c(1:3)){
    if(l1[i]==TRUE){
      l2[i]="reject"
    }else{
      l2[i]="not to reject"
    }
  }
  cat("AD6: ",l2[1],"\n","AD7",l2[2],"\n","AD7p",l2[3],"\n")
  #return(list(AD6=l2[1],AD7=l2[2],AD7p=l2[3]))
}

filepath = paste0("data/Region/South/male/",
                  as.character(2016),".csv")
bac_dat=read.csv(filepath,header=TRUE)
age=bac_dat$AGE
alc_res=bac_dat$ALC_RES/1000
round(alc_res,2)

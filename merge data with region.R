Northeast = c(9,23,25,33,34,36,42,44,50)
Midwest = c(17,18,19,20,26,27,29,31,38,39,46,55)
South = c(1,5,10,12,13,21,22,24,28,37,40,45,47,48,51,54)
West = c(2,4,6,8,15,16,30,32,35,41,49,53,56)
US = list(Northeast = Northeast,Midwest=Midwest,South=South,West=West)


for(region in US){
  for(state in region){
    for(sex in c(1,2)){
      filename = paste0("data/STATE",as.character(state),"/",
                        as.character(ifelse(sex==1,"male","female")),
                        "/",as.character(year),"data.csv")
      d = read_csv(filename,col_names=T, skip = 0)
      
    }
    
   
  }
}





merge <- function(state,sex){
  
}


library(tidyverse)
#male

for(sex in c(1,2)){
  #dir.create(paste0("data/Region/",names(US)[4],"/",
             #as.character(ifelse(sex==1,"male","female"))))
  for(year in c(2015:2019)){
    temp = {}
    for(state in US$West){
      filename = paste0("data/STATE/",as.character(state),"/",
                        as.character(ifelse(sex==1,"male","female")),
                        "/",as.character(year),"data.csv")
      d = read_csv(filename,col_names=T, skip = 0)
      temp = rbind(temp,d)
    }
    output = paste0("data/Region/",names(US)[4],"/",
                    as.character(ifelse(sex==1,"male","female")),"/",
                    as.character(year),".csv")
    write_csv(temp, file=output)
  }
}

#find nrow = 260~500
#northeast
names(US)
for(year in c(2009:2019)){
  sex = 2
  filepath = paste0("data/Region/",names(US)[4],"/",
                    as.character(ifelse(sex==1,"male","female")),"/",
                    as.character(year),".csv")
  bac_dat=read.csv(filepath,header=TRUE)
  age=bac_dat$AGE
  alc_res=bac_dat$ALC_RES/1000 
  young=alc_res[age<30]
  old=alc_res[age>=30]
  if(length(young)>=130 && length(old)>=130){
    print(year)
  }
}


#file list test
filename = list.files(path="data/match_condition_data")
print(filename)
yearlist = list.files(path=paste0("data/match_condition_data/",filename[1]))
yearlist
year = strsplit(yearlist,".csv")
year
i = year[[1]]
library(Hmisc)
r=c(0.5,1)
nboot=1000
plotdata <- function(year,group1,group2){
  
  
  wtcdf1=wtd.Ecdf(sort(group1), weights=1/sort(group1^0.5)/(sum(1/sort(group1^0.5))), 
                  normwt=FALSE, na.rm=TRUE)
  wtcdf2=wtd.Ecdf(sort(group2), weights=1/sort(group2)/(sum(1/sort(group2))), 
                  normwt=FALSE, na.rm=TRUE)
  
  plot(stepfun(wtcdf1$x[-1],wtcdf1$ecdf),verticals = T,main=as.character(year),
       do.points = F, xlim=c(0,max(wtcdf1$x,wtcdf2$x)),ylab="")
  lines(stepfun(wtcdf2$x[-1],wtcdf2$ecdf), verticals = T, do.points = FALSE, col="red")
  
}
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
  
  p6 <- mean(AD6_star > AD6,na.rm=TRUE)
  p7 <- mean(AD7_star > AD7,na.rm=TRUE)
  p7p <- mean(AD7p_star > AD7p, na.rm=TRUE)
  
  l1 = list(rej.rate6,rej.rate7,rej.rate7p)
  star = list(AD6_star,AD7_star,AD7p_star)
  
  l2 = {}
  for(i in c(1:3)){
    if(l1[i]==TRUE){
      if(length(which(is.na(star[i])))!=0){
        l2[i]="reject**"
      }else{
        l2[i]="reject"
      }
      
    }else{
      if(length(which(is.na(star[i])))!=0){
        l2[i]="not to reject**"
      }else{
        l2[i]="not to reject"
      }
    }
  }
  cat("** represent there have na in star")
  cat("AD6:",l2[1],"\n",
      "AD7:",l2[2],"\n",
      "AD7p:",l2[3],"\n",
      "p6:",p6,"\n",
      "p7:",p7,"\n",
      "p7p:",p7p,"\n",
      "pcompare:",ifelse(p7<p6,"p7<p6"," "))
}
out <- function(a){
  for(i in year){
    i = as.integer(i)
    filepath = paste0("data/match_condition_data/",
                      filename[1],"/",as.character(i),".csv")
    
    
    #filepath = paste0("data/Region/Midwest/malefemale/",
    #                  as.character(i),".csv")
    #filepath1 = paste0("data/Region/West/malefemale/",
    #                  as.character(i),".csv")
    bac_dat=read.csv(filepath,header=TRUE)
    #bac_dat1=read.csv(filepath1,header=TRUE)
    age=bac_dat$AGE
    #age1=bac_dat1$AGE
    if(i<2015){
      alc_res=bac_dat$ALC_RES/100  
      #alc_res1=bac_dat1$ALC_RES/100
    }else{
      alc_res=bac_dat$ALC_RES/1000
      alc_res=round(alc_res,2)
      #alc_res1=bac_dat1$ALC_RES/1000
      #alc_res1=round(alc_res1,2)
    }
    young=alc_res[age<a]
    old=alc_res[age>=a]
    #young1=alc_res1[age1<a]
    #old1=alc_res1[age1>=a]
    
    #young=c(young,young1)
    #old=c(old,old1)
    data = list(young,old)
    plotdata(i,young,old)
    cat("len young: ", length(young),"\n")
    cat("len old: ", length(old),"\n")
    BS_2sample(data,r,nboot)
  }
}
out(30)

  
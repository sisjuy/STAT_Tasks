library(tidyverse)
#dir create combine gender
#Northeast = Northeast,Midwest=Midwest,South=South,West=West
dir.create(paste0("data/Region/",names(US)[4],"/malefemale"))
US$Northeast
temp = {}
for(year in c(2009:2019)){
  temp={}
  for(sex in c(1:2)){
    for(state in US$West){
      filename = paste0("data/STATE/",as.character(state),"/",
                        as.character(ifelse(sex==1,"male","female")),
                        "/",as.character(year),"data.csv")
      d = read_csv(filename,col_names=T, skip = 0)
      temp = rbind(temp,d)
    }
  }
  output = paste0("data/Region/",names(US)[4],"/malefemale/",
                  as.character(year),".csv")
  write_csv(temp, file=output)
}

#check sample size
for(year in c(2009:2019)){
  
  filepath = paste0("data/Region/",names(US)[4],"/malefemale/",
                    
                    as.character(year),".csv")
  bac_dat=read.csv(filepath,header=TRUE)
  age=bac_dat$AGE
  if (year<2015){
    alc_res = bac_dat$ALC_RES/100
  }else{
    alc_res=bac_dat$ALC_RES/1000
  }
  young=alc_res[age<40]
  old=alc_res[age>=40]
  if(length(young)>=130 && length(old)>=130 && length(young)<=250 && length(old)<=250){
    print(year)
  }
}

#combine region
x = seq(from=2009,to=2019)
a=30
for(i in x){
  
  filepath1 = paste0("data/Region/",names(US)[1],"/malefemale/",
                      as.character(i),".csv")
  filepath2 = paste0("data/Region/",names(US)[2],"/malefemale/",
                     as.character(i),".csv")
  filepath3 = paste0("data/Region/",names(US)[3],"/malefemale/",
                     as.character(i),".csv")
  filepath4 = paste0("data/Region/",names(US)[4],"/malefemale/",
                     as.character(i),".csv")
  bac_dat1=read.csv(filepath1,header=TRUE)
  bac_dat2=read.csv(filepath2,header=TRUE)
  bac_dat3=read.csv(filepath3,header=TRUE)
  bac_dat4=read.csv(filepath4,header=TRUE)
  
  bac_dat=rbind(bac_dat2,bac_dat3,bac_dat4)
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
  ysum = length(young)
  osum = length(old)
  if(ysum>=130 && osum>=130 && ysum<=250 && osum<=250){
    print(i)
    print(length(young))
    print(length(old))
  }
  #cat("year: ",i,"\n")
  #cat("len young: ", length(young+young1),"\n")
  #cat("len old: ", length(old+old1),"\n")
  #cat(length(old1))
}
yearfun <- function(region,a){
  x={}
  for(year in c(2009:2019)){
    
    filepath = paste0("data/Region/",names(US)[region],"/malefemale/",
                      
                      as.character(year),".csv")
    bac_dat=read.csv(filepath,header=TRUE)
    age=bac_dat$AGE
    if (year<2015){
      alc_res = bac_dat$ALC_RES/100
    }else{
      alc_res=bac_dat$ALC_RES/1000
    }
    young=alc_res[age<a]
    old=alc_res[age>=a]
    if(length(young)>=130 && length(old)>=130 && length(young)<=250 && length(old)<=250){
      x=append(x,year)
    }
  }
  return(x)
}


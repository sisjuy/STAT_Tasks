x = seq(from=2009,to=2019)
Northeast = c(9,23,25,33,34,36,42,44,50)
Midwest = c(17,18,19,20,26,27,29,31,38,39,46,55)
South = c(1,5,10,12,13,21,22,24,28,37,40,45,47,48,51,54)
West = c(2,4,6,8,15,16,30,32,35,41,49,53,56)
US = list(Northeast = Northeast,Midwest=Midwest,South=South,West=West)
a=35
sex=3
combinetext = paste0(names(US)[1],"+",names(US)[4])
dir.create(paste0("data/Region/",combinetext))
dir.create(paste0("data/Region/",combinetext,"/",as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale")))))
for(i in x){
  
  filepath1 = paste0("data/Region/",names(US)[1],"/",
                    as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale"))),"/",
                    as.character(i),".csv")
  filepath2 = paste0("data/Region/",names(US)[2],"/",
                    as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale"))),"/",
                    as.character(i),".csv")
  filepath3 = paste0("data/Region/",names(US)[3],"/",
                    as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale"))),"/",
                    as.character(i),".csv")
  filepath4 = paste0("data/Region/",names(US)[4],"/",
                     as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale"))),"/",
                     as.character(i),".csv")
  bac_dat1=read.csv(filepath1,header=TRUE)
  bac_dat2=read.csv(filepath2,header=TRUE)
  bac_dat3=read.csv(filepath3,header=TRUE)
  bac_dat4=read.csv(filepath4,header=TRUE)
  
  bac_dat=rbind(bac_dat1,bac_dat4)
  
  output = paste0("data/Region/",combinetext,"/",
                  as.character(ifelse(sex==1,"male",ifelse(sex==2,"female","malefemale"))),
                  "/",
                  as.character(i),".csv")
  write.csv(bac_dat, file=output)
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
    
  }
  
  #cat("year: ",i,"\n")
  #cat("len young: ", length(young+young1),"\n")
  #cat("len old: ", length(old+old1),"\n")
  #cat(length(old1))
}


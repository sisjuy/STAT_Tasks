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
  dir.create(paste0("data/Region/",names(US)[4],"/",
             as.character(ifelse(sex==1,"male","female"))))
  for(year in c(2009:2019)){
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
  sex = 1
  filepath = paste0("data/Region/",names(US)[3],"/",
                    as.character(ifelse(sex==1,"male","female")),"/",
                    as.character(year),".csv")
  bac_dat=read.csv(filepath,header=TRUE)
  age=bac_dat$AGE
  alc_res=bac_dat$ALC_RES/100 
  young=alc_res[age<30]
  old=alc_res[age>=30]
  if(length(young)>=130 && length(old)>=130){
    print(year)
  }
}





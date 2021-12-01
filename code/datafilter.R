# restrictions
# sex = male
# state = California(CA)
# interstate highway accidents = Yes
# alcohol test = Blood

##
rm(list=ls(all=TRUE))

library(readr)
library(dplyr)
#library(sas7bdat)

beforedata1  = read_csv("data/2014/ACCIDENT.CSV", col_names=T, skip = 0)
beforedata1 = filter(beforedata1, STATE == 6) #California(CA)

beforedata2  = read_csv("data/2014/PERSON.CSV", col_names=T, skip = 0)
beforedata2 = filter(beforedata2, STATE == 6) #California(CA)

afterdata   = read_csv("data.csv", col_names=T, skip = 0)

#step1 = filter(beforedata1, ROUTE == 1) #interstate highway #without? OK 0802
step2a = filter(beforedata1, ROAD_FNC == 1) #Rural Principal Arterial-Interstate
step2b = filter(beforedata1, ROAD_FNC == 11) #Urban Principal Arterial-Interstate
step2 = rbind(step2a,step2b)
step2$ST_CASE
step2 = step2[order(step2$ST_CASE),] #order by state-case number(ST_CASE)

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

step4 = filter(step3, SEX == 1) #sex = male
step5 = filter(step4, ATST_TYP == 1) #alcohol test = Blood
step6 = filter(step5, ALC_RES >= 1, ALC_RES < 96) #BAC>=0.01%
#filter(step3, HARM_EV != 23) #appendex: if harmful event=23 then interstate highway should not =1
step6


step6$AGE == afterdata$AGE
step6$AGE
##same? TRUE
step6$AGE == afterdata$AGE
step6$ALC_RES == afterdata$ALC_RES

##
# for testing
# kkk = data.frame()
# for ( j in c(1:nrow(step6))){
#   for (i in c(1:nrow(step2))){
#     if (step2$ST_CASE[i] == step6$ST_CASE[j]){
#       kkk = rbind(kkk,step2[i,])
#     }
#   }
# }
## choose h_hat's  jump points and

fun <- function(data, r, nboot){
  data = data1
  T1 = data[[1]]
  T2 = data[[2]]
  T12 = sort(c(T1,T2))
  
  nboot=nboot
  n1 = length(T1)
  n2 = length(T2)
  n = n1 + n2 
  
  weights1=T1 ^ r[1]
  weights2=T2 ^ r[2]
  
  WWeights1 = n1/(sum(1/weights1))
  WWeights2 = n2/(sum(1/weights2))
  
  p_hat1 = WWeights1/(n1*weights1)
  p_hat2 = WWeights2/(n2*weights2)

  #F hat
  ECDF1 <- apply((matrix(rep(p_hat1,n-1), n-1, n1, byrow=TRUE) * (matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE))), 1, sum) #F_hat1
  ECDF2 <- apply((matrix(rep(p_hat2,n-1), n-1, n2, byrow=TRUE) * (matrix(rep(T2,n-1), n-1, n2, byrow=TRUE) <= matrix(rep(T12[-n],n2), n-1, n2, byrow=FALSE))), 1, sum) #F_hat2
  #matrix(rep(p_hat1,n-1), n-1, n1, byrow=TRUE)
  #matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE)
  
  #calculate AD6
  H_hat =  ( n1*ECDF1 + n2*ECDF2 )/n
  
  
  #store jump point in h_hat
  H_hat1={}
  i = 2
  while(i<=length(H_hat)){
    if((H_hat[i])!=(H_hat[i-1])){
      H_hat1 = append(H_hat1,H_hat[i-1])
    }
    if(i==length(H_hat) && H_hat[i]!=H_hat[i-1]){
      H_hat1 = append(H_hat1, H_hat[i])
    }
    i = i + 1 
  }
  
  
  dH_hatvec = c(0, H_hat)
  dH_hatvec1 = c(0,H_hat1)
  
  AD6<- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat*(1-H_hat)) * diff(dH_hatvec))
  #AD61 <- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat1*(1-H_hat1)) * diff(dH_hatvec1))
  AD6
  b1 <- ECDF1 - ECDF2 
  b2 <- b1      

  a11 <- H_hat*(1-H_hat)   
  a11_1 <- H_hat1*(1-H_hat1)
  
  
  a21 <- H_hat%*%t(1-H_hat)
  a21_1 <- H_hat1%*%t(1-H_hat1)
  
  a22 <- a11
  a22_1 <- a11_1
  
  AD7_unit <- matrix(rep(NA,(n-1)*(n-1)),n-1,n-1)
  
  #kj1={}
  #k1={}
  #j1={}
  nn = length(a11_1)+1

  for(k in 2:(nn-1)){ #t k=2
    for(j in 1:(k-1)){ #s j=1
      #if (det(matrix(c(a11[j], a21[j,k], a21[j,k], a22[k]),2,2))==0){
      #  cat("k:",k,"j:",j,"\n")
      #  kj1 = append(kj1,paste(k,j,sep=","))
      #  k1 = append(k1,k)
      #  j1 = append(j1,j)
      #}
      Vinv <- solve(matrix(c(a11_1[j], a21_1[j,k], a21_1[j,k], a22_1[k]),2,2)) #A2 second matrix inv
      dH_hat_t_vec <- c(0, H_hat[k])
      dH_hat_s_vec <- c(0, H_hat[j])
      
      AD7_unit[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    }
  }
  
  #observe a11 a21:
  #for(i in c(1:length(k1))){
  #  print(matrix(c(a11[j1[i]], a21[j1[i],k1[i]], a21[j1[i],k1[i]], a22[k1[i]]),2,2))
  #}
  
  
  
  
  
  
  AD7 <- (n1*n2/n)*sum(AD7_unit, na.rm=TRUE)  
  AD7p <- AD6 + AD7
  AD7
  #nboot=1000
  AD6_star <- rep(NA,nboot)
  AD7_star <- rep(NA,nboot)
  for(b in 1:nboot){
    xi1<- rnorm(n1, 0, 1)
    xi2<- rnorm(n2, 0, 1)
    
    D1_star <- apply(matrix(rep((xi1 * p_hat1),n-1), n-1, n1, byrow=TRUE) * ((matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE)) - H_hat[-n]), 1, sum)
    D2_star <- apply(matrix(rep((xi2 * p_hat2),n-1), n-1, n2, byrow=TRUE) * ((matrix(rep(T2,n-1), n-1, n2, byrow=TRUE) <= matrix(rep(T12[-n],n2), n-1, n2, byrow=FALSE)) - H_hat[-n]), 1, sum)
    #AD6_star
    AD_s_unit_6 <- (D1_star - D2_star) ^ 2
    AD6_star[b] <- (n1*n2/n)*sum(AD_s_unit_6/(H_hat[-n]*(1-H_hat[-n]))*diff(dH_hatvec), na.rm = TRUE)
    #AD7
    AD7_unit_star <- matrix(rep(NA,(n-1)*(n-1)),n-1,n-1)
    
    b1s <- (D1_star - D2_star)
    
    b2s<- b1s # v2
    
    for(k in 2:(nn-1)){   #t
      for(j in 1:(k-1)){ #s
        dH_hat_t_vec = c(0, H_hat1[k])
        dH_hat_s_vec = c(0, H_hat1[j])
        
        Vinv<- solve(matrix(c(a11_1[j], a21_1[j,k], a21_1[j,k], a22_1[k]),2,2))
        AD7_unit_star[j,k]<- (t(c(b1s[j], b2s[k]))%*%Vinv%*%c(b1s[j], b2s[k]))*diff(dH_hat_s_vec)*diff(dH_hat_t_vec)
      }
    }
    AD7_star[b]<- (n1*n2/n)*sum(AD7_unit_star,na.rm=TRUE)
    
    
  }
  AD7p_star  <- AD6_star + AD7_star
  
  rej.rate6<- AD6  > quantile(AD6_star,0.95,na.rm=TRUE)
  rej.rate7<- AD7  > quantile(AD7_star,0.95,na.rm=TRUE)
  rej.rate7p<- AD7p > quantile(AD7p_star,0.95,na.rm=TRUE)
  
  l1 = c(rej.rate6,rej.rate7,rej.rate7p)
  l2 = {}
  for(i in c(1:3)){
    if(l1[i]==TRUE){
      l2[i]="reject"
    }else{
      l2[i]="not to reject"
    }
  }
  
  list(AD6=l2[1],AD7=l2[2],AD7p=l2[3])
  return(list(AD6=l2[1],AD7=l2[2],AD7p=l2[3]))
  
}
bac_dat=read.csv("data.csv",header=TRUE)

age=bac_dat$AGE
alc_res=bac_dat$ALC_RES/100 
young=alc_res[age<30]
old=alc_res[age>=30]
T1 = young
T2 = old
data1 = list(T1,T2)
r = c(0.5,1)
nboot = 1000
fun(data1, r, nboot)


##########
#data2
T1 <- rbeta(90, 6.5, 7)
T2 <- rbeta(100, 3, 3)
T1
max_length <- max(c(length(T1), length(T2)))
max_length
data2 <- data.frame(col1 = c(T1,rep(NA, max_length - length(T1))),
                    col2 = c(T2,rep(NA, max_length - length(T2))))
data1
r = c(0.5,1)
nboot = 10

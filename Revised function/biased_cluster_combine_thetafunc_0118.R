## 20220106 theta function
## 20210805 biased
## 20210708 try cluster
## 20210705 for loop to matrix
## 20210624 bootstrap 1sample to 2sample
## 20210610 AD6, AD7 1sample to 2sample

## 20210420 updating paper_2t_simu_20210415's 1 sample to 2-sample

## 20180228 
## Codes for: An improved integral-typed statistic for testing crossing survival functions
## section 3.1 size simulation
## (1) Table 1
## (2) Figure 1,2,3

##(1)
## Table 1 (code for the clustering)
# t20180213 size
# 20171121 AD7 prime
# connect AD6 + AD7 = AD7' (Original idea of Prof Chang for two t)
# ref: AD7_size.r, 20171115_AD6check.r
# t20171101.r
# 20171101 AD7

rm(list=ls(all=TRUE))
library(parallel)
#install.packages('foreach')
library(foreach)
#install.packages('doParallel')
library(doParallel)

#dir_path = "/home/whsiang/NA/sh"
# dir_path="/home/whsiang/NA/sh" #for sourcing
dir_path2=paste("/home/whsiang/NA/theta",sep="") #data output
#/home/whsiang/NA/rdata 
#C:\\Hsiang\\paper\\twoSamples\\bias_computertest_data

parameters = list()
# starting core index for parallel computing
parameters$core_ind = 0
# calculate the number of cores for parallel computing
no_cores <- 1
# nsplit: number of parallel tasks for computing the nrep results
parameters$split_set = (parameters$core_ind * no_cores) + (1:no_cores)
split = parameters$split_set #for testing
# number of datasets to be generated per parallel task
parameters$nrep_sub = 1
# number of bootstrap samples
parameters$nboot = 1000 

# first entry is main alpha level of interest; second entry is the secondary alpha level
parameters$alpha_vec = c(0.05, 0.01)
# starting seed for each replication of data generation
parameters$seedstart = 0

parameters$n=c(3, 4)
#nrep_sub<- 10 #for testing
parameters$aa<- c(3, 4)
parameters$bb<- c(4, 7)
parameters$r<- c(0.5,1) #bias

#####function#####
theta_func = function(s, t, TT, T12, WWeights, weights, ns, ns_sum, nn, H_hat){
  #testing
  # s = "t"
  # t = "t"
  # TT = list(T1,T2)
  # T12=T12
  # WWeights = c(WWeights1, WWeights2)
  # weights = list(weights1,weights2)
  # ns = c(n1,n2)
  # ns_sum = n
  # nn = nn
  # H_hat = H_hat
  #####
  
  kk = c(ns[1]/ns_sum, ns[2]/ns_sum)
  if(s==t){
    theta_hat = 0
    for(j in 1:2){ #(t,t) #test j=1
      Indicatorf = matrix(rep(TT[[j]],nn-1), nn-1, ns[j], byrow=TRUE) <= matrix(rep(T12[-nn],ns[j]), nn-1,ns[j], byrow=FALSE)
      denominator = matrix(rep((ns_sum * kk[j]^2 * weights[[j]]^2),nrow(Indicatorf)),nrow(Indicatorf),ns[j],byrow = TRUE)
      theta_hat_temp = ((WWeights[j])^2) * ((Indicatorf-H_hat)^2)/denominator
      theta_hat_temp = apply(theta_hat_temp, 1, sum)
      theta_hat = theta_hat + theta_hat_temp
    }
  }else{ #(s,t)or(t,s)
    theta_hat = 0
    matrix_temp = list(matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1),matrix(rep(NA,(nn-1)*(nn-1)),nn-1,nn-1))

    for(j in 1:2){ #j=1
      Indicatorf = matrix(rep(TT[[j]],nn-1), nn-1, ns[j], byrow=TRUE) <= matrix(rep(T12[-nn],ns[j]), nn-1, ns[j], byrow=FALSE)
      s_matrix = Indicatorf - H_hat
      t_matrix = Indicatorf - H_hat
      WW = WWeights[j]^2
      denominator = ns_sum * kk[j]^2 * weights[[j]]^2
      # s_matrix = s_matrix[-nrow(s_matrix),]
      # t_matrix = t_matrix[-1,]
      # denominator = matrix(rep((ns_sum * kk[j]^2 * weights[[j]]^2),nrow(s_matrix)),nrow(s_matrix),ns[j],byrow = TRUE)
      # matrix_temp = WW * s_matrix * t_matrix / denominator
      # theta_hat_temp = apply(matrix_temp, 1, sum)
      for(t in 2:(nn-1)){ #t=2
        for(s in 1:(t-1)){ #s=1
          sum_n = 0
          theta_hat_temp = c()
          nj_vec = c()
          for(ns_temp in 1:ns[j]){ #ns_temp=1
            nj_temp = WW * s_matrix[s,ns_temp] * t_matrix[t,ns_temp] / denominator[ns_temp]
            nj_vec = c(nj_vec,nj_temp)
          }
          matrix_temp[[j]][s,t] = sum(nj_vec)
        }
      }
      theta_hat = matrix_temp[[1]] + matrix_temp[[2]]
    }
  }
  return(theta_hat)
}
##################

powerfn = function(parameters) {
  nboot = parameters$nboot
  nrep_sub = parameters$nrep_sub
  #alpha_vec = parameters$alpha_vec
  seedstart = parameters$seedstart
  aa = parameters$aa
  bb = parameters$bb
  
  n = sum(parameters$n)
  n1 = parameters$n[1]
  n2 = parameters$n[2]
  r = parameters$r
  
  # data cdf and pdfs
  # Under H0: Beta(1,2)
  # shape1s = list(aa, bb) # for F0
  # F0<- function(t){pbeta(t,shape1=shape1s[[1]][1],shape2=shape1s[[2]][2])}
  # f0<- function(t){dbeta(t,shape1=shape1s[[1]][1],shape2=shape1s[[2]][2])}
  # Under H1
  # a_beta<- 5
  # F1<- function(t){t^a_beta}
  # f1<- function(t){a_beta*t^(a_beta-1)}
  
#nrep_sub = 1#for test
  
  rej.rate6_t <- rep(NA, nrep_sub)
  rej.rate7_t <- rep(NA, nrep_sub)
  rej.rate7p_t <- rep(NA, nrep_sub)
  rej.rate6_o <- rep(NA, nrep_sub)
  rej.rate7_o <- rep(NA, nrep_sub)
  rej.rate7p_o <- rep(NA, nrep_sub)
  
#d=1#for test
  for(d in 1:nrep_sub){
    
    set.seed(0+nrep_sub*(split-1)+d)
    # tic<-Sys.time()
    T1 <- rbeta(n1, aa[1]+r[1], bb[1]) #X1
    T2 <- rbeta(n2, aa[2]+r[2], bb[2]) #X2
    T12 <- sort(unique(c(T1,T2)))
    nn = length(T12)
    
    weights1=T1 ^ r[1] # #w1
    weights2=T2 ^ r[2] # #w2
    weights = c(weights1,weights2)
    WWeights1 = n1/(sum(1/weights1)) #W1
    WWeights2 = n2/(sum(1/weights2)) #W2
    WWeights  = c(WWeights1,WWeights2)
    
    kk = c(n1/n,n2/n)
    
    p_hat1 = WWeights1/(n1*weights1)
    p_hat2 = WWeights2/(n2*weights2)
    
    #for loop check
    # ECDF1test = data.frame()
    # TF = data.frame()
    # for (tt in 1:(n-1)){ #t
    #   for (i in 1:n1){
    #     ECDF1test[tt,i] = (T1[i] <= T12[tt]) * p_hat1[i]
    #     TF[tt,i] = (T1[i] <= T12[tt])
    # }}
    # ECDF1testsum = apply(ECDF1test, 1, sum)
    
    # ECDF1 <- ecdf(T1)
    # ECDF2 <- ecdf(T2)
    
    ECDF1 <- apply((matrix(rep(p_hat1,nn-1), nn-1, n1, byrow=TRUE) * (matrix(rep(T1,nn-1), nn-1, n1, byrow=TRUE) <= matrix(rep(T12[-nn],n1), nn-1, n1, byrow=FALSE))), 1, sum) #F_hat1
    ECDF2 <- apply((matrix(rep(p_hat2,nn-1), nn-1, n2, byrow=TRUE) * (matrix(rep(T2,nn-1), nn-1, n2, byrow=TRUE) <= matrix(rep(T12[-nn],n2), nn-1, n2, byrow=FALSE))), 1, sum) #F_hat2
    
    ##0811 ok    
    ##Error in { : 
    ##    task 1 failed - "Lapack routine dgesv: system is exactly singular: U[2,2] = 0"
    ##https://statisticsglobe.com/r-error-in-solve-system-is-exactly-singular
    
    # H_hat=-Inf? Yes  0716
    # bias H_hat=Inf? Yes 0812
    # H_hat =  ( n1*ECDF1(T12[-n]) + n2*ECDF2(T12[-n]) )/n
    # H_hat =  ( n1*ECDF1(T12) + n2*ECDF2(T12) )/n 
    H_hat =  ( n1*ECDF1 + n2*ECDF2 )/n
    
    #dH_hatvec = c(0, H_hat[-n])
    dH_hatvec = c(0, H_hat)
    
    theta_hat_t_t = theta_func("t","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)

    #-Inf without[-n]
    AD6_o<- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat*(1-H_hat)) * diff(dH_hatvec))
    AD6_t<- n * sum( ((ECDF1 - ECDF2)^2) / (theta_hat_t_t) * diff(dH_hatvec))
    
    #AD7
    theta_hat_s_t = theta_func("s","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)
    
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
    for(k in 2:(nn-1)){  #t k=2
      for(j in 1:(k-1)){ #s j=1
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
    rej.rate6_t[d]<- AD6_t  > quantile(AD6_star_t,0.95)
    rej.rate7_t[d]<- AD7_t  > quantile(AD7_star_t,0.95)
    rej.rate7p_t[d]<- AD7p_t > quantile(AD7p_star_t,0.95)
    
    AD7p_star_o  <- AD6_star_o + AD7_star_o
    rej.rate6_o[d]<- AD6_o  > quantile(AD6_star_o,0.95)
    rej.rate7_o[d]<- AD7_o  > quantile(AD7_star_o,0.95)
    rej.rate7p_o[d]<- AD7p_o > quantile(AD7p_star_o,0.95)
    
    # toc<-Sys.time()
    # as.numeric(toc-tic,units="mins") # 0.82
    print(d)
    print(Sys.time())
  }
  
  out_matrix = matrix(c(mean(rej.rate6_o),mean(rej.rate6_t),
                        mean(rej.rate7_o),mean(rej.rate7_t),
                        mean(rej.rate7p_o),mean(rej.rate7p_t)),2,3)
  out_array=list(
    rej.rate6_t=rej.rate6_t,
    rej.rate7_t=rej.rate7_t,
    rej.rate7p_t=rej.rate7p_t,
    rej.rate6_o=rej.rate6_o,
    rej.rate7_o=rej.rate7_o,
    rej.rate7p_o=rej.rate7p_o
  )
  return(list(out_matrix,out_array))
} # EnD powerfn

############# calculating coverage given n_subject,  mu,  sigma (utilizing R parallel computing)
##.libPaths("C:\\Users\\Administrator\\Dropbox\\computing\\R_package\\Rlib")
# initiate cluster
cl = makeCluster(no_cores)
#cl = makeCluster(no_cores, type = "PSOCK") # type of cluster
registerDoParallel(cl)
##clusterEvalQ(cl, .libPaths("C:\\Users\\Administrator\\Dropbox\\computing\\R_package\\Rlib"))
### compute the results:
time1 = Sys.time()
foreach(split = parameters$split_set,
        .combine = c  ,
        .packages = c("parallel", "foreach", "iterators", "doParallel", "rootSolve", "Rcpp", "plyr", "coda", "numDeriv", "bbmle", "fdANOVA", "tryCatchLog", "futile.logger")
) %dopar% {
  out = powerfn(parameters)
  # above: to-do
  setwd(dir_path2)
  save(out, file=paste("both_n_", paste(parameters$n, collapse = "_"),"_a_",paste(parameters$aa, collapse = "_"),"_b_",paste(parameters$bb, collapse = "_"),"_r_",paste(parameters$r, collapse = "_"),"_split_",split,".Rdata",sep=""))
}  # END foreach
time2 = Sys.time()
time2-time1
stopCluster(cl)


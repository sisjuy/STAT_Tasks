#AD6 code for 2 sample function
theta_func = function(s, t, TT, T12, WWeights, weights, ns, ns_sum, nn, H_hat){
  #testing
  # s = "t"
  # t = "t"
  # TT = list(T1,T2)
  # T12=T12
  # WWeights = c(WWeights1, WWeights2)
  # weights = list(weights1,weights2)
  # ns = c(n1,n2,n3....)
  # ns = nj
  # ns_sum = n
  # nn = nn
  # H_hat = H_hat
  #####
  
  kk = c(ns[1]/ns_sum, ns[2]/ns_sum) #kk=kappa
  theta_hat = 0
  if(s==t){
    
    for(j in 1:2){ #(t,t) #test j=1
      Indicatorf = matrix(rep(TT[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn,ns[j], byrow=FALSE)
      
      denominator = matrix(rep((ns_sum * kk[j]^2 * weights[[j]]^2),nrow(Indicatorf)),nrow(Indicatorf),ns[j],byrow = TRUE)
      
      theta_hat_temp = ((WWeights[j])^2) * ((Indicatorf-H_hat)^2)/denominator
      theta_hat_temp_sum = apply(theta_hat_temp, 1, sum)
      theta_hat = theta_hat + theta_hat_temp_sum
    }
  }else{ #(s,t)or(t,s)
    matrix_temp = list(matrix(rep(NA,(nn)*(nn)),nn,nn),matrix(rep(NA,(nn)*(nn)),nn,nn))
    
    for(j in 1:2){ #j=1
      Indicatorf = matrix(rep(TT[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn, ns[j], byrow=FALSE)
      s_matrix = Indicatorf - H_hat
      t_matrix = Indicatorf - H_hat
      WW = WWeights[j]^2
      denominator = ns_sum * kk[j]^2 * weights[[j]]^2
      
      for(t in 2:(nn)){ #t=2
        for(s in 1:(t)){ #s=1
          sum_n = 0
          theta_hat_temp = c()
          nj_vec = c()
          for(ns_temp in 1:ns[j]){ #ns_temp=3
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


T1 = c(1,2)
T2 = c(2,3)
data = list(T1,T2)
r = c(0.5,1)
nboot = 10
T1 = data[[1]]
T2 = data[[2]]
n1 = length(T1)
n2 = length(T2)
n = n1+n2
lower_bound=max(min(T1),min(T2))
upper_bound=min(rev(unique(sort(T1)))[2],rev(unique(sort(T2)))[2]) 
T12_temp <- sort(unique(c(T1,T2)))
nn_temp = length(T12_temp)
lowerbindx=(1:nn_temp)[abs(T12_temp-lower_bound)==min(abs(T12_temp-lower_bound))]
upperbindx=(1:nn_temp)[abs(T12_temp-upper_bound)==min(abs(T12_temp-upper_bound))]
T12 = T12_temp[lowerbindx:upperbindx]
nn = length(T12)

weights1=T1 ^ r[1] # #w1
weights2=T2 ^ r[2] # #w2
weights = c(weights1,weights2)
WWeights1 = n1/(sum(1/weights1)) #W1
WWeights2 = n2/(sum(1/weights2)) #W2
WWeights  = c(WWeights1,WWeights2)

kappa = c(n1/n,n2/n)

p_hat1 = WWeights1/(n1*weights1)
p_hat2 = WWeights2/(n2*weights2)
ECDF1 <- apply((matrix(rep(p_hat1,nn), nn, n1, byrow=TRUE) * (matrix(rep(T1,nn), nn, n1, byrow=TRUE) <= matrix(rep(T12,n1), nn, n1, byrow=FALSE))), 1, sum) #F_hat1
ECDF2 <- apply((matrix(rep(p_hat2,nn), nn, n2, byrow=TRUE) * (matrix(rep(T2,nn), nn, n2, byrow=TRUE) <= matrix(rep(T12,n2), nn, n2, byrow=FALSE))), 1, sum) #F_hat2
H_hat =  ( n1*ECDF1 + n2*ECDF2 )/n
#k sample weighted sum sum(nj*Fj/n)


#dH_hatvec = c(0, H_hat[-n])
dH_hatvec = c(0, H_hat)

theta_hat_t_t = theta_func("t","t",list(T1,T2),T12,c(WWeights1, WWeights2),list(weights1,weights2),c(n1,n2),n,nn,H_hat)

#-Inf without[-n]
AD6_o<- ((n1*n2)/n) * sum( ((ECDF1 - ECDF2)^2) / (H_hat*(1-H_hat)) * diff(dH_hatvec))
AD6_t<- n * sum( ((ECDF1 - ECDF2)^2) / (theta_hat_t_t) * diff(dH_hatvec))


#bootstrap
AD6_star_t <- rep(NA,nboot)
AD6_star_o <- rep(NA,nboot)
for(b in 1:nboot){
  b=1
  xi1<- rnorm(n1, 0, 1)
  xi2<- rnorm(n2, 0, 1)
  
  D1_star <- apply(matrix(rep((xi1 * p_hat1),nn), nn, n1, byrow=TRUE) * ((matrix(rep(T1,nn), nn, n1, byrow=TRUE) <= matrix(rep(T12,n1), nn, n1, byrow=FALSE)) - H_hat), 1, sum)
  D2_star <- apply(matrix(rep((xi2 * p_hat2),nn), nn, n2, byrow=TRUE) * ((matrix(rep(T2,nn), nn, n2, byrow=TRUE) <= matrix(rep(T12,n2), nn, n2, byrow=FALSE)) - H_hat), 1, sum)
  
  #AD6
  # AD_s_unit_6_I1 <- apply(matrix(rep(xi1,n-1), n-1, n1, byrow=TRUE) * ((matrix(rep(T1,n-1), n-1, n1, byrow=TRUE) <= matrix(rep(T12[-n],n1), n-1, n1, byrow=FALSE)) - H_hat[-n]), 1, mean)
  # AD_s_unit_6_I2 <- apply(matrix(rep(xi2,n-1), n-1, n2, byrow=TRUE) * ((matrix(rep(T2,n-1), n-1, n2, byrow=TRUE) <= matrix(rep(T12[-n],n2), n-1, n2, byrow=FALSE)) - H_hat[-n]), 1, mean)
  
  #new dinominater for bootstrap 20220406
  dinominater_n1_data  <- matrix(rep((xi1 * WWeights1 / (sqrt(kappa[1])*weights1)),nn), nn, n1, byrow=TRUE) * ((matrix(rep(T1,nn), nn, n1, byrow=TRUE) <= matrix(rep(T12,n1), nn, n1, byrow=FALSE)) - H_hat)
  dinominater_n1_data_mean  <- apply(dinominater_n1_data,1,sum)/n1
  #data-mean
  
  dinominater_n1_temp1 <- dinominater_n1_data - dinominater_n1_data_mean
  #var
  dinominater_n1_sum   <- apply(dinominater_n1_temp1^2,1,sum)/n1
  
  dinominater_n2_data  <- matrix(rep((xi2 * WWeights2 / (sqrt(kappa[2])*weights2)),nn), nn, n2, byrow=TRUE) * ((matrix(rep(T2,nn), nn, n2, byrow=TRUE) <= matrix(rep(T12,n2), nn, n2, byrow=FALSE)) - H_hat)
  dinominater_n2_data_mean  <- apply(dinominater_n2_data,1,sum)/n2
  dinominater_n2_temp1 <- dinominater_n2_data - dinominater_n2_data_mean
  dinominater_n2_sum   <- apply(dinominater_n2_temp1^2,1,sum)/n2
  
  dinominater_bootstrap <- dinominater_n1_sum + dinominater_n2_sum
  
  a11_t_dinominater <- dinominater_bootstrap #A2 second matrix[1,1]
  #0413 COV
  a21_t_dinominater <- COV_func(dinominater_n1_data,dinominater_n1_data_mean,n1,nn) + COV_func(dinominater_n2_data,dinominater_n2_data_mean,n2,nn) #A2 second matrix[2,1]&[1,2]
  a22_t_dinominater <- a11_t_dinominater #A2 second matrix[2,2]
  
  #for testing: s==t dinominater_bootstrap==a21_t_dinominater?
  #test=c()
  #for(i in c(1:nn)){
  #  test = c(test,a21_t_dinominater[i,i])
  #}
  #dinominater_bootstrap == test
  
  AD_s_unit_6 <- (D1_star - D2_star) ^ 2
  AD6_star_o[b] <- (n1*n2/n)*sum(AD_s_unit_6/(H_hat*(1-H_hat))*diff(dH_hatvec), na.rm = TRUE)
  AD6_star_t[b] <- n*sum(AD_s_unit_6/(dinominater_bootstrap)*diff(dH_hatvec), na.rm = TRUE)
}
p6_t <- mean(AD6_star_t >= AD6_t,na.rm=TRUE)
p6_o <- mean(AD6_star_o >= AD6_o,na.rm=TRUE)
rej.rate6_t <- AD6_t  > quantile(AD6_star_t,0.95)
out_array=list(
  rej.rate6_o=rej.rate6_o,
  rej.rate6_t=rej.rate6_t
)
p = list(
  p6_o = p6_o,
  p6_t = p6_t
)
return(c(out_array,p))
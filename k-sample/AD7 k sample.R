#big_theta_hat_j
a11=list()
a22=list()
a21=list()

for(i in 1:length(data)){
  a11[[i]] = theta_hat_t_t$theta_hat_j[[i]]
  a21[[i]] = theta_hat_s_t$theta_hat_j[[i]]
}
a22 = a11

big_theta_hat_j[[i]] = matrix(c(theta_hat_t_t$theta_hat_j[[i]][j],theta_hat_s_t$theta_hat_j[[i]][j,k],theta_hat_s_t$theta_hat_j[[i]][j,k],theta_hat_t_t$theta_hat_j[[i]][k]),2,2)

#big_theta_check
s=0
for(i in 1:length(data)){
  s = s + solve(big_theta_hat_j[[i]])
}
big_theta_check_j = 1/s
#problem: big_theta_hat_j can not inverse


#Nj
Nj=list()
Nj[[i]] = (big_theta_check_j^(0.5))*(1/big_theta_hat_j[[i]])*big_theta_check_j^(0.5)

#spsi_hat_j
#D_hat_j[[i]]
spsi_hat_j=list()
spsi_hat_j[[i]] = sqrt(n)*(big_theta_check_j^(-0.5))*t(c(D_hat_j[[i]][j],D_hat_j[[i]][k]))

#spsi_check_j
spsi_check_j=list()
Nj[[i]]



#1:k
for(i in 1:length(data)){
  a11[[i]] = theta_hat_t_t$theta_hat_j[[i]]
  a21[[i]] = theta_hat_s_t$theta_hat_j[[i]]
}
a22 = a11



for(k in 2:(nn)){  #t k=2
  for(j in 1:(k-1)){ #s j=1
    Vinv_t <- solve(matrix(c(a11[[i]][j], a21[[i]][j,k], a21[[i]][j,k], a22[[i]][k]),2,2))
    
    dH_hat_t_vec <- c(0, H_hat[k])
    dH_hat_s_vec <- c(0, H_hat[j])
    
    AD7_unit_t[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv_t%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    
  }
}



########################
library(expm)

AD7_unit <- matrix(rep(NA,(nn)*(nn)),nn,nn)
for(k in 2:(nn)){  #t k=3
  for(j in 1:(k-1)){ #s j=2
    
    big_theta_hat_j = list()
    #big_theta_check_j = list()
    cal = 0 
    for(i in 1:length(data)){
      #i=1
      a11 = theta_hat_t_t$theta_hat_j[[i]]
      a21 = theta_hat_s_t$theta_hat_j[[i]]
      a22 = a11
      big_theta_hat_j[[i]] = matrix(c(a11[j], a21[j,k], a21[j,k], a22[k]),2,2)
      cal = cal + solve(big_theta_hat_j[[i]])
      #big_theta_check_j[[i]] = solve(big_theta_hat_j[[i]]) #A2 second matrix inv
    }
    big_theta_check_j = solve(cal)
    
    
    Nj = list()
    spsi_hat_j = list()
    for(i in 1:length(data)){
      Nj[[i]] = sqrtm(big_theta_check_j)%*%solve(big_theta_hat_j[[i]])%*%sqrtm(big_theta_check_j)
      spsi_hat_j[[i]] = sqrt(n)*solve(sqrtm(big_theta_check_j))%*%c(D_hat_j[[i]][j],D_hat_j[[i]][k])
    }
    
    
    #spsi_check
    spsi_check = 0
    for(i in 1:length(data)){
      i=2
      spsi_check = spsi_check + Nj[[i]]%*%spsi_hat_j[[i]]
    }
    
    dH_hat_t_vec <- c(0, H_hat[k])
    dH_hat_s_vec <- c(0, H_hat[j])
    
    #BSSB
    BSSB = 0
    for(i in 1:length(data)){
      BSSB = BSSB + t(spsi_hat_j[[i]]-spsi_check)%*%Nj[[i]]%*%(spsi_hat_j[[i]]-spsi_check) 
    }
    
    AD7_unit[j,k] = BSSB * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    
  }
}
AD7_ksample = sum(AD7_unit, na.rm=TRUE)

#####AD7 bootstrap#####
big_theta_hat_j_star = list()
AD7_unit_star <- matrix(rep(NA,(nn)*(nn)),nn,nn)
for(k in 2:(nn)){  #t k=2
  for(j in 1:(k-1)){ #s j=1
    
    cal = 0
    for(i in 1:length(data)){
      #i=1
      a11_star = theta_hat_j_star_tt[[i]]
      a21_star = theta_hat_j_star_st[[i]]
      a22_star = a11_star
      big_theta_hat_j_star[[i]] = matrix(c(a11_star[j], a21_star[j,k], a21_star[j,k], a22_star[k]),2,2)
      cal = cal + solve(big_theta_hat_j_star[[i]]) #A2 second matrix inv
    }
    big_theta_check_j_star = solve(cal)

    Nj_star = list()
    spsi_hat_j_star = list()
    for(i in 1:length(data)){
      Nj_star[[i]] = sqrtm(big_theta_check_j_star)*(solve(big_theta_hat_j_star[[i]]))*sqrtm(big_theta_check_j_star)
      spsi_hat_j_star[[i]] = sqrt(n)*(solve(sqrtm(big_theta_check_j_star))%*%c(Dj_star[[i]][j],Dj_star[[i]][k]))
    }
    
    
    #spsi_check
    spsi_check_star = 0
    for(i in 1:length(data)){
      spsi_check_star = spsi_check_star + Nj_star[[i]]%*%spsi_hat_j_star[[i]]
    }
    
    
    dH_hat_t_vec <- c(0, H_hat[k])
    dH_hat_s_vec <- c(0, H_hat[j])
    
    #BSSB
    BSSB_star = 0
    for(i in 1:length(data)){
      BSSB_star = BSSB_star + t(spsi_hat_j_star[[i]]-spsi_check_star)%*%Nj_star[[i]]%*%(spsi_hat_j_star[[i]]-spsi_check_star)
    }
    
    AD7_unit_star[j,k] = BSSB_star * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
  }
}

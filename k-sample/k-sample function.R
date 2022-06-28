


#test data
a = c(1,2,3,4,2,3,6,7,3)
c = c(1,1,1,1,2,2,2,2,2,3,4,7,10)
d = c(1,3,4,5,6,11,12)
data = list(a,c,d)
r=c(0.5,0.5,1)
nboot = 10


#theta function:
theta_func = function(s, t, TT, T12, WWeights, weights, ns, ns_sum, nn, H_hat){
  #testing
  # s = "s"
  # t = "t"
  # TT = data (T1,T2,T3)
  # T12=T12
  # WWeights = WWeights
  # weights = weights
  # ns = ns = c(n1,n2,n3....)
  # ns = nj
  # ns_sum = n
  # nn = nn
  # H_hat = H_hat
  #####
  
  # TT = data
  
  
  #kk: kappa
  kk = c()
  for(i in 1:length(data)){
    kk = c(kk,ns[i]/ns_sum)
  }
  
  theta_hat = 0
  theta_hat_j = list()
  if(s==t){
    
    for(j in 1:length(data)){ #(t,t) #test j=1
      
      Indicatorf = matrix(rep(TT[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn,ns[j], byrow=FALSE)
      
      denominator = matrix(rep((ns_sum * kk[j]^2 * weights[[j]]^2),nrow(Indicatorf)),nrow(Indicatorf),ns[j],byrow = TRUE)
      
      theta_hat_temp = ((WWeights[j])^2) * ((Indicatorf-H_hat)^2)/denominator
      theta_hat_temp_sum = apply(theta_hat_temp, 1, sum)
      theta_hat_j[[j]] = theta_hat_temp_sum
      theta_hat = theta_hat + theta_hat_temp_sum
    }
    return(list(theta_hat=theta_hat,theta_hat_j=theta_hat_j))
  }else{ #(s,t)or(t,s)
    matrix_temp = list()
    
    for(i in 1:length(data)){
      matrix_temp[[i]] = matrix(rep(NA,(nn)*(nn)),nn,nn)
    }

    
    for(j in 1:length(data)){ #j=1
      Indicatorf = matrix(rep(TT[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn, ns[j], byrow=FALSE)
      s_matrix = Indicatorf - H_hat
      t_matrix = Indicatorf - H_hat
      WW = WWeights[j]^2
      denominator = ns_sum * kk[j]^2 * weights[[j]]^2
      
      for(t in 2:(nn)){ #t=2
        for(s in 1:(t-1)){ #s=1
          sum_n = 0
          theta_hat_temp = c()
          nj_vec = c()
          for(ns_temp in 1:ns[j]){ #ns_temp=3
            nj_temp = WW * s_matrix[s,ns_temp] * t_matrix[t,ns_temp] / denominator[ns_temp]
            nj_vec = c(nj_vec,nj_temp)
          }
          matrix_temp[[j]][s,t] = sum(nj_vec)
          #print(matrix_temp)
        }
        
      }
      #print(matrix_temp)
      
      theta_hat = Reduce("+", matrix_temp)
      
    }
    return(list(theta_hat=theta_hat,theta_hat_j=matrix_temp))
  }
  
}

#cov function:
COV_func = function(data, mean, ns, nn){
  ###test
  # data = dinominater_n1_data
  # mean = dinominater_n1_data_mean
  # ns = n1
  # nn = nn
  ###
  matrix_temp = matrix(rep(NA,(nn)*(nn)),nn,nn)
  
  s_matrix = data - mean
  t_matrix = data - mean
  
  for(t in 2:(nn)){ #t=2
    for(s in 1:(t)){ #s=1
      nj_vec = c()
      for(ns_temp in 1:ns){ #ns_temp=3
        nj_temp = s_matrix[s,ns_temp] * t_matrix[t,ns_temp]
        nj_vec = c(nj_vec,nj_temp)
      }
      matrix_temp[s,t] = sum(nj_vec)
    }
  }
  return(matrix_temp/ns)  
}


###input format: 
#a = c(1,2,3,4,2,3,6,7,3)
#c = c(1,1,1,1,2,2,2,2,2,3,4,7,10)
#d = c(1,3,4,5,6,11,12)
#r = c(0.5,0.5,1)
#data = list(a,c,d)
#nboot = 10



####main function
BS_ksample <- function(data, r, nboot){
  # n = n1 + n2 + ...
  # ns = c(n1,n2,n3)
  # Tj = data[[j]]
  # alldata = c(T1,T2,...)
  alldata = c()
  mindata = c()
  ns = c()
  cal = c()
  weights = list()
  WWeights = c()
  
  for(i in 1:length(data)){
    #i=1
    alldata = c(alldata,data[[i]])
    mindata = c(mindata,min(data[[i]]))
    ns = c(ns,length(data[[i]]))
    cal = c(cal,rev(unique(sort(data[[i]])))[2])
    weights[[i]] = data[[i]]^r[i]
    WWeights = c(WWeights,ns[i]/(sum(1/weights[[i]])))
  }
  n = sum(ns)
  lower_bound = max(mindata)
  upper_bound = min(cal)
  T12_temp <- sort(unique(alldata))
  nn_temp = length(T12_temp)
  lowerbindx=(1:nn_temp)[abs(T12_temp-lower_bound)==min(abs(T12_temp-lower_bound))]
  upperbindx=(1:nn_temp)[abs(T12_temp-upper_bound)==min(abs(T12_temp-upper_bound))]
  T12 = T12_temp[lowerbindx:upperbindx]
  nn = length(T12)
  
  kappa = c()
  p_hat = list()
  ECDF = list()
  H_hat = 0
  for(i in 1:length(data)){
    kappa = c(kappa,ns[i]/n)
    p_hat[[i]] = WWeights[i]/(ns[i]*weights[[i]])
    ECDF[[i]] = apply((matrix(rep(p_hat[[i]],nn), nn, ns[i], byrow=TRUE) * (matrix(rep(data[[i]],nn), nn, ns[i], byrow=TRUE) <= matrix(rep(T12,ns[i]), nn, ns[i], byrow=FALSE))), 1, sum)
    H_hat = H_hat + ns[i]*ECDF[[i]]/n
  }
  dH_hatvec = c(0, H_hat)
  theta_hat_t_t = theta_func("t","t",data,T12,WWeights,weights,ns,n,nn,H_hat)
  theta_hat_s_t = theta_func("s","t",data,T12,WWeights,weights,ns,n,nn,H_hat)
  
  #####AD6#####
  
  #nu_hat_j 
  nu_hat_j_sum = 0
  for(j in 1:length(data)){
    nu_hat_j_sum = nu_hat_j_sum + 1/theta_hat_t_t$theta_hat_j[[j]]
  }
  nu_hat_j = list()
  for(j in 1:length(data)){
    #print(1/theta_hat_t_t$theta_hat_j[[j]]/nu_hat_j_sum)
    nu_hat_j[[j]] = 1/theta_hat_t_t$theta_hat_j[[j]]/nu_hat_j_sum
  }
  
  #psi_hat_j
  #psi_hat_j = sqrt(n)*(ECDF[[j]]-H_hat)/sqrt(theta_hat_t_t$theta_hat_j[[j]]*nu_hat_j[[j]])
  
  psi_check = 0 
  for(j in 1:length(data)){
    psi_hat_j = sqrt(n)*(ECDF[[j]]-H_hat)/sqrt(theta_hat_t_t$theta_hat_j[[j]]*nu_hat_j[[j]])
    psi_check = psi_check + nu_hat_j[[j]]*psi_hat_j
  }
  
  #D_hat_j for original version
  D_hat_j = list()
  for(j in 1:length(data)){
    D_hat_j[[j]] = ECDF[[j]] - H_hat
  }
  
  #SSO: original version
  #SSB
  SSO = 0
  SSB = 0 
  for(i in 1:length(data)){
    SSO = SSO + ns[[j]]*(D_hat_j[[j]])^2/(H_hat*(1-H_hat))
    psi_hat_j = sqrt(n)*(ECDF[[j]]-H_hat)/sqrt(theta_hat_t_t$theta_hat_j[[j]]*nu_hat_j[[j]])
    SSB = SSB + nu_hat_j[[j]]*(psi_hat_j - psi_check)^2
  }
  
  AD6_o = sum(SSO*diff(dH_hatvec))
  AD6_ksample = sum(SSB*diff(dH_hatvec))
  
  #####AD7#####
  library(expm)
  a11_o <- H_hat*(1-H_hat) 
  a21_o <- H_hat%*%t(1-H_hat)
  a22_o <- a11_o
  AD7_unit_o <- matrix(rep(NA,(nn)*(nn)),nn,nn)
  AD7_unit_k <- matrix(rep(NA,(nn)*(nn)),nn,nn)
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
      
      
      #original version
      Vinv_o <- solve(matrix(c(a11_o[j], a21_o[j,k], a21_o[j,k], a22_o[k]),2,2))
      
      #MU
      MU = 0
      for(i in 1:length(data)){
        #i=1
        #j=1
        #k=2
        MU = MU + ns[i]*t(c(D_hat_j[[i]][j],D_hat_j[[i]][k]))%*%Vinv_o%*%c(D_hat_j[[i]][j],D_hat_j[[i]][k])
      }
      
      
      AD7_unit_k[j,k] = BSSB * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
      AD7_unit_o[j,k] = MU * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    }
  }
  AD7_ksample = sum(AD7_unit_k, na.rm=TRUE)
  AD7_o = sum(AD7_unit_o, na.rm = TRUE)
  
  #######bootstrap#######
  AD6_star_k <- rep(NA,nboot)
  AD6_star_o <- rep(NA,nboot)
  AD7_star_k <- rep(NA,nboot)
  AD7_star_o <- rep(NA,nboot)
  #nboot=200
  for(b in 1:nboot){
    #b=1
    xi_list = list()
    for(j in 1:length(data)){
      xi_list[[j]] = rnorm(ns[j],0,1)
    }
    
    
    #j=1
    
    #theta_hat_j(t,t)
    Dj_star = list()
    dinominater_bootstrap = 0
    nu_hat_star_collect = list()
    nu_hat_star_sum = 0
    theta_hat_j_star_tt = list()
    theta_hat_j_star_st = list()
    for(j in 1:length(data)){
      
      Dj_star[[j]] <- apply(matrix(rep((xi_list[[j]] * p_hat[[j]]),nn), nn, ns[j], byrow=TRUE) * ((matrix(rep(data[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn, ns[j], byrow=FALSE)) - H_hat), 1, sum)
      #V*(t)
      dinominater_data_j  <- matrix(rep((xi_list[[j]] * WWeights[j] / (sqrt(kappa[j])*weights[[j]])),nn), nn, ns[j], byrow=TRUE) * ((matrix(rep(data[[j]],nn), nn, ns[j], byrow=TRUE) <= matrix(rep(T12,ns[j]), nn, ns[j], byrow=FALSE)) - H_hat)
      #V*(t) mean
      dinominater_data_j_mean  <- apply(dinominater_data_j,1,sum)/ns[j]
      #(data-mean)^2
      dinominater_temp_j <- dinominater_data_j - dinominater_data_j_mean
      
      dinominater_sum_j <- apply(dinominater_temp_j^2,1,sum)/ns[j]
      theta_hat_j_star_tt[[j]] = dinominater_sum_j
      
      #dinominater_bootstrap = dinominater_bootstrap + dinominater_sum_j
      
      theta_hat_j_star_st[[j]] = COV_func(dinominater_data_j,dinominater_data_j_mean,ns[j],nn)
      nu_hat_star_collect[[j]] = 1/theta_hat_j_star_tt[[j]]
      nu_hat_star_sum = nu_hat_star_sum + 1/theta_hat_j_star_tt[[j]]
      
    }
    #nu_hat_j_star -> sum = 1
    nu_hat_j_star = list()
    
    #x = 0 -> check if sum up equal to 1
    for(j in 1:length(data)){
      nu_hat_j_star[[j]] = nu_hat_star_collect[[j]]/nu_hat_star_sum
      #x = x + nu_hat_j_star[[j]] #ok
    }
    
    
    #psi_check_star
    psi_check_star = 0
    psi_hat_j_star = list()
    for(j in 1:length(data)){
      psi_hat_j_star[[j]] = sqrt(n)*Dj_star[[j]]/sqrt(theta_hat_j_star_tt[[j]]*nu_hat_j_star[[j]])
      psi_check_star = psi_check_star + nu_hat_j_star[[j]]*psi_hat_j_star[[j]]
    }
    
    SSB_star = 0
    SSO_star = 0
    for(j in 1:length(data)){
      D = psi_hat_j_star[[j]] - psi_check
      cal = nu_hat_j_star[[j]]*D^2
      SSB_star = SSB_star + cal
      
      #SSO_star
      SSO_star = SSO_star + ns[[j]]*Dj_star[[j]]/(H_hat*(1-H_hat))
    }
    AD6_star_k[b] = sum(SSB_star*diff(dH_hatvec),na.rm = TRUE)
    AD6_star_o[b] = sum(SSO_star*diff(dH_hatvec),na.rm = TRUE)
    
    
    #####bootstrap AD7#####
    big_theta_hat_j_star = list()
    AD7_unit_star_k <- matrix(rep(NA,(nn)*(nn)),nn,nn)
    AD7_unit_star_o <- matrix(rep(NA,(nn)*(nn)),nn,nn)
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
          Nj_star[[i]] = sqrtm(big_theta_check_j_star)%*%(solve(big_theta_hat_j_star[[i]]))%*%sqrtm(big_theta_check_j_star)
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
        
        Vinv_o_star <- solve(matrix(c(a11_o[j], a21_o[j,k], a21_o[j,k], a22_o[k]),2,2))
        
        #MU
        MU = 0
        for(i in 1:length(data)){
          #i=1
          #j=1
          #k=2
          MU = MU + ns[i]*t(c(Dj_star[[i]][j],Dj_star[[i]][k]))%*%Vinv_o%*%c(Dj_star[[i]][j],Dj_star[[i]][k])
        }
        
        
        AD7_unit_star_o[j,k] = MU * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
        AD7_unit_star_k[j,k] = BSSB_star * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
      }
    }
    AD7_star_o[b] = sum(AD7_unit_star_o, na.rm=TRUE)
    AD7_star_k[b] = sum(AD7_unit_star_k, na.rm=TRUE)
  }
  p6_k <- mean(AD6_star_k >= AD6_ksample,na.rm=TRUE)
  p6_o <- mean(AD6_star_o >= AD6_o,na.rm=TRUE)
  p7_k <- mean(AD7_star_k >= AD7_ksample,na.rm=TRUE)
  p7_o <- mean(AD7_star_o >= AD7_o,na.rm=TRUE)
  rej.rate6_k <- AD6_ksample  > quantile(AD6_star_k,0.95)
  rej.rate6_o <- AD6_o  > quantile(AD6_star_o,0.95)
  rej.rate7_k <- AD7_ksample  > quantile(AD7_star_k,0.95)
  rej.rate7_o <- AD7_o  > quantile(AD7_star_o,0.95)
  
  out_array=list(
    rej.rate6_o=rej.rate6_o,
    rej.rate7_o=rej.rate7_o,
    rej.rate6_k=rej.rate6_k,
    rej.rate7_k=rej.rate7_k
  )
  p = list(
    p6_o = p6_o,
    p7_o = p7_o,
    p6_k = p6_k,
    p7_k = p7_k
    
  )
  return(c(out_array,p))
}


BS_ksample(data, r, nboot)


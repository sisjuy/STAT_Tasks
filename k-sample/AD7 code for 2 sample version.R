b1  <- ECDF1 - ECDF2 #A2 first matrix[1,1]
b2  <- b1 #A2 first matrix[1,2]
a11_t <- theta_hat_t_t #A2 second matrix[1,1]
a11_o <- H_hat*(1-H_hat) 
a21_t <- theta_hat_s_t #A2 second matrix[2,1]&[1,2]
a21_o <- H_hat%*%t(1-H_hat) 
a22_t <- a11_t #A2 second matrix[2,2]
a22_o <- a11_o #A2 second matrix[2,2]

AD7_unit_t <- matrix(rep(NA,(nn)*(nn)),nn,nn)

for(k in 2:(nn)){  #t k=2
  for(j in 1:(k-1)){ #s j=1
    Vinv_t <- solve(matrix(c(a11_t[j], a21_t[j,k], a21_t[j,k], a22_t[k]),2,2)) #A2 second matrix inv
    
    dH_hat_t_vec <- c(0, H_hat[k])
    dH_hat_s_vec <- c(0, H_hat[j])
    
    AD7_unit_t[j,k] <- (t(c(b1[j], b2[k]))%*%Vinv_t%*%c(b1[j], b2[k])) * diff(dH_hat_s_vec) * diff(dH_hat_t_vec)
    
  }
}


AD7_o <- (n1*n2/n)*sum(AD7_unit_o, na.rm=TRUE)  # NA bc for j=1, AD7_unit[1, 1] not defined
AD7_t <- n*sum(AD7_unit_t, na.rm=TRUE)  

AD7p_t <- AD6_t + AD7_t
AD7p_o <- AD6_o + AD7_o
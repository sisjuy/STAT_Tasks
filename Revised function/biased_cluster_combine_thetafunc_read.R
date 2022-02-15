##T/F reading problem: solve 0712

rm(list=ls(all=TRUE))

#rdata output folder
dir_path2=paste("C:\\for_work\\paper\\twoSamples\\theta",sep="")
#C:\\for_work\\paper\\twoSamples\\0722finalver_code

parameters = list()
# starting core index for parallel computing
parameters$core_ind = 0
# calculate the number of cores for parallel computing
no_cores <- 2
# nsplit: number of parallel tasks for computing the nrep results
parameters$split_set = (parameters$core_ind * no_cores) + (1:no_cores)
Set = parameters$split_set
# split = parameters$split_set #for testing
# number of datasets to be generated per parallel task
parameters$nrep_sub = 10
nrep_sub = parameters$nrep_sub
# number of bootstrap samples
parameters$nboot = 1000 #

# first entry is main alpha level of interest; second entry is the secondary alpha level
parameters$alpha_vec = c(0.05, 0.01)
# starting seed for each replication of data generation
parameters$seedstart = 0

parameters$n=c(2, 3)
parameters$aa<- c(1, 1)
parameters$bb<- c(2, 2)
parameters$r<- c(0.5,1) #bias

setwd(dir_path2)

mean_rej_rate6_o = rep(NA, no_cores)
mean_rej_rate7_o = rep(NA, no_cores)
mean_rej_rate7p_o = rep(NA, no_cores)
rej_rate6_o = rep(NA, no_cores*nrep_sub)
rej_rate7_o = rep(NA, no_cores*nrep_sub)
rej_rate7p_o = rep(NA, no_cores*nrep_sub)
mean_rej_rate6_t = rep(NA, no_cores)
mean_rej_rate7_t = rep(NA, no_cores)
mean_rej_rate7p_t = rep(NA, no_cores)
rej_rate6_t = rep(NA, no_cores*nrep_sub)
rej_rate7_t = rep(NA, no_cores*nrep_sub)
rej_rate7p_t = rep(NA, no_cores*nrep_sub)

for (split in Set) {

  load(paste("both_n_", paste(parameters$n, collapse = "_"),"_a_",paste(parameters$aa, collapse = "_"),"_b_",paste(parameters$bb, collapse = "_"),"_r_",paste(parameters$r, collapse = "_"),"_split_",split,".Rdata",sep=""))
  mean_rej_rate6_o[split] = out[[1]][1,1]
  mean_rej_rate7_o[split] = out[[1]][1,2]
  mean_rej_rate7p_o[split]= out[[1]][1,3]
  
  mean_rej_rate6_t[split] = out[[1]][2,1]
  mean_rej_rate7_t[split] = out[[1]][2,2]
  mean_rej_rate7p_t[split]= out[[1]][2,3]
  
  rej_rate6_o[((split-1)*nrep_sub+1):(split*nrep_sub)] = out[[2]][["rej.rate6_o"]]
  rej_rate7_o[((split-1)*nrep_sub+1):(split*nrep_sub)] = out[[2]][["rej.rate7_o"]]
  rej_rate7p_o[((split-1)*nrep_sub+1):(split*nrep_sub)]= out[[2]][["rej.rate7p_o"]]
  
  rej_rate6_t[((split-1)*nrep_sub+1):(split*nrep_sub)] = out[[2]][["rej.rate6_t"]]
  rej_rate7_t[((split-1)*nrep_sub+1):(split*nrep_sub)] = out[[2]][["rej.rate7_t"]]
  rej_rate7p_t[((split-1)*nrep_sub+1):(split*nrep_sub)]= out[[2]][["rej.rate7p_t"]]
  
}

#rej_rate6
#rej_rate7
#rej_rate7p

# mean(mean_rej_rate6)
# mean(mean_rej_rate7)
# mean(mean_rej_rate7p)
# mean_rej_rate6 == mean_rej_rate7

#AD7_data
#AD7_star_data

##hist 0722
# AD7_hist <- hist(AD7_data, 
#      prob=TRUE, 
#      breaks=10,
# )
# AD7_star_hist <- hist(AD7_star_data, 
#      prob=TRUE, 
#      breaks=10,
# )
# 
# plot(AD7_hist, col = "lightblue")
# plot(AD7_star_hist, col = "gray")#, add=TRUE) 
output_rate = matrix(c(mean(rej_rate6_o),mean(rej_rate6_t),
                       mean(rej_rate7_o),mean(rej_rate7_t),
                       mean(rej_rate7p_o),mean(rej_rate7p_t)),2,3)
colnames(output_rate) <- c("AD6", "AD7","AD7p")
rownames(output_rate) <- c("Original Denominator", "New Denominator")
output_rate
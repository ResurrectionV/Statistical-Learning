#Clear the buff
rm(list = ls())

#Import
#File
NG <- read.csv("NG.csv")
#show first 6 rows in NG
head(NG)
#Summarize the csv file
summary(NG)
#linear model to original data
lmod <- lm(delta.temp ~ X.disaster, data = NG)
lmod
#Summarize the linear model
summary(lmod)

#1.
#Build a confidence interval estimator for the mean based on the bootstrap (use 10000 =nboot)
ConfItv.boot <- function(b_vec,nboot=10000,alpha=0.1){
  #Find the number of data, mean of the data, and the standard deviation of the data
  r_num <- length(b_vec)
  mean_0 <- mean(b_vec)
  s.d_0 <- sqrt(var(b_vec))
  #Initiate a NULL variable that will contains a list of cutoff based on sampled b_vec
  boot_vec <- NULL
  #run nboot times with sampling distribution, and each time stores the cutoff in boot_vec
  for (i in 1:nboot){
    #generate new sample of b_vec with replacement
    new_vec <- sample(b_vec,replace=T)
    #Check is var(new_vec) is 0, because it will cause Inf value in the future calculation in cutoff
    #if it is 0, resample the b_vec
    while (var(new_vec) == 0) new_vec <- sample(b_vec)
    #Find the mean and standard deviation of new sample vec
    mean_new <- mean(new_vec)
    s.d_new <- sqrt(var(new_vec))
    #Add new cutoff into the list
    boot_vec <- c(boot_vec, (mean_new - mean_0)/(s.d_new/sqrt(r_num)))
  }
  #Create the lower bound and upper bound of confidence interval based on bot_vec and CL = alpha
  lq <- quantile(boot_vec, alpha/2)
  uq <- quantile(boot_vec, 1 - alpha/2)
  lb <- mean_0 - uq*(s.d_0/sqrt(r_num))
  ub <- mean_0 - lq*(s.d_0/sqrt(r_num))
  #Create the lower bound and upper bound based on t distribution
  nlb <- mean_0 - qt(1-alpha/2,r_num-1)*(mean_0/sqrt(r_num))
  nub <- mean_0 + qt(1-alpha/2,r_num-1)*(mean_0/sqrt(r_num))
  #Output the confidence interval from bootstrap and t distribution respectively
  list(my.boot_interval = c(lb,ub), normal_interval = c(nlb,nub))
}

#2. 
#Build a simulator that draws n samples form a lognormal distribution (rlnorm) and builds both the central limit theorem based confidence interval, and compares it to the coverage rate for the bootstrap (confidence interval based on the 1st program). (1000 simulation runs minimum)
###Extra Credit: Use apply function to replace the for loop in bootstrap function. Fix the problems that occur same as sdb at N = 3

LogNorm.sim <- function(mu.val = 3, n=30, nsim = 1000, cl = 0.1){
  #Initialize the variables that will iteratively stores the logical value (either 0 or 1) in each simulation that will compare the lognormal mean and coverage of bootstrap
  cvec.boot <- NULL
  cvec.norm <- NULL
  
  #Create a variable that stores the lognormal mean
  mu_lnorm <- exp(mu.val+1/2)
  
  #Do iterations nsim times
  for (i in 1:nsim){
    
    #print i whenever i is a multiple of 10
    #if (i/10 == floor(i/10)) print(i)
    
    #Generate n random sample in lognormal distribution based on 3 as the mean of lognormal distribution
    rand.lnorm.sample <- rlnorm(n,mu.val)
    
    #Use bootstrap function we created above to get the boot confidence interval and t-distribution confidence interval
    boot.list <- ConfItv.boot(b_vec = rand.lnorm.sample, alpha = cl)
    confboot <- boot.list$my.boot_interval
    conft <- boot.list$normal_interval
    
    
    #logical value list
    #For each entry, if each confidence interval covers the mean of lognormal distribution, it is 1
    cvec.boot <- c(cvec.boot, (confboot[1] < mu_lnorm)*(confboot[2] > mu_lnorm))
    cvec.norm <- c(cvec.norm, (conft[1] < mu_lnorm)*(conft[2] > mu_lnorm))
  }
  
  #returns the probability that confidence interval(boot and t) that covers the mean of lognormal distribution
  list(boot.cov = sum(cvec.boot)/nsim, t.cov = sum(cvec.norm)/nsim)
}



#3. 
#Compare the coverage rates for the bootstrap confidence interval and the central limit theorem based confidence interval. For sample sizes 3, 10,30, and 100, alpha= .1 (90% ci), and .05 (95% CI)

#Implement
#Try NOAA and GISS csv data
Output.1 <- ConfItv.boot(NG$X.disaster)
names(Output.1$my.boot_interval)[1] <- "5%"
names(Output.1$my.boot_interval)[2] <- "95%"
Output.1

tb_temp <- matrix(ncol = 2, nrow = 8, 0)
rownames(tb_temp) <- c("n=3,alpha=0.1","n=3,alpha=0.05","n=10,alpha=0.1","n=10,alpha=0.05","n=30,alpha=0.1","n=30,alpha=0.05","n=100,alpha=0.1","n=100,alpha=0.05")
colnames(tb_temp) <- c("boot","normal")
prob_tb <- as.table(tb_temp)
#n = 3,10,30,100, alpha = 0.1,0.05
row_count <- 1
for (i in c(3,10,30,100)){
  for (j in c(0.1,0.05)){
    print(paste("For n = ", toString(i), " alpha = ", toString(j)))
    temp_prob <- LogNorm.sim(n = i, cl = j)
    prob_tb[row_count,1] <- temp_prob$boot.cov
    prob_tb[row_count,2] <- temp_prob$t.cov
    row_count <- row_count + 1
  }
}
prob_tb

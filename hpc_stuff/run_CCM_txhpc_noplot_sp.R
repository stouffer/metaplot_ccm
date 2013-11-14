#source("CCM_multi_Library.R")
#Should take 100hrs for 1000 iterations on 1 procesor
source("CCM_multi_Library_tx.R")

require(doParallel)
require(deSolve)

#Parameters that will change among simulations
cl <- makeCluster(1) #Number of CPU cores to use for simluations
iter_max<-10 #Number of iterations
number_of_chains=1 #Number of plots
filename<-"datout" #Name of output file

#Register cores
registerDoParallel(cl)

#Names for output data
datnames<-c(paste("T_", c("sl_025", "sl_160", "sl_500", "sl_840", "sl_975", "rho", "vrho", "prho", "cause"), sep=""),
            paste("M_", c("sl_025", "sl_160", "sl_500", "sl_840", "sl_975", "rho", "vrho", "prho", "cause"), sep=""),
            "xstr")

#Set up parameters
seednum=1114
times=seq(1, 501, by = 1)
x_mean_sd = c(0, 1)
y_mean_sd = c(0, 1)
Wrates=c(1, 1)
OUrates=c(0.1, 0.1)

#Make system
compsimout<-make_comp_data_boot(seednum=seednum, times=times, x_mean_sd = x_mean_sd, y_mean_sd = y_mean_sd, Wrates=Wrates, OUrates=OUrates)
xstr_list<-seq(0, 1, length=iter_max)

system.time({
  #dat_out<-foreach(z_iter_counter=1:iter_max, .combine=rbind) %dopar% {
  for(z_iter_counter in 1:iter_max) {
    program_init_bootstrap_short()
    #Initialize and make data
    xstr<-xstr_list[z_iter_counter] #Pick interaction strength (0 is 100% y, 1 is 100% x)
    
    ode_result<-FALSE
    ode_result<-make_comp_sim(compsimout, seednum=seednum, xstr=xstr,  times=times, number_of_chains=number_of_chains, x_mean_sd = x_mean_sd, y_mean_sd = y_mean_sd, Wrates=Wrates, OUrates=OUrates)
    #plot_output_boot(ode_result)
    
    #fix xlist and ylist
    times<-ode_result$pars$times
    if(times[1]!=0) {
      times<-c(0, ode_result$pars$times)
    }
    cutlist<-cut(1:max(times), times)
    ode_result$pars$xlist<-tapply(ode_result$pars$xlist, cutlist, mean)
    ode_result$pars$ylist<-tapply(ode_result$pars$ylist, cutlist, mean)
    
    ccm_out<-doCCM_environment(ode_result=ode_result, target_sp="all", predstep=1, tau=1, maxE=7, iterations=100, twoway=FALSE)
    
    test_out<-testccm_boot(ccm_out, iter=1000)
    c(test_out$T_test, test_out$M_test, xstr)
  }
  colnames(dat_out)<-datnames
  write.csv(dat_out, paste(filename, ".csv", sep=""), row.names=FALSE)
})
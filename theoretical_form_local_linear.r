setwd("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/inference/theoretical variance")
library(MASS)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(doSNOW)
###############################################
rho   = 1.0
sigma = 0.030
mu=6
R=0.18
DT=30
m=100
h.max = 0.016
h=0.010
#
# pcf grid
n.grid= 500
grid=seq(0,R,length.out = n.grid)
gridt=R/n.grid
######   apprximate the Q1&Q2 ###############
# calculus grid
n_grid_c = 8000
bin_grid_c   = DT/n_grid_c 
grid_c= seq(bin_grid_c,DT,by = bin_grid_c)
calculus_bin_area = bin_grid_c^2

lam_fun = function(t){0.28*rho*mu*(sin(2*pi*t)+sin(4*pi*t)+1.811256)}
#
grid_uv=expand.grid(grid_c,grid_c)
grid_uv = grid_uv %>% mutate(dis=abs(Var1-Var2))
grid_uv_mat =grid_uv[grid_uv$dis <=R+h.max,]

cl = makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

ptm<-proc.time()
Q_mat_cal = foreach(i = 1:50,.packages = c("dplyr")) %dopar% {
  r=grid[i] # each time point that pcf being estimated 
  grid_uv =grid_uv_mat %>% mutate(w_h = ifelse(abs(dis-r)<=h,1/(2*h*(DT-dis)),0))
  grid_uv = grid_uv[grid_uv$w_h>0,]
  grid_uv = grid_uv %>% mutate(lam_prod = lam_fun(Var1)*lam_fun(Var2),a12=(dis-r)/h)
  grid_uv = grid_uv %>% mutate(Q1_A11=lam_prod*w_h,Q1_A12 = lam_prod*w_h*a12,Q1_A22 = lam_prod*w_h*(a12)^2)
  grid_uv = grid_uv %>% mutate(Q2_A11=(DT*h)*Q1_A11*w_h,Q2_A12 = (DT*h)*Q1_A12*w_h,Q2_A22 = (DT*h)*Q1_A22*w_h)
  Q1_vec =c(calculus_bin_area*sum(grid_uv$Q1_A11),calculus_bin_area*sum(grid_uv$Q1_A12),
            calculus_bin_area*sum(grid_uv$Q1_A12),calculus_bin_area*sum(grid_uv$Q1_A22))
  Q2_vec =c(calculus_bin_area*sum(grid_uv$Q2_A11),calculus_bin_area*sum(grid_uv$Q2_A12),
            calculus_bin_area*sum(grid_uv$Q2_A12),calculus_bin_area*sum(grid_uv$Q2_A22))
  Q1=matrix(Q1_vec,nrow = 2)
  Q2=matrix(Q2_vec,nrow = 2)
  list(Q1,Q2)
}
stopCluster(cl)
Q1_list=Q2_list=NULL
for(i in 1:50){
  Q1_list[[i]]=Q_mat_cal[[i]][[1]]
  Q2_list[[i]]=Q_mat_cal[[i]][[2]]
}

Q1_list2=lapply(Q1_list,function(x) ginv(x, tol = sqrt(.Machine$double.eps)))

# pcf part
pcf=function(t){
  1+(exp(-(t^2)/((2*sigma)^2))/(rho*sqrt(4*pi)*(sigma)))
}
pcfGrid=pcf(grid)
pcf_part =2*pcfGrid*(1+pcfGrid/(m-1))
pcf_part =2*pcfGrid*(1+pcfGrid/(m-1))

variance_vector=rep(0,n.grid)

for(k in 1:50){
  e_Q_part =  Q1_list2[[k]]%*%Q2_list[[k]]%*%Q1_list2[[k]]
  variance_vector[k]=pcf_part[k]*e_Q_part[1,1]
}
variance_vector50=variance_vector[1:50]/(m*DT*h)

#
variance_vector2 =rep(0,n.grid)
for(k in 1:50){
  cal_mat = matrix(calculus_cal[k,],byrow = T,ncol = 2)
  e_Q_part2 =  Q1_list2[[k]]%*% cal_mat %*%Q1_list2[[k]]
  variance_vector2[k]=e_Q_part2[1,1]*6/m
}
variance_vector2[1:50]

variance_vector3 =rep(0,n.grid)
for(k in 1:50){
  cal_mat = matrix(calculus_cal_z1[k,],byrow = T,ncol = 2)
  e_Q_part3 =  Q1_list2[[k]]%*% cal_mat %*%Q1_list2[[k]]
  variance_vector3[k]=e_Q_part3[1,1]*6/m
}
variance_vector3[1:50]

load("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/simulation with selection/setting1/rho2_sigma0.030/rho=1.0/exp local linear/exp local linear estimator result.RData")
estResult=Simu.array[5,,]
var_ll = apply(estResult,1,var)

result_50_ll =cbind(variance_vector50,var_ll[1:50])
save(result_50_ll,file="result_50_ll_jan24.RData")


result_compare =cbind(variance_vector50,variance_vector2[1:50],variance_vector3[1:50],var_ll[1:50])
colnames(result_compare) = c("theoretical_var","third_z2_cov","third_z1","simulation_var")
save(result_compare,file = "result_compare_Feb12.RData")
result_compare2 = data.frame(result_compare) %>% mutate(var_theory=theoretical_var+third_z2_cov+third_z1)
save(result_compare2,file = "result_compare2_Feb13.RData")
library(MASS)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(doSNOW)
###############################################
setwd("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/inference/theoretical variance/variance_Z_one/variance_Z_one")
#
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
#lam_fun = function(t){rho*mu*(t-t+1)}
#
grid_uv=expand.grid(grid_c,grid_c)
grid_uv = grid_uv %>% mutate(dis=abs(Var1-Var2))
grid_uv_mat =grid_uv[grid_uv$dis <=R+h.max,]

cl = makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

ptm<-proc.time()
Q_mat_cal = foreach(i = 1:n.grid,.combine = 'rbind',.packages = c("dplyr")) %dopar% {
  r=grid[i] # each time point that pcf being estimated 
  grid_uv =grid_uv_mat %>% mutate(w_h = ifelse(abs(dis-r)<h,1/(2*h*(DT-dis)),0))
  grid_uv = grid_uv[grid_uv$w_h>0,]
  grid_uv = grid_uv %>% mutate(lam_prod = lam_fun(Var1)*lam_fun(Var1))
  grid_uv = grid_uv %>% mutate(q1=lam_prod*w_h,q2 = (DT*h)*lam_prod*w_h^2)
  Q1 =calculus_bin_area*sum(grid_uv$q1)
  Q2 =calculus_bin_area*sum(grid_uv$q2)
  c(Q1,Q2)
}
stopCluster(cl)
Q1_vec=Q_mat_cal[,1]
Q2_vec=Q_mat_cal[,2]

# pcf part
pcf=function(t){
  1+(exp(-(t^2)/((2*sigma)^2))/(rho*sqrt(4*pi)*(sigma)))
}
pcfGrid=pcf(grid)
pcf_part =2*pcfGrid*(1+pcfGrid/(m-1))
#pcf_part =2*pcfGrid
variance_proof=variance_triple_z1=variance_triple_z2_cov=rep(0,n.grid)

for(k in 1:n.grid){
  e_Q_part =  Q2_vec[k]/((Q1_vec[k])^2)
  variance_proof[k]=pcf_part[k]*e_Q_part
}
variance_proof=variance_proof[1:n.grid]/(m*DT*h)

load("cov_var_calculus_local_constant_feb14.RData")
for(k in 1:n.grid){
  e_Q_part2 =  6*calculus_cal[k]/(m*(Q1_vec[k])^2)
  variance_triple_z2_cov[k]=e_Q_part2
}
variance_triple_z2_cov[1:n.grid]

load("triple_calculus_Z1_local_constant_feb14.RData")
for(k in 1:n.grid){
  e_Q_part3 =  4*calculus_cal[k]/(m*(Q1_vec[k])^2)
  variance_triple_z1[k]=e_Q_part3
}
variance_triple_z1[1:n.grid]
# load estimator result
load("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/simulation with selection/setting1/rho2_sigma0.030/rho=1.0/local constant/local constant estimator result.RData")
estResult=Simu.array[5,,]
var_lc = apply(estResult,1,var)

result= cbind(var_lc[1:n.grid],variance_proof,variance_triple_z1[1:n.grid],variance_triple_z2_cov[1:n.grid])
colnames(result) = c("variance_simu","variance_proof","variance_triple_z1","variance_triple_z2_and_cov")
result_combine_lc =data.frame(result) %>% mutate(variance_proof_and_triple =variance_proof+variance_triple_z1+variance_triple_z2_and_cov )
save(result_combine_lc,file="result_combine_lc_Feb14.RData")
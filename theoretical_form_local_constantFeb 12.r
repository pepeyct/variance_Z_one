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
#lam_fun = function(t){rho*mu*(t-t+1)}
#
grid_uv=expand.grid(grid_c,grid_c)
grid_uv = grid_uv %>% mutate(dis=abs(Var1-Var2))
grid_uv_mat =grid_uv[grid_uv$dis <=R+h.max,]

cl = makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

ptm<-proc.time()
Q_mat_cal = foreach(i = 1:50,.combine = 'rbind',.packages = c("dplyr")) %dopar% {
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
variance_vector=variance_vector2=rep(0,n.grid)

for(k in 1:50){
  e_Q_part =  Q2_vec[k]/((Q1_vec[k])^2)
  variance_vector[k]=pcf_part[k]*e_Q_part
}
variance_vector50=variance_vector[1:50]/(m*DT*h)

for(k in 1:50){
  e_Q_part2 =  2*calculus_cal[k]/(m*(Q1_vec[k])^2)
  variance_vector2[k]=e_Q_part2
}
variance_vector2[1:50]
# load estimator result

estResult=Simu.array[5,,]
var_lc = apply(estResult,1,var)

result50= cbind(variance_vector50,var_lc[1:50])

# 
set.seed(123)
#
library(dplyr)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)
###############################################
rho = 0.8
mu = 3.0
sigma = 0.030
R=0.18
DT=30
m=200
h.max = 0.016
h=0.010
#
# pcf grid
n.grid= 500
grid=seq(0,R,length.out = n.grid)
gridt=R/n.grid
######   apprximate the Q1&Q2 ###############
# calculus grid
n_grid_c = 5000
bin_grid_c   = DT/n_grid_c 
grid_c= seq(bin_grid_c,DT,by = bin_grid_c)
calculus_bin_area = bin_grid_c^3
## theoretical pcf
# pcf part
pcf=function(t){
  1+(exp(-(t^2)/((2*sigma)^2))/(rho*sqrt(4*pi)*(sigma)))
}
pcfGrid=pcf(grid)
### all combinations of three grid vectors that can be used to appximate the calculas
nn1 = ceiling((R+h.max)/bin_grid_c)
grid_list =list(n_grid_c)
for( i in 1:n_grid_c){
  point_u = grid_c[i]
  #
  if(i >=nn1 & i <= (n_grid_c -nn1)){
    vec_v= grid_c[(i-nn1):(i+nn1)]
  }else{ if(i < nn1){
    vec_v = grid_c[1:(i+nn1)]
  } else {
    vec_v = grid_c[(i-nn1):n_grid_c]
  }}
  grid_list[[i]]= t(expand.grid(point_u,vec_v,vec_v))
}
###################################################
grid_list_frame = data.frame(matrix(unlist(grid_list),byrow = T,ncol = 3))

grid_triple_point = grid_list_frame %>% mutate(dis1 = abs(X1-X2),dis2 = abs(X1-X3),dis3 = abs(X2-X3))
# grid_uv=expand.grid(grid_c,grid_c)
# grid_uv = grid_uv %>% mutate(dis=abs(Var1-Var2))
# grid_uv_mat =grid_uv[grid_uv$dis <=R+h.max,]
# grid_double_uv = expand.grid(grid_uv$Var1,grid_c)

#
cl = makeCluster(4, type = "SOCK")
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











# 
set.seed(123)
#
setwd("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/inference/theoretical variance/variance_Z_one/variance_Z_one")
library(dplyr)
library(data.table)
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
######   apprximate the third order calculus in variance z1 ###############
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
nn1 = ceiling((R+h+0.002)/bin_grid_c)
grid_list =list(n_grid_c)
for( i in 1:n_grid_c){
  point_u = grid_c[i]
  #
  if(i >nn1 & i < (n_grid_c -nn1)){
    vec_v= grid_c[(i-nn1):(i+nn1)]
  }else{ if(i <= nn1){
    vec_v = grid_c[1:(i+nn1)]
  } else {
    vec_v = grid_c[(i-nn1):n_grid_c]
  }}
  grid_list[[i]]= t(expand.grid(point_u,vec_v,vec_v))
}
###################################################
grid_list_frame = data.frame(matrix(unlist(grid_list),byrow = T,ncol = 3))
grid_triple_point_fill = grid_list_frame %>% mutate(dis1 = abs(X1-X2),dis2 = abs(X1-X3),dis3 = abs(X2-X3))
rm(grid_list_frame,grid_list)
gc()
# calculate the calculus 
normal_fun = function(x){exp(-x^2/(2*sigma^2))/(sigma*sqrt(2*pi))}
#
cl = makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

ptm<-proc.time()
calculus_cal_z1 = foreach(i = 1:50,.combine = 'rbind',.packages = c("dplyr")) %dopar% {
  r=grid[i] # each time point that pcf being estimated 
  grid_triple_point =grid_triple_point_fill %>% mutate(w_h_1 = ifelse(abs(dis1-r)<h,1/(2*h*(DT-dis1)),0),
                                                  w_h_2= ifelse(abs(dis2-r)<h,1/(2*h*(DT-dis2)),0))
  grid_triple_point = grid_triple_point[grid_triple_point$w_h_1>0 & grid_triple_point$w_h_2>0,]
  # calculate the third order pair correlation fun
  g3=apply(grid_triple_point,1,function(row_Vec){integrate(f=function(p){normal_fun(p+row_Vec[1])*normal_fun(p+row_Vec[2])*normal_fun(p+row_Vec[3])},lower=-Inf,upper=Inf)$value})
  grid_triple_point = grid_triple_point %>% mutate(g3_cls = g3/(rho^2)+pcf(dis1)+pcf(dis2)+pcf(dis3)-2,a12=(dis2-r)/h,a21=(dis1-r)/h,a22=((dis2-r)*(dis1-r))/(h^2))
  grid_triple_point = grid_triple_point %>% mutate(c11=w_h_1*w_h_2*g3_cls,c12 = a12*c11,c21 = a21*c11,c22=a22*c11)
  C_vec= colSums(grid_triple_point[,13:16])*bin_grid_c^3
  C_vec
}
save(calculus_cal_z1,file = "calculus_ca_feb11.RData")
proc.time()-ptm
# three order calculus



